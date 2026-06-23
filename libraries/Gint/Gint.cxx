#include <iostream>

#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TEnv.h>

#include <fstream>
#include <sstream>
#include <vector>
#include <TH1D.h>
#include <TMemFile.h>
#include <TParameter.h>
#include <TObjString.h>

#include <Gint.h>
#include <Gtypes.h>
#include <argParser.h>
#include <Histomatic.h>

#include <globals.h>
#include <GGlobals.h>

Gint *Gint::fGint = 0;

//Gint::Gint(int argc, char **argv) : TRint("gint",&argc,argv,0,0,true,false) {
Gint::Gint(int argc, char **argv) : TRint("gint",0,0,0,0,true,false), 
  fCalC0(0), fCalC1(1), fCalC2(0), fCalUnit("keV"), fHasCalibration(false),
  fRootFilesOpened(0), fTabLock(false), fMainThreadId(std::this_thread::get_id())  {

  LoadOptions(argc,argv);
  LoadStyle();
  SetPrompt("groot [%d] ");

}

Gint *Gint::Get(int argc,char **argv) {
  if(!fGint)
    fGint = new Gint(argc,argv);
  return fGint;
}

Gint::~Gint() { }

void Gint::Terminate(int status) {
  printf("\nbye,bye\n\n");
  SetPrompt("");
  TRint::Terminate(status);

}

void Gint::LoadStyle() {
  // Load the ROOT style file
  //gStyle->SetPalette(kVisibleSpectrum);
  gStyle->SetPalette(gEnv->GetValue("Gint.Style",kVisibleSpectrum));
  gStyle->SetHistLineWidth(2);
  gStyle->SetHistFillStyle(0);
 
  gStyle->SetFrameBorderMode(1);
  //gStyle->SetFrameFillColor(1);
  gStyle->SetCanvasBorderMode(1);
  gStyle->SetPadBorderMode(1);
  gStyle->SetPadBorderSize(5);
  //gStyle->SetPadColor(2);


//canvas.SetBorderSize(6);
//canvas.SetHighLightColor(kBlue); //or whatever color


  //Hists.Stats: "nemri"
  gStyle->SetOptStat(gEnv->GetValue("Hists.Stats","nemri"));

  //gEnv->Print();


  gROOT->ForceStyle();
}


void Gint::LoadOptions(int argc, char **argv) {
  //check the grutrc file for set preset optrions....
  

  argParser parser;

  std::vector<std::string> input_files;
  bool doHelp,doGui,doVersion;

  parser.default_option(&input_files)
    .description("Input file(s)");
  parser.option("h help ?",&doHelp)
    .description("Show this help Message")
    .default_value(false);
  parser.option("g gui",&doGui)
    .description("Start the GUI")
    .default_value(false);
  parser.option("v version",&doVersion)
    .description("Show version")
    .default_value(false);
 
  // Do the parsing...
  try{
    parser.parse(argc, argv);
  } catch (ParseError& e){
    std::cerr << "ERROR: " << e.what() << "\n"
              << parser << std::endl;
    //fShouldExit = true;
  }

 

  // Print help if requested.
  if(doHelp){
    //Version();
    std::cout << parser << std::endl;
    //fShouldExit = true;
  }
  if(doVersion) {
    //Version();
    printf("version not available.\n");
    //fShouldExit = true;
  }
  if(doGui) {
    printf("starting gui...\n"); 
    Histomatic::Get(); 
    //gHistomatic = new Histomatic;
  }

  
  for(auto& file : input_files){
    switch(DetermineFileType(file)){
      case kFileType::CALIBRATION:
	LoadCalibrationFile(file);
        break;
      case kFileType::ROOTFILE:
        {
          TFile *rfile = OpenRootFile(file);
          if(rfile && doGui && gHistomatic) 
            gHistomatic->AddRootFile(rfile);
        }
        break;
      case kFileType::TXT3:
        {
	  TH1D *hist = OpenTxt3File(file);
	  if(hist && doGui && gHistomatic) {
     		std::string name = file.substr(file.find_last_of("/\\") + 1);
     		name += ".root";

     		TMemFile *memfile = new TMemFile(name.c_str(), "RECREATE");
     		memfile->cd();
     		hist->Write();

      		gHistomatic->AddRootFile(memfile);
	    }
	}
        break;
      case kFileType::MACRO:
        break;
      case kFileType::CUTS:
        break;
      default:
        printf("\tDiscarding unknown file: %s\n",file.c_str());
      break;
    };
  }



}

kFileType Gint::DetermineFileType(const std::string& filename) const {
  size_t dot = filename.find_last_of('.');
  if(dot==std::string::npos) return kFileType::UNKNOWN;
  std::string ext = filename.substr(dot+1);

  if((ext=="gz") || (ext=="bz2") || (ext=="zip")) {
    std::string remaining = filename.substr(0,dot);
    ext = remaining.substr(remaining.find_last_of('.')+1);
  }
 
  if(ext == "cal") {
    return kFileType::CALIBRATION;
  } else if(ext == "root") {
    return kFileType::ROOTFILE;
  } else if((ext=="c") || (ext=="C") 
            || (ext=="c+") || (ext=="C+") 
            || (ext=="c++") || (ext=="C++")) {
    return kFileType::MACRO;
  } else if(ext == "cuts") {
    return kFileType::CUTS;
  } else if(ext == "txt3") {
    return kFileType::TXT3;
  } else {
    return kFileType::UNKNOWN;
  }
};

/*
bool Gint::FileAutoDetect(const std::string& filename) {
  switch(DetermineFileType(filename)){
    case kFileType::CALIBRATION:
      break;
    case kFileType::ROOTFILE:
      OpenRootFile(filename);
      break;
    case kFileType::MACRO:
      break;
    case kFileType::CUTS:
    default:
      printf("\tDiscarding unknown file: %s\n",filename.c_str());
      return false;    
    break;
  };
  return true;
}
*/

bool Gint::LoadCalibrationFile(const std::string& filename) {
  std::ifstream input(filename);
  if(!input.is_open()) {
    std::cout << "Could not open calibration file: " << filename << std::endl;
    return false;
  }




  std::string line;
  while(std::getline(input, line)) {
    if(line.empty()) {
      continue;
    }

    if(line.find("C0") != std::string::npos) {
      sscanf(line.c_str(), "C0 = %lf", &fCalC0);
    } else if(line.find("C1") != std::string::npos) {
      sscanf(line.c_str(), "C1 = %lf", &fCalC1);
    } else if(line.find("C2") != std::string::npos) {
      sscanf(line.c_str(), "C2 = %lf", &fCalC2);
    } else if(line.find("unit") != std::string::npos) {
      size_t unitPos = line.find("=");
      if(unitPos != std::string::npos) {
        fCalUnit = line.substr(unitPos + 1);
        while(!fCalUnit.empty() && fCalUnit[0] == ' ') {
          fCalUnit.erase(0, 1);
        }
      }
    }
  }

  fHasCalibration = true;

  std::cout << "\tloaded calibration: "
            << "C0=" << fCalC0
            << ", C1=" << fCalC1
            << ", C2=" << fCalC2
            << ", unit=" << fCalUnit
            << std::endl;

  return true;
}

auto parseTimeSeconds = [](const std::string& line) {
  size_t equals = line.find('=');

  if(equals == std::string::npos)
    return 0.0;

  std::string value = line.substr(equals + 1);

  int hours = 0;
  int minutes = 0;
  double seconds = 0;

  if(sscanf(value.c_str(),
            " %d:%d:%lf",
            &hours,
            &minutes,
            &seconds) == 3) {
    return hours * 3600.0 +
           minutes * 60.0 +
           seconds;
  }

  try {
    return std::stod(value);
  } catch(...) {
    return 0.0;
  }
};




TH1D* Gint::OpenTxt3File(const std::string& filename) {
  std::ifstream input(filename);
  if(!input.is_open()) {
    std::cout << "Could not open txt3 file: " << filename << std::endl;
    return nullptr;
  }

  std::vector<double> channels;
  std::vector<double> counts;
  std::vector<double> energies;

  double c0 = 0;
  double c1 = 1;
  double c2 = 0;
  std::string unit = "keV";
  double liveTime = 0;
  double realTime = 0;




  std::string line;
  while(std::getline(input, line)) {
    if(line.empty()) {
      continue;
    }

    if(line.find("C0") != std::string::npos) {
      sscanf(line.c_str(), "C0 = %lf; C1 = %lf; C2 = %lf;", &c0, &c1, &c2);
      
	size_t unitPos = line.find("unit =");
	if(unitPos != std::string::npos) {
	  unit = line.substr(unitPos + 6);
	  while(!unit.empty() && unit[0] == ' ') {
	  unit.erase(0,1);
	}
	}
	continue;
    }

if(line.find("RealTime") != std::string::npos) {
  realTime = parseTimeSeconds(line);
  continue;
}

if(line.find("LiveTime") != std::string::npos) {
  liveTime = parseTimeSeconds(line);
  continue;
}  

    std::stringstream ss(line);
    double channel = 0;
    double count = 0;
    double energy = 0;

    if(ss >> channel >> count >> energy) {
      channels.push_back(channel);
      counts.push_back(count);
      energies.push_back(energy);
    }
  }

  if(counts.empty()) {
    std::cout << "No spectrum data found in txt3 file: " << filename << std::endl;
    return nullptr;
  }

  //int nBins = counts.size();
  //double xMin = energies.front() - 0.5;
  //double xMax = energies.back() + 0.5;

  //TH1D* hist = new TH1D("txt3_spectrum", filename.c_str(), nBins, xMin, xMax);

  if(fHasCalibration) {
	c0 = fCalC0;
	c1 = fCalC1;
	c2 = fCalC2;
	unit = fCalUnit;
  }

  int nBins = counts.size();

  std::vector<double> binEdges;
  binEdges.reserve(nBins + 1);

  for(int i = 0; i <= nBins; ++i) { 
    double channelEdge = channels.front() - 0.5 + i;
    double energyEdge = c0 + c1 * channelEdge + c2 * channelEdge * channelEdge;
    binEdges.push_back(energyEdge);
  }

  TH1D* hist = new TH1D("TXT3 -_- Spectrum", filename.c_str(), nBins, binEdges.data());
  
  hist->GetListOfFunctions()->Add(new TParameter<double>("C0", c0));
  hist->GetListOfFunctions()->Add(new TParameter<double>("C1", c1));
  hist->GetListOfFunctions()->Add(new TParameter<double>("C2", c2));
  hist->GetListOfFunctions()->Add(new TParameter<double>("LiveTime", liveTime));
  hist->GetListOfFunctions()->Add(new TParameter<double>("RealTime", realTime));
  hist->GetListOfFunctions()->Add(new TObjString(("unit=" + unit).c_str()));

  hist->GetXaxis()->SetTitle(unit.c_str());
  hist->GetYaxis()->SetTitle("Counts");

  for(int i = 0; i < nBins; ++i) {
    hist->SetBinContent(i + 1, counts[i]);
  }

  std::cout << "\topened txt3 spectrum: " << filename << std::endl;
  return hist;
}







TFile *Gint::OpenRootFile(const std::string& filename, Option_t* opt) {
  TString sopt(opt);
  sopt.ToLower();

  TFile *file = NULL;
  if(sopt.Contains("recreate") || sopt.Contains("new")) {
    file = new TFile(filename.c_str(),"recreate");
    if(!file->IsOpen()) { delete file; file = NULL; }
    if(file) {
      const char* command = Form("TFile* _file%i = (TFile*)%luL",
                                 fRootFilesOpened,
                                 (unsigned long)file);
      TRint::ProcessLine(command);
      fRootFilesOpened++;
    } else {
      std::cout << "Could not create " << filename << std::endl;
    }
  } else {
    //file = TFile::Open(filename.c_str(),opt);
    file = new TFile(filename.c_str(),opt);
    if(!file->IsOpen()) { delete file; file = NULL; }
    if(file) {
      const char* command = Form("TFile* _file%i = (TFile*)%luL",
                                 fRootFilesOpened,
                                 (unsigned long)file);
      TRint::ProcessLine(command);
      std::cout << "\tfile " << BLUE << file->GetName() << RESET_COLOR
                <<  " opened as " << BLUE <<  "_file" << fRootFilesOpened << RESET_COLOR <<  std::endl;

      fRootFilesOpened++;
    } else {
      std::cout << "Could not create " << filename << std::endl;
    }
  }

  return file;
}

int Gint::TabCompletionHook(char* buf, int* pLoc, std::ostream& out) {
  fTabLock = true;
  int result = TRint::TabCompletionHook(buf,pLoc,out);
  fTabLock = false;
  return result;
}


long Gint::ProcessLine(const char* line, bool sync, int* error) {
  long retval = 0;
  if(fTabLock) {
    return TRint::ProcessLine(line,sync,error);
  }
  TString sline(line);
  if(!sline.Length()) {
    return 0;
  }
  sline.ReplaceAll("TCanvas","GCanvas");

  if(std::this_thread::get_id() != fMainThreadId){
    printf("Not the main thread...  ");
    fflush(stdout);
  }           
  
  if(!sline.CompareTo("clear")) {
    retval = TRint::ProcessLine(".! clear");
  } else {
    retval = TRint::ProcessLine(sline.Data(),sync,error);
  }

  if(retval < 0) {
    //std::cerr << "Error processing line: " << line << std::endl;
  }
  return retval;
}



