#include <GCalibration.h>

#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>

#include <fstream>
#include <iostream>
#include <sstream>

ClassImp(GCalibration)

GCalibration::GCalibration()
  : TObject(),
    fOrder(1),
    fC0(0),
    fC1(1),
    fC2(0),
    fUnit("keV") {
}

GCalibration::~GCalibration() {
}

void GCalibration::Clear(Option_t* opt) {
  TObject::Clear(opt);
  fPoints.clear();
  fResiduals.clear();
  fOrder = 1;
  fC0 = 0;
  fC1 = 1;
  fC2 = 0;
  fUnit = "keV";
}

void GCalibration::AddPoint(double channel,
                            double energy,
                            double channelErr,
                            double energyErr,
                            bool enabled) {
  GCalibrationPoint point;
  point.channel = channel;
  point.energy = energy;
  point.channelErr = channelErr;
  point.energyErr = energyErr;
  point.enabled = enabled;

  fPoints.push_back(point);
}

bool GCalibration::LoadPoints(const char* pointsFile) {
  std::ifstream input(pointsFile);
  if(!input.is_open()) {
    std::cout << "Could not open calibration points file: "
              << pointsFile << std::endl;
    return false;
  }

  Clear();

  std::string line;
  while(std::getline(input, line)) {
    if(line.empty() || line[0] == '#')
      continue;

    std::stringstream ss(line);

    double channel = 0;
    double energy = 0;
    double channelErr = 0;
    double energyErr = 0;

    if(!(ss >> channel >> energy))
      continue;

    ss >> channelErr >> energyErr;
    AddPoint(channel, energy, channelErr, energyErr);
  }

  if(fPoints.size() < 2) {
    std::cout << "Need at least two calibration points." << std::endl;
    return false;
  }

  return true;
}

bool GCalibration::Fit(int order, const char* unit) {
  if(order < 1) order = 1;
  if(order > 2) order = 2;

  fOrder = order;
  fUnit = unit ? unit : "keV";

  std::vector<GCalibrationPoint> enabled;
  for(const auto& point : fPoints) {
    if(point.enabled)
      enabled.push_back(point);
  }

  const size_t minPoints = order == 2 ? 3 : 2;
  if(enabled.size() < minPoints) {
    std::cout << "Need at least " << minPoints
              << " enabled calibration points for order "
              << order << "." << std::endl;
    return false;
  }

  TGraphErrors graph(enabled.size());
  graph.SetName("calibration_graph");
  graph.SetTitle("Calibration;Channel;Energy");

  for(size_t i = 0; i < enabled.size(); ++i) {
    graph.SetPoint(i, enabled[i].channel, enabled[i].energy);
    graph.SetPointError(i, enabled[i].channelErr, enabled[i].energyErr);
  }

  TF1 fit("calibration_fit", order == 1 ? "pol1" : "pol2");
  TFitResultPtr result = graph.Fit(&fit, "QSEX0");

  if(!result.Get() || !result->IsValid() || int(result) != 0) {
    std::cout << "Energy calibration fit failed."
              << " Status = " << int(result)
              << std::endl;
    return false;
  }

  fC0 = fit.GetParameter(0);
  fC1 = fit.GetParameter(1);
  fC2 = order == 2 ? fit.GetParameter(2) : 0;

  fResiduals.clear();
  fResiduals.reserve(fPoints.size());

  for(const auto& point : fPoints) {
    if(point.enabled)
      fResiduals.push_back(Eval(point.channel) - point.energy);
    else
      fResiduals.push_back(0);
  }

  return true;
}

double GCalibration::Eval(double channel) const {
  return fC0 + fC1 * channel + fC2 * channel * channel;
}

double GCalibration::Residual(size_t index) const {
  if(index >= fResiduals.size())
    return 0;
  return fResiduals[index];
}


bool GCalibration::Load(const char* calFile) {
  std::ifstream input(calFile);
  if(!input.is_open()) {
    std::cout << "Could not open calibration file: "
              << calFile << std::endl;
    return false;
  }

  Clear();

  std::string line;
  while(std::getline(input, line)) {
    if(line.empty() || line[0] == '#')
      continue;

    if(line.find("C0") != std::string::npos) {
      sscanf(line.c_str(), "C0 = %lf", &fC0);
    } else if(line.find("C1") != std::string::npos) {
      sscanf(line.c_str(), "C1 = %lf", &fC1);
    } else if(line.find("C2") != std::string::npos) {
      sscanf(line.c_str(), "C2 = %lf", &fC2);
    } else if(line.find("unit") != std::string::npos) {
      size_t unitPos = line.find("=");
      if(unitPos != std::string::npos) {
        fUnit = line.substr(unitPos + 1);
        while(!fUnit.empty() && fUnit[0] == ' ')
          fUnit.erase(0, 1);
      }
    }
  }

  fOrder = fC2 == 0 ? 1 : 2;
  return true;
}


bool GCalibration::Save(const char* calFile) const {
  std::ofstream output(calFile);
  if(!output.is_open()) {
    std::cout << "Could not write calibration file: "
              << calFile << std::endl;
    return false;
  }

  output << "C0 = " << fC0 << "\n";
  output << "C1 = " << fC1 << "\n";
  output << "C2 = " << fC2 << "\n";
  output << "unit = " << fUnit << "\n";

  return true;
}
