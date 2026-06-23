#include <GLearning.h>

#include <GCommands.h>
#include <GF1.h>
#include <GFunctions.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCollection.h>
#include <TList.h>
#include <TMath.h>
#include <TParameter.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

ClassImp(GLearning)

GLearning::GLearning(TH1* spectrum)
  : TObject(),
    fSpectrum(nullptr),
    fReferenceActivity(0),
    fHalfLifeSeconds(0),
    fElapsedSeconds(0),
    fLiveTime(0),
    fEfficiencyGraph(nullptr),
    fEfficiencyFit(nullptr) {
  SetSpectrum(spectrum);
}

GLearning::~GLearning() {
  ClearGraphs();
}

void GLearning::ClearGraphs() {
  delete fEfficiencyFit;
  fEfficiencyFit = nullptr;

  delete fEfficiencyGraph;
  fEfficiencyGraph = nullptr;
}

void GLearning::Clear(Option_t* opt) {
  TObject::Clear(opt);

  ClearGraphs();

  fSpectrum = nullptr;
  fSource.clear();
  fLines.clear();
  fPeaks.clear();
  fCalibration.Clear();

  fReferenceActivity = 0;
  fHalfLifeSeconds = 0;
  fElapsedSeconds = 0;
  fLiveTime = 0;
}

void GLearning::SetSpectrum(TH1* spectrum) {
  fSpectrum = spectrum;

  if(!fSpectrum || !fSpectrum->GetListOfFunctions())
    return;

  TObject* object =
    fSpectrum->GetListOfFunctions()->FindObject("LiveTime");

  TParameter<double>* liveTime =
    dynamic_cast<TParameter<double>*>(object);

  if(liveTime)
    fLiveTime = liveTime->GetVal();
}

void GLearning::SetHalfLifeYears(double years) {
  fHalfLifeSeconds = years * 365.25 * 24.0 * 60.0 * 60.0;
}

void GLearning::SetElapsedDays(double days) {
  fElapsedSeconds = days * 24.0 * 60.0 * 60.0;
}

bool GLearning::SetSource(const char* source) {
  if(!source)
    return false;

  fSource = source;
  fLines.clear();
  fPeaks.clear();
  fCalibration.Clear();

  std::string name = source;

  std::transform(
    name.begin(),
    name.end(),
    name.begin(),
    [](unsigned char c) { return std::tolower(c); }
  );

  if(name == "eu-152" || name == "eu152") {
    fSource = "Eu-152";
    SetHalfLifeYears(13.5376);

    // Values copied from Groot's current eu152.sou.
    fLines = {
      {121.782,  0.001, 28.579488 / 100.0, 0.12 / 100.0, true},
      {244.697,  0.001,  7.583340 / 100.0, 0.03 / 100.0, true},
      {344.279,  0.001, 26.522720 / 100.0, 1.50 / 100.0, true},
      {443.976,  0.005,  3.156705 / 100.0, 0.027 / 100.0, true},
      {778.904,  0.002, 12.940970 / 100.0, 0.09 / 100.0, true},
      {867.373,  0.003,  4.245233 / 100.0, 0.029 / 100.0, true},
      {964.079,  0.018, 14.606016 / 100.0, 0.05 / 100.0, true},
      {1085.869, 0.024, 10.453510 / 100.0, 0.04 / 100.0, true},
      {1112.069, 0.003, 13.642824 / 100.0, 0.05 / 100.0, true},
      {1408.006, 0.003, 21.003336 / 100.0, 0.006 / 100.0, true}
    };

    return true;
  }

  if(name == "co-60" || name == "co60") {
    fSource = "Co-60";
    SetHalfLifeYears(5.2714);

    // Values copied from Groot's current co60.sou.
    fLines = {
      {1173.228, 0.003, 99.8500 / 100.0, 0.0300 / 100.0, true},
      {1332.492, 0.004, 99.9826 / 100.0, 0.0006 / 100.0, true}
    };

    return true;
  }

  std::cout << "Unknown source: " << source << "\n"
            << "Supported sources are Eu-152 and Co-60."
            << std::endl;

  fSource.clear();
  return false;
}

const GLearningLine* GLearning::FindLine(double energy) const {
  const GLearningLine* closest = nullptr;
  double closestDistance = 1.0;

  for(const auto& line : fLines) {
    double distance = std::abs(line.energy - energy);

    if(distance < closestDistance) {
      closest = &line;
      closestDistance = distance;
    }
  }

  return closest;
}

double GLearning::GetCurrentActivity() const {
  if(fReferenceActivity <= 0 || fHalfLifeSeconds <= 0)
    return 0;

  double decayConstant = TMath::Log(2.0) / fHalfLifeSeconds;

  return fReferenceActivity *
         std::exp(-decayConstant * fElapsedSeconds);
}

void GLearning::Start() const {
  std::cout << "\nGLearning\n"
            << "Source: " << fSource << "\n"
            << "Live time: " << fLiveTime << " s\n"
            << "\n"
            << "Place two markers around a peak and press g.\n"
            << "Then accept the fit with:\n"
            << "  learn->AcceptLastFit(knownEnergy);\n"
            << std::endl;

  PrintSourceLines();
}

void GLearning::PrintSourceLines() const {
  std::cout << "\nKnown " << fSource << " lines:\n";

  for(const auto& line : fLines) {
    std::cout << "  "
              << line.energy << " keV"
              << "  intensity=" << line.intensity * 100.0 << "%"
              << std::endl;
  }
}

bool GLearning::AcceptLastFit(double knownEnergy) {
  if(!fSpectrum || !fSpectrum->GetListOfFunctions()) {
    std::cout << "No spectrum is attached." << std::endl;
    return false;
  }

  GF1* lastFit = nullptr;

  TIter next(fSpectrum->GetListOfFunctions());

  while(TObject* object = next()) {
    GF1* fit = dynamic_cast<GF1*>(object);

    if(fit)
      lastFit = fit;
  }

  if(!lastFit) {
    std::cout << "No Gaussian fit was found." << std::endl;
    return false;
  }

  return AcceptFit(lastFit, knownEnergy);
}

bool GLearning::AcceptFit(GF1* fit, double knownEnergy) {
  if(!fit)
    return false;

  const GLearningLine* line = FindLine(knownEnergy);

  if(!line) {
    std::cout << "No source line near "
              << knownEnergy << " keV."
              << std::endl;
    return false;
  }

  int centroidParameter = fit->GetParNumber("centroid");

  if(centroidParameter < 0) {
    std::cout << "Fit has no centroid parameter." << std::endl;
    return false;
  }

  GLearningPeak point;

  point.channel = fit->GetParameter(centroidParameter);
  point.channelError = fit->GetParError(centroidParameter);

  point.energy = line->energy;
  point.energyError = line->energyError;

  point.intensity = line->intensity;
  point.intensityError = line->intensityError;

  point.area = fit->GetArea();
  point.areaError = fit->GetAreaErr();

  // Current GF1 area errors may be zero. Use the summed-count
  // uncertainty as a simple first-version fallback.
  if(point.areaError <= 0)
    point.areaError = fit->GetSumErr();

  point.chi2Ndf =
    fit->GetNdf() > 0 ? fit->GetChi2() / fit->GetNdf() : 0;

  fPeaks.push_back(point);

  std::cout << "\nAccepted " << point.energy << " keV\n"
            << "  centroid: " << point.channel
            << " +/- " << point.channelError << "\n"
            << "  area:     " << point.area
            << " +/- " << point.areaError << "\n"
            << "  chi2/NDF: " << point.chi2Ndf
            << std::endl;

  return true;
}

void GLearning::PrintPeaks() const {
  std::cout << "\nAccepted peaks:\n";

  for(size_t i = 0; i < fPeaks.size(); ++i) {
    const auto& point = fPeaks[i];

    std::cout << i
              << "  E=" << point.energy
              << "  channel=" << point.channel
              << "  area=" << point.area
              << "  efficiency=" << point.efficiency
              << "  enabled=" << point.enabled
              << std::endl;
  }
}

bool GLearning::FitEnergyCalibration(int order) {
  fCalibration.Clear();

  for(const auto& point : fPeaks) {
    if(!point.enabled)
      continue;

    fCalibration.AddPoint(
      point.channel,
      point.energy,
      point.channelError,
      point.energyError
    );
  }

  if(fSource == "Co-60")
    order = 1;

  bool success = fCalibration.Fit(order, "keV");

  if(success) {
    std::cout << "\nEnergy calibration:\n"
              << "  C0 = " << fCalibration.GetC0() << "\n"
              << "  C1 = " << fCalibration.GetC1() << "\n"
              << "  C2 = " << fCalibration.GetC2()
              << std::endl;
  }

  return success;
}

TGraphErrors* GLearning::DrawEnergyCalibration() {
  auto* graph = new TGraphErrors();

  int index = 0;

  for(const auto& point : fPeaks) {
    if(!point.enabled)
      continue;

    graph->SetPoint(index, point.channel, point.energy);
    graph->SetPointError(
      index,
      point.channelError,
      point.energyError
    );

    ++index;
  }

  graph->SetName("glearning_energy_calibration");
  graph->SetTitle(
    "Energy Calibration;Centroid Channel;Energy (keV)"
  );

  graph->SetMarkerStyle(20);

  auto* canvas =
    new TCanvas("glearning_calibration_canvas",
                "GLearning Energy Calibration",
                900,
                600);

  graph->Draw("AP");

  const char* formula =
    fCalibration.GetOrder() == 2 ? "pol2" : "pol1";

  auto* fit = new TF1(
    "glearning_calibration_fit",
    formula,
    graph->GetXaxis()->GetXmin(),
    graph->GetXaxis()->GetXmax()
  );

  fit->SetParameter(0, fCalibration.GetC0());
  fit->SetParameter(1, fCalibration.GetC1());

  if(fCalibration.GetOrder() == 2)
    fit->SetParameter(2, fCalibration.GetC2());

  fit->Draw("same");
  canvas->Update();

  return graph;
}

TH1D* GLearning::MakeCalibratedSpectrum(const char* name) {
  if(!fSpectrum)
    return nullptr;

  if(fCalibration.GetPoints().size() < 2) {
    std::cout << "No valid energy calibration exists."
              << std::endl;
    return nullptr;
  }

  int bins = fSpectrum->GetNbinsX();
  std::vector<double> edges(bins + 1);

  for(int i = 0; i < bins; ++i) {
    double channelEdge =
      fSpectrum->GetXaxis()->GetBinLowEdge(i + 1);

    edges[i] = fCalibration.Eval(channelEdge);
  }

  edges[bins] = fCalibration.Eval(
    fSpectrum->GetXaxis()->GetBinUpEdge(bins)
  );

  for(int i = 1; i <= bins; ++i) {
    if(!std::isfinite(edges[i]) ||
       edges[i] <= edges[i - 1]) {
      std::cout << "Invalid calibration: energy edges "
                   "are not increasing."
                << std::endl;
      return nullptr;
    }
  }

  auto* calibrated = new TH1D(
    name,
    fSpectrum->GetTitle(),
    bins,
    edges.data()
  );

  calibrated->SetDirectory(nullptr);

  for(int bin = 1; bin <= bins; ++bin) {
    calibrated->SetBinContent(
      bin,
      fSpectrum->GetBinContent(bin)
    );

    calibrated->SetBinError(
      bin,
      fSpectrum->GetBinError(bin)
    );
  }

  calibrated->GetXaxis()->SetTitle("Energy (keV)");
  calibrated->GetYaxis()->SetTitle("Counts");

  return calibrated;
}

bool GLearning::BuildEfficiencyPoints() {
  double activity = GetCurrentActivity();

  if(activity <= 0) {
    std::cout << "Set the source activity first." << std::endl;
    return false;
  }

  if(fLiveTime <= 0) {
    std::cout << "Set the live time first." << std::endl;
    return false;
  }

  for(auto& point : fPeaks) {
    double emittedRate = activity * point.intensity;
    double denominator = fLiveTime * emittedRate;

    if(denominator <= 0) {
      point.efficiency = 0;
      point.efficiencyError = 0;
      continue;
    }

    point.efficiency = point.area / denominator;

    double relativeAreaError =
      point.area > 0 ? point.areaError / point.area : 0;

    double relativeIntensityError =
      point.intensity > 0
        ? point.intensityError / point.intensity
        : 0;

    point.efficiencyError =
      point.efficiency *
      std::sqrt(
        relativeAreaError * relativeAreaError +
        relativeIntensityError * relativeIntensityError
      );
  }

  PrintPeaks();
  return true;
}

TGraphErrors* GLearning::DrawEfficiency() {
  delete fEfficiencyGraph;

  fEfficiencyGraph = new TGraphErrors();

  int index = 0;

  for(const auto& point : fPeaks) {
    if(!point.enabled || point.efficiency <= 0)
      continue;

    fEfficiencyGraph->SetPoint(
      index,
      point.energy,
      point.efficiency
    );

    fEfficiencyGraph->SetPointError(
      index,
      point.energyError,
      point.efficiencyError
    );

    ++index;
  }

  fEfficiencyGraph->SetName("glearning_efficiency");
  fEfficiencyGraph->SetTitle(
    "Absolute Efficiency;Energy (keV);Absolute Efficiency"
  );

  fEfficiencyGraph->SetMarkerStyle(20);

  auto* canvas =
    new TCanvas("glearning_efficiency_canvas",
                "GLearning Absolute Efficiency",
                900,
                600);

  fEfficiencyGraph->Draw("AP");
  canvas->Update();

  return fEfficiencyGraph;
}

TF1* GLearning::FitEfficiency() {
  if(!fEfficiencyGraph)
    DrawEfficiency();

  if(!fEfficiencyGraph || fEfficiencyGraph->GetN() < 4) {
    std::cout << "Need at least four efficiency points."
              << std::endl;
    return nullptr;
  }

  double minimum = fEfficiencyGraph->GetX()[0];
  double maximum = fEfficiencyGraph->GetX()[0];

  for(int i = 1; i < fEfficiencyGraph->GetN(); ++i) {
    minimum = std::min(minimum, fEfficiencyGraph->GetX()[i]);
    maximum = std::max(maximum, fEfficiencyGraph->GetX()[i]);
  }

  delete fEfficiencyFit;

  fEfficiencyFit = new TF1(
    "glearning_efficiency_fit",
    GFunctions::Efficiency,
    minimum,
    maximum,
    4
  );

  fEfficiencyFit->SetParNames("p0", "p1", "p2", "p3");
  fEfficiencyFit->SetParameters(0, -1, 0, 0);

  TFitResultPtr result =
    fEfficiencyGraph->Fit(fEfficiencyFit, "SEX0R");

  if(!result.Get() || !result->IsValid() || int(result) != 0) {
    std::cout << "Efficiency fit failed."
              << " Status = " << int(result)
              << std::endl;

    delete fEfficiencyFit;
    fEfficiencyFit = nullptr;
    return nullptr;
  }

  return fEfficiencyFit;
}

TF1* GLearning::FitEfficiencyPolynomial(int order) {
  if(order < 2)
    order = 2;

  if(order > 5)
    order = 5;

  if(!fEfficiencyGraph)
    DrawEfficiency();

  if(!fEfficiencyGraph ||
     fEfficiencyGraph->GetN() < order + 1) {
    std::cout << "Not enough points for polynomial order "
              << order << "."
              << std::endl;
    return nullptr;
  }

  double minimum = fEfficiencyGraph->GetX()[0];
  double maximum = fEfficiencyGraph->GetX()[0];

  for(int i = 1; i < fEfficiencyGraph->GetN(); ++i) {
    minimum = std::min(minimum, fEfficiencyGraph->GetX()[i]);
    maximum = std::max(maximum, fEfficiencyGraph->GetX()[i]);
  }

  delete fEfficiencyFit;

  fEfficiencyFit = new TF1(
    "glearning_efficiency_polynomial",
    Form("pol%d", order),
    minimum,
    maximum
  );

  TFitResultPtr result =
    fEfficiencyGraph->Fit(fEfficiencyFit, "SEX0R");

  if(!result.Get() || !result->IsValid() || int(result) != 0) {
    std::cout << "Efficiency polynomial fit failed."
              << " Status = " << int(result)
              << std::endl;

    delete fEfficiencyFit;
    fEfficiencyFit = nullptr;
    return nullptr;
  }

  return fEfficiencyFit;
}

bool GLearning::EnablePoint(size_t index, bool enabled) {
  if(index >= fPeaks.size())
    return false;

  fPeaks[index].enabled = enabled;
  return true;
}

bool GLearning::Save(const char* filename) const {
  std::ofstream output(filename);

  if(!output.is_open())
    return false;

  output << "# GLearning\n";
  output << "source " << fSource << "\n";
  output << "activity " << fReferenceActivity << "\n";
  output << "half_life_seconds " << fHalfLifeSeconds << "\n";
  output << "elapsed_seconds " << fElapsedSeconds << "\n";
  output << "live_time " << fLiveTime << "\n";

  output << "calibration "
         << fCalibration.GetC0() << " "
         << fCalibration.GetC1() << " "
         << fCalibration.GetC2() << "\n";

  output << "# energy channel channel_error area area_error "
            "intensity efficiency efficiency_error enabled\n";

  for(const auto& point : fPeaks) {
    output << "peak "
           << point.energy << " "
           << point.channel << " "
           << point.channelError << " "
           << point.area << " "
           << point.areaError << " "
           << point.intensity << " "
           << point.efficiency << " "
           << point.efficiencyError << " "
           << point.enabled << "\n";
  }

  return true;
}
