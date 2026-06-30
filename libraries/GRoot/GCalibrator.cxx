#include <GCalibrator.h>

#include <GCanvas.h>
#include <GCommands.h>
#include <GDecayNucleus.h>
#include <GFunctions.h>
#include <GGaus.h>
#include <GTransition.h>

#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TSpectrum.h>
#include <TVirtualPad.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

ClassImp(GCalibrator)

namespace {
int NextFitId() {
  static int id = 0;
  return id++;
}

void BuildCombinations(const std::vector<double>& values,
                       size_t count,
                       size_t start,
                       std::vector<double>& current,
                       std::vector<std::vector<double>>& output) {
  if(current.size() == count) {
    output.push_back(current);
    return;
  }

  const size_t remaining = count - current.size();
  for(size_t i = start; i + remaining <= values.size(); ++i) {
    current.push_back(values[i]);
    BuildCombinations(values, count, i + 1, current, output);
    current.pop_back();
  }
}

std::vector<std::vector<double>> Combinations(const std::vector<double>& values, size_t count) {
  std::vector<std::vector<double>> output;
  std::vector<double> current;
  if(count <= values.size())
    BuildCombinations(values, count, 0, current, output);
  return output;
}

bool LinearFit(const std::vector<double>& x,
               const std::vector<double>& y,
               double& intercept,
               double& slope,
               double& chi2) {
  if(x.size() != y.size() || x.size() < 2)
    return false;

  double sumX = 0.0;
  double sumY = 0.0;
  double sumXX = 0.0;
  double sumXY = 0.0;
  for(size_t i = 0; i < x.size(); ++i) {
    sumX += x[i];
    sumY += y[i];
    sumXX += x[i] * x[i];
    sumXY += x[i] * y[i];
  }

  const double n = static_cast<double>(x.size());
  const double denominator = n * sumXX - sumX * sumX;
  if(std::abs(denominator) < std::numeric_limits<double>::epsilon())
    return false;

  slope = (n * sumXY - sumX * sumY) / denominator;
  intercept = (sumY - slope * sumX) / n;

  chi2 = 0.0;
  for(size_t i = 0; i < x.size(); ++i) {
    const double residual = y[i] - (intercept + slope * x[i]);
    chi2 += residual * residual;
  }
  return true;
}

bool HasStableIntercept(const std::vector<double>& peaks,
                        const std::vector<double>& source,
                        double maxIntercept = 10.0) {
  if(peaks.size() != source.size())
    return false;
  if(peaks.size() <= 2)
    return true;

  for(size_t skippedPoint = 0; skippedPoint < source.size(); ++skippedPoint) {
    std::vector<double> testPeaks;
    std::vector<double> testSource;
    testPeaks.reserve(peaks.size() - 1);
    testSource.reserve(source.size() - 1);
    for(size_t i = 0; i < source.size(); ++i) {
      if(i == skippedPoint)
        continue;
      testPeaks.push_back(peaks[i]);
      testSource.push_back(source[i]);
    }

    double intercept = 0.0;
    double slope = 0.0;
    double chi2 = 0.0;
    if(!LinearFit(testPeaks, testSource, intercept, slope, chi2) ||
       std::abs(intercept) > maxIntercept)
      return false;
  }
  return true;
}

bool IsLocallyConsistentWithSourceLines(const std::vector<double>& peaks,
                                        const std::vector<double>& matchedSource,
                                        const std::vector<double>& allSource,
                                        double intercept,
                                        double slope) {
  if(peaks.size() != matchedSource.size())
    return false;

  for(size_t i = 0; i < peaks.size(); ++i) {
    const double calculatedEnergy = intercept + slope * peaks[i];
    const double assignedDistance = std::abs(calculatedEnergy - matchedSource[i]);
    double nearestDistance = std::numeric_limits<double>::max();
    for(const double sourceEnergy : allSource)
      nearestDistance = std::min(nearestDistance, std::abs(calculatedEnergy - sourceEnergy));

    const double tolerance = std::max(2.0, 0.002 * std::abs(matchedSource[i]));
    if(assignedDistance > nearestDistance + tolerance &&
       assignedDistance > 1.5 * nearestDistance)
      return false;
  }

  return true;
}

double MatchScore(const std::vector<double>& peaks,
                  const std::vector<double>& source,
                  double intercept,
                  double slope,
                  double chi2) {
  double relativeResidual = 0.0;
  for(size_t i = 0; i < peaks.size(); ++i) {
    if(source[i] == 0.0)
      continue;
    const double calculated = intercept + slope * peaks[i];
    const double residual = (calculated - source[i]) / source[i];
    relativeResidual += residual * residual;
  }

  // The raw chi2 keeps the old TCalibrator spirit, while the relative
  // residual term prevents pathological assignments such as matching a
  // candidate near 963 keV to the 867 keV Eu line when the 964 keV line exists.
  return chi2 + 1.0e6 * relativeResidual;
}

std::map<double, double> MatchPeaksToSource(std::vector<double> peaks,
                                            std::vector<double> source) {
  std::map<double, double> match;
  std::sort(peaks.begin(), peaks.end());
  std::sort(source.begin(), source.end());

  const size_t maxPoints = std::min(peaks.size(), source.size());
  for(size_t numPoints = maxPoints; numPoints > 0; --numPoints) {
    double bestScore = std::numeric_limits<double>::max();
    std::map<double, double> bestMatch;

    const auto peakCombinations = Combinations(peaks, numPoints);
    const auto sourceCombinations = Combinations(source, numPoints);
    for(auto peakValues : peakCombinations) {
      for(auto sourceValues : sourceCombinations) {
        if(!HasStableIntercept(peakValues, sourceValues))
          continue;

        if(numPoints >= 2 && peaks.size() > 3) {
          const double peakRatio = peakValues.front() / peakValues.back();
          const double sourceRatio = sourceValues.front() / sourceValues.back();
          if(std::abs(peakRatio - sourceRatio) > 0.02)
            continue;
        }

        const std::vector<double> matchedPeaks = peakValues;
        const std::vector<double> matchedSource = sourceValues;
        peakValues.push_back(0.0);
        sourceValues.push_back(0.0);

        double intercept = 0.0;
        double slope = 0.0;
        double chi2 = 0.0;
        if(!LinearFit(peakValues, sourceValues, intercept, slope, chi2))
          continue;
        if(!IsLocallyConsistentWithSourceLines(matchedPeaks, matchedSource, source, intercept, slope))
          continue;

        const double score = MatchScore(matchedPeaks, matchedSource, intercept, slope, chi2);
        if(score < bestScore) {
          bestMatch.clear();
          for(size_t i = 0; i < numPoints; ++i)
            bestMatch[matchedPeaks[i]] = matchedSource[i];
          bestScore = score;
        }
      }
    }

    if(!bestMatch.empty())
      return bestMatch;
  }

  return match;
}
}

GCalibrator::GCalibrator()
  : TNamed("calibrator", "Energy calibrator") {
}

GCalibrator::GCalibrator(TH1* spectrum,
                         const char* source,
                         bool autoFit,
                         double sigma,
                         double threshold,
                         bool removeBackground)
  : GCalibrator() {
  AddData(spectrum, source, autoFit, sigma, threshold, removeBackground);
}

GCalibrator::GCalibrator(TH1* spectrum,
                         const char* source,
                         double sigma,
                         double threshold,
                         bool removeBackground)
  : GCalibrator(spectrum, source, true, sigma, threshold, removeBackground) {
}

GCalibrator::GCalibrator(TH1* (*histogramGetter)(int),
                         const char* source,
                         bool autoFit,
                         double sigma,
                         double threshold,
                         bool removeBackground)
  : GCalibrator(histogramGetter ? histogramGetter(0) : nullptr,
                source,
                autoFit,
                sigma,
                threshold,
                removeBackground) {
}

GCalibrator::GCalibrator(TH1* (*histogramGetter)(int),
                         const char* source,
                         double sigma,
                         double threshold,
                         bool removeBackground)
  : GCalibrator(histogramGetter, source, true, sigma, threshold, removeBackground) {
}

GCalibrator::~GCalibrator() {
  ClearFit();
  ClearEfficiencyFit();
}

void GCalibrator::ClearFit() {
  delete fFitFunction;
  fFitFunction = nullptr;
  fHasFit = false;
  fFitOrder = 1;
  fParameters[0] = 0.0;
  fParameters[1] = 1.0;
  fParameters[2] = 0.0;
  fParameterErrors[0] = 0.0;
  fParameterErrors[1] = 0.0;
  fParameterErrors[2] = 0.0;
  fChi2 = 0.0;
  fNdf = 0;
}

void GCalibrator::ClearEfficiencyFit() {
  delete fEfficiencyFit;
  fEfficiencyFit = nullptr;
  fEfficiencyGraph.Set(0);
  for(double& parameter : fEfficiencyParameters)
    parameter = 0.0;
}

void GCalibrator::Clear(Option_t* opt) {
  TNamed::Clear(opt);
  ClearFit();
  ClearEfficiencyFit();
  fPeaks.clear();
  fCandidates.clear();
  fCalibrationGraph.Set(0);
  fReferenceSource.clear();
}

int GCalibrator::FindCandidates(TH1* spectrum, double sigma, double threshold) {
  fCandidates.clear();
  if(!spectrum) {
    std::cerr << "GCalibrator::AddData: no histogram was provided.\n";
    return 0;
  }

  TSpectrum search;
  constexpr int maxAutoCandidates = 10;
  int found = search.Search(spectrum, sigma, "nodraw", threshold);
  while(found > maxAutoCandidates) {
    search.Clear();
    threshold += 0.01;
    found = search.Search(spectrum, sigma, "nodraw", threshold);
  }

  if(found <= 0) {
    std::cerr << "GCalibrator::AddData: no candidate peaks found in "
              << spectrum->GetName() << " with sigma = " << sigma
              << " and threshold = " << threshold << ".\n";
    return 0;
  }

  const double* positions = search.GetPositionX();
  for(int i = 0; i < found; ++i) {
    const int bin = spectrum->GetXaxis()->FindBin(positions[i]);
    Candidate candidate;
    candidate.channel = positions[i];
    candidate.channelError = 0.5 * spectrum->GetXaxis()->GetBinWidth(bin);
    candidate.counts = spectrum->GetBinContent(bin);
    fCandidates.push_back(candidate);
  }

  std::sort(fCandidates.begin(), fCandidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              return lhs.channel < rhs.channel;
            });
  return static_cast<int>(fCandidates.size());
}

int GCalibrator::AutoMatchAndFit(TH1* spectrum, const GDecayNucleus& source) {
  if(!spectrum || fCandidates.empty())
    return 0;
  if(source.GetNTransitions() == 0) {
    std::cerr << "GCalibrator::AddData: no decay radiation loaded for "
              << source.GetName() << ".\n";
    return 0;
  }

  std::vector<double> peakPositions;
  peakPositions.reserve(fCandidates.size());
  for(const auto& candidate : fCandidates)
    peakPositions.push_back(candidate.channel);

  std::vector<double> sourceEnergies;
  std::map<double, const GTransition*> transitionForEnergy;
  TIter iter(source.GetTransitionList());
  while(const auto* transition = static_cast<const GTransition*>(iter.Next())) {
    const double energy = transition->GetEnergy();
    if(!std::isfinite(energy))
      continue;
    sourceEnergies.push_back(energy);
    transitionForEnergy[energy] = transition;
  }

  const auto matches = MatchPeaksToSource(peakPositions, sourceEnergies);
  if(matches.empty()) {
    std::cerr << "GCalibrator::AddData: found " << fCandidates.size()
              << " candidate peaks but could not match them to "
              << source.GetName() << " lines.\n";
    return 0;
  }

  const size_t firstNewPeak = fPeaks.size();
  for(const auto& match : matches) {
    const double peakPosition = match.first;
    const double energy = match.second;
    double range = 0.01 * peakPosition;
    if(range < 8.0)
      range = 8.0;

    GGaus* fit = GausFit(spectrum, peakPosition - range, peakPosition + range / 2.0, "no-print");

    double centroid = peakPosition;
    double centroidError = 0.0;
    double area = 0.0;
    double areaError = 0.0;
    double resolution = 0.0;
    if(fit) {
      centroid = fit->GetCentroid();
      centroidError = fit->GetCentroidErr();
      area = fit->GetSum();
      areaError = fit->GetSumErr();
      resolution = centroid != 0.0 ? fit->GetFWHM() / centroid : 0.0;
    }

    double energyError = 0.0;
    double intensity = 0.0;
    double intensityError = 0.0;
    const auto transition = transitionForEnergy.find(energy);
    if(transition != transitionForEnergy.end() && transition->second) {
      energyError = transition->second->GetEnergyUncertainty();
      intensity = transition->second->GetIntensity();
      intensityError = transition->second->GetIntensityUncertainty();
    }

    AddPeak(centroid,
            energy,
            source.GetName(),
            centroidError,
            energyError,
            area,
            areaError,
            intensity,
            intensityError,
            resolution);
  }

  const int added = static_cast<int>(fPeaks.size() - firstNewPeak);
  if(added >= 2 && Fit(fFitOrder))
    Print();
  return added;
}

int GCalibrator::AddData(TH1* spectrum,
                         const GDecayNucleus& source,
                         bool autoFit,
                         double sigma,
                         double threshold,
                         bool removeBackground) {
  fReferenceSource = source.GetName();

  std::unique_ptr<TH1> searchSpectrum;
  TH1* candidateSpectrum = spectrum;
  if(removeBackground && spectrum) {
    TSpectrum backgroundEstimator;
    if(TH1* background = backgroundEstimator.Background(spectrum, 20, "Compton")) {
      searchSpectrum.reset(static_cast<TH1*>(spectrum->Clone(
        Form("%s_gcal_bgsub", spectrum->GetName()))));
      searchSpectrum->SetDirectory(nullptr);
      searchSpectrum->Add(background, -1.0);
      candidateSpectrum = searchSpectrum.get();
      delete background;
    }
  }

  const int candidates = FindCandidates(candidateSpectrum, sigma, threshold);
  if(!autoFit)
    return candidates;
  return AutoMatchAndFit(spectrum, source);
}

int GCalibrator::AddData(TH1* spectrum,
                         const GDecayNucleus& source,
                         double sigma,
                         double threshold,
                         bool removeBackground) {
  return AddData(spectrum, source, true, sigma, threshold, removeBackground);
}

int GCalibrator::AddData(TH1* spectrum,
                         const char* source,
                         bool autoFit,
                         double sigma,
                         double threshold,
                         bool removeBackground) {
  if(!source || !source[0]) {
    std::cerr << "No calibration source provided.\n";
    fCandidates.clear();
    fReferenceSource.clear();
    return 0;
  }

  GDecayNucleus nucleus(source);
  return AddData(spectrum, nucleus, autoFit, sigma, threshold, removeBackground);
}

int GCalibrator::AddData(TH1* spectrum,
                         const char* source,
                         double sigma,
                         double threshold,
                         bool removeBackground) {
  return AddData(spectrum, source, true, sigma, threshold, removeBackground);
}

void GCalibrator::AddPeak(double centroid,
                          double energy,
                          const std::string& source,
                          double centroidError,
                          double energyError,
                          double area,
                          double areaError,
                          double intensity,
                          double intensityError,
                          double resolution) {
  Peak peak;
  peak.centroid = centroid;
  peak.centroidError = centroidError;
  peak.energy = energy;
  peak.energyError = energyError;
  peak.area = area;
  peak.areaError = areaError;
  peak.intensity = intensity;
  peak.intensityError = intensityError;
  peak.resolution = resolution;
  peak.source = source.empty() ? fReferenceSource : source;
  fPeaks.push_back(peak);
  ClearFit();
}

bool GCalibrator::SetPeakEnabled(size_t index, bool enabled) {
  if(index >= fPeaks.size())
    return false;
  fPeaks[index].enabled = enabled;
  ClearFit();
  return true;
}

bool GCalibrator::Fit(int order, bool includeOrigin) {
  order = std::clamp(order, 1, 2);
  const size_t minimumPoints = static_cast<size_t>(order + 1);

  std::vector<const Peak*> points;
  for(const auto& peak : fPeaks) {
    if(peak.enabled)
      points.push_back(&peak);
  }
  if(points.size() + (includeOrigin ? 1 : 0) < minimumPoints) {
    std::cerr << "Need at least " << minimumPoints << " enabled calibration points for order "
              << order << ".\n";
    return false;
  }

  ClearFit();
  fCalibrationGraph.Set(static_cast<int>(points.size() + (includeOrigin ? 1 : 0)));
  int graphIndex = 0;
  if(includeOrigin) {
    fCalibrationGraph.SetPoint(graphIndex, 0.0, 0.0);
    fCalibrationGraph.SetPointError(graphIndex, 0.0, 0.0);
    ++graphIndex;
  }
  for(const auto* peak : points) {
    fCalibrationGraph.SetPoint(graphIndex, peak->centroid, peak->energy);
    // Keep the accepted peak uncertainties on fPeaks for reporting, but fit
    // the calibration graph unweighted by default, matching the old
    // TCalibrator behavior. ROOT's default TGraphErrors fit path does not use
    // x-coordinate errors in this chi2 calculation, so passing centroid errors
    // here only produces warnings without changing the fit in the intended way.
    fCalibrationGraph.SetPointError(graphIndex, 0.0, 0.0);
    ++graphIndex;
  }
  fCalibrationGraph.SetName("calibration_graph");
  fCalibrationGraph.SetTitle("Energy calibration;Channel;Energy");

  const int id = NextFitId();
  fFitFunction = new TF1(Form("calibration_fit_%d", id),
                         order == 1 ? "pol1" : "pol2");
  TFitResultPtr result = fCalibrationGraph.Fit(fFitFunction, "QSN");
  if(!result.Get() || result->Status() != 0) {
    std::cerr << "Calibration fit failed with status "
              << (result.Get() ? result->Status() : -1) << ".\n";
    ClearFit();
    return false;
  }

  fFitOrder = order;
  fHasFit = true;
  for(int i = 0; i <= order; ++i) {
    fParameters[i] = fFitFunction->GetParameter(i);
    fParameterErrors[i] = fFitFunction->GetParError(i);
  }
  fChi2 = result->Chi2();
  fNdf = result->Ndf();
  UpdateResiduals();
  return true;
}

void GCalibrator::UpdateResiduals() {
  for(auto& peak : fPeaks)
    peak.residual = peak.enabled && fHasFit ? Eval(peak.centroid) - peak.energy : 0.0;
}

bool GCalibrator::Check() const {
  if(!fHasFit || fNdf < 0)
    return false;
  return std::all_of(fPeaks.begin(), fPeaks.end(), [](const Peak& peak) {
    return !peak.enabled || std::isfinite(peak.residual);
  });
}

double GCalibrator::Eval(double channel) const {
  return fParameters[0] + fParameters[1] * channel + fParameters[2] * channel * channel;
}

double GCalibrator::GetParameter(int index) const {
  return index >= 0 && index < 3 ? fParameters[index] : std::numeric_limits<double>::quiet_NaN();
}

double GCalibrator::GetParameterError(int index) const {
  return index >= 0 && index < 3 ? fParameterErrors[index] : std::numeric_limits<double>::quiet_NaN();
}

double GCalibrator::GetEffParameter(int index) const {
  return index >= 0 && index < 4 ?
         fEfficiencyParameters[index] :
         std::numeric_limits<double>::quiet_NaN();
}

std::string GCalibrator::PrintEfficency(const char* filename) {
  std::string outputText;
  outputText.append("line\teng\tcounts\tt1/2\tactivity\n");
  outputText.append("--------------------------------------\n");

  int counter = 1;
  for(const auto& peak : fPeaks) {
    outputText.append(Form("%i\t%.02f\t%.02f\t%i\t%.02f\n",
                           counter++,
                           peak.energy,
                           peak.area,
                           100,
                           (peak.intensity / 100.0) * 1e5));
  }
  outputText.append("--------------------------------------\n");

  if(filename && *filename) {
    std::ofstream output(filename);
    output << outputText;
  }

  printf("%s\n", outputText.c_str());
  return outputText;
}

TGraphErrors& GCalibrator::MakeEffGraph(double seconds, double activity, Option_t* opt) {
  TString option(opt ? opt : "");
  TString fitOptions;
  if(option.Contains("Q")) {
    fitOptions.Append("Q");
    option.ReplaceAll("Q", "");
  }

  std::vector<double> energy;
  std::vector<double> energyError;
  std::vector<double> observed;
  std::vector<double> observedError;
  for(const auto& peak : fPeaks) {
    if(!peak.enabled ||
       !std::isfinite(peak.energy) ||
       !std::isfinite(peak.area) ||
       !std::isfinite(peak.intensity) ||
       peak.area <= 0.0 ||
       peak.intensity <= 0.0 ||
       seconds <= 0.0 ||
       activity <= 0.0)
      continue;

    energy.push_back(peak.energy);
    // Match the old TCalibrator behavior: source-energy uncertainties are
    // kept on the peaks, but the efficiency fit does not use x-errors.
    // Passing nonzero coordinate errors to ROOT's default chi2 graph fit only
    // produces "coordinates are not used" warnings.
    energyError.push_back(0.0);
    observed.push_back(peak.area / ((peak.intensity / 100.0) * activity * seconds));
    const double areaError = peak.areaError > 0.0 && std::isfinite(peak.areaError) ?
                             peak.areaError :
                             std::sqrt(peak.area);
    observedError.push_back(observed.back() * (areaError / peak.area));
  }

  ClearEfficiencyFit();
  fEfficiencyGraph = TGraphErrors(static_cast<int>(energy.size()),
                                  energy.data(),
                                  observed.data(),
                                  energyError.data(),
                                  observedError.data());
  fEfficiencyGraph.SetName("efficiency_graph");
  fEfficiencyGraph.SetTitle("Gamma efficiency;Energy;Efficiency");

  if(fEfficiencyGraph.GetN() > 0) {
    static int counter = 0;
    fEfficiencyFit = new TF1(Form("eff_fit_%i", counter++), GFunctions::GammaEff, 0.1, 4000.0, 4);
    fEfficiencyFit->SetParName(0, "zeroth");
    fEfficiencyFit->SetParName(1, "first");
    fEfficiencyFit->SetParName(2, "second");
    fEfficiencyFit->SetParName(3, "inverse");
    fEfficiencyFit->SetParameter(0, -1.0);
    fEfficiencyFit->SetParameter(1, 0.5);
    fEfficiencyFit->SetParameter(2, -0.2);
    fEfficiencyFit->SetParameter(3, -50.0);
    fEfficiencyGraph.Fit(fEfficiencyFit, fitOptions.Data());
    for(int i = 0; i < 4; ++i)
      fEfficiencyParameters[i] = fEfficiencyFit->GetParameter(i);
  }

  if(option.Contains("draw", TString::kIgnoreCase)) {
    TVirtualPad* current = gPad;
    new GCanvas;
    fEfficiencyGraph.Draw("AP");
    if(fEfficiencyFit)
      fEfficiencyFit->Draw("same");
    if(current)
      current->cd();
  }

  if(fEfficiencyFit) {
    for(int i = 0; i < fEfficiencyGraph.GetN(); ++i) {
      double x = 0.0;
      double y = 0.0;
      fEfficiencyGraph.GetPoint(i, x, y);
      printf("[%.1f] Observed  = %.04f  | Calculated = %.04f  |  per diff = %.2f\n",
             x,
             y,
             fEfficiencyFit->Eval(x),
             y != 0.0 ? (std::abs(y - fEfficiencyFit->Eval(x)) / y) * 100.0 : 0.0);
    }
  }

  return fEfficiencyGraph;
}

bool GCalibrator::SaveEffGraph(std::string datafile, std::string fitfile) {
  if(!fEfficiencyFit)
    return false;

  if(datafile.empty())
    datafile = Form("%s_data.dat", GetName());
  if(fitfile.empty())
    fitfile = Form("%s_fit.dat", GetName());

  std::ofstream dataOutput(datafile);
  if(!dataOutput)
    return false;
  for(int i = 0; i < fEfficiencyGraph.GetN(); ++i) {
    double x = 0.0;
    double y = 0.0;
    fEfficiencyGraph.GetPoint(i, x, y);
    dataOutput << x << "\t" << y << "\t" << fEfficiencyGraph.GetErrorY(i) << '\n';
  }
  dataOutput << '\n';

  std::ofstream fitOutput(fitfile);
  if(!fitOutput)
    return false;
  double xmin = 0.0;
  double xmax = 0.0;
  fEfficiencyFit->GetRange(xmin, xmax);
  const double range = xmax - xmin;
  while(xmin < xmax) {
    fitOutput << xmin << "\t" << fEfficiencyFit->Eval(xmin) << '\n';
    xmin += range / 10000.0;
  }
  fitOutput << '\n';
  return true;
}

TH1D* GCalibrator::ApplyCalibration(const TH1* source, const char* name, const char* title) const {
  if(!source || !fHasFit)
    return nullptr;

  std::vector<double> edges;
  edges.reserve(source->GetNbinsX() + 1);
  for(int bin = 1; bin <= source->GetNbinsX() + 1; ++bin)
    edges.push_back(Eval(source->GetXaxis()->GetBinLowEdge(bin)));
  if(!std::is_sorted(edges.begin(), edges.end())) {
    std::cerr << "Calibration is not monotonic across this spectrum.\n";
    return nullptr;
  }

  const std::string outputName = name && *name ? name : std::string(source->GetName()) + "_cal";
  const std::string outputTitle = title && *title ? title : source->GetTitle();
  auto* calibrated = new TH1D(outputName.c_str(), outputTitle.c_str(),
                              source->GetNbinsX(), edges.data());
  calibrated->SetDirectory(nullptr);
  calibrated->Sumw2();
  for(int bin = 0; bin <= source->GetNbinsX() + 1; ++bin) {
    calibrated->SetBinContent(bin, source->GetBinContent(bin));
    calibrated->SetBinError(bin, source->GetBinError(bin));
  }
  calibrated->SetEntries(source->GetEntries());
  calibrated->GetXaxis()->SetTitle("Energy");
  calibrated->GetYaxis()->SetTitle(source->GetYaxis()->GetTitle());
  calibrated->GetListOfFunctions()->Add(new TParameter<double>("C0", fParameters[0]));
  calibrated->GetListOfFunctions()->Add(new TParameter<double>("C1", fParameters[1]));
  calibrated->GetListOfFunctions()->Add(new TParameter<double>("C2", fParameters[2]));
  calibrated->GetListOfFunctions()->Add(new TObjString(("calibrator=" + std::string(GetName())).c_str()));
  return calibrated;
}

bool GCalibrator::Save(const char* filename) const {
  if(!filename || !*filename || !fHasFit)
    return false;
  std::ofstream output(filename);
  if(!output)
    return false;
  output << std::setprecision(16)
         << "C0 = " << fParameters[0] << '\n'
         << "C1 = " << fParameters[1] << '\n'
         << "C2 = " << fParameters[2] << '\n'
         << "order = " << fFitOrder << '\n'
         << "chi2 = " << fChi2 << '\n'
         << "ndf = " << fNdf << '\n'
         << "source = " << fReferenceSource << '\n';
  return true;
}

bool GCalibrator::Load(const char* filename) {
  if(!filename || !*filename)
    return false;
  std::ifstream input(filename);
  if(!input)
    return false;

  ClearFit();
  std::string line;
  while(std::getline(input, line)) {
    std::stringstream stream(line);
    std::string key;
    char equal = 0;
    if(!(stream >> key >> equal) || equal != '=')
      continue;
    if(key == "C0") stream >> fParameters[0];
    else if(key == "C1") stream >> fParameters[1];
    else if(key == "C2") stream >> fParameters[2];
    else if(key == "order") stream >> fFitOrder;
    else if(key == "chi2") stream >> fChi2;
    else if(key == "ndf") stream >> fNdf;
    else if(key == "source") {
      std::getline(stream, fReferenceSource);
      if(!fReferenceSource.empty() && fReferenceSource.front() == ' ')
        fReferenceSource.erase(0, 1);
    }
  }
  fFitOrder = std::clamp(fFitOrder, 1, 2);
  fHasFit = true;
  return true;
}

void GCalibrator::Draw(Option_t* option) {
  if(fCalibrationGraph.GetN() == 0)
    return;
  fCalibrationGraph.Draw(option && *option ? option : "AP");
  if(fFitFunction)
    fFitFunction->Draw("same");
}

void GCalibrator::Print(Option_t*) const {
  printf("\t%2senergy%10scent%10scalc%10sarea%7snuc%8sintensity\n",
         "", "", "", "", "", "");

  int counter = 0;
  for(const auto& peak : fPeaks) {
    const double calculatedEnergy = fHasFit ? Eval(peak.centroid)
                                            : std::numeric_limits<double>::quiet_NaN();
    const double percentDifference = peak.energy != 0.0 ?
                                     std::abs(calculatedEnergy - peak.energy) / peak.energy * 100.0 :
                                     std::numeric_limits<double>::quiet_NaN();
    printf("%i:\t%7.02f%16.02f%8.2f%3s[%%%3.2f]%16.02f%8s%16.04f%s\n",
           counter++,
           peak.energy,
           peak.centroid,
           calculatedEnergy,
           "",
           percentDifference,
           peak.area,
           peak.source.c_str(),
           peak.intensity,
           peak.enabled ? "" : " [disabled]");
  }
  printf("-------------------------------\n");
}
