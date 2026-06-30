#include <GDecayNucleus.h>

#include <GTransition.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

ClassImp(GDecayNucleus)

namespace {
std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return value;
}

bool IsReadable(const std::string& filename) {
  std::ifstream input(filename);
  return input.good();
}
}

std::string& GDecayNucleus::DecayDataDirectory() {
  static std::string directory;
  return directory;
}

std::string GDecayNucleus::DefaultDecayDataDirectory() {
  if(const char* configured = std::getenv("GROOT_DECAY_DATA_DIR"))
    return configured;

  std::vector<std::string> candidates;
  if(const char* gsys = std::getenv("GSYS")) {
    const std::string root = gsys;
    candidates.push_back(root + "/data/SourceData");
    candidates.push_back(root + "/build/data/SourceData");
    candidates.push_back(root + "/libraries/SourceData");
  }
  candidates.emplace_back("libraries/SourceData");

  for(const auto& candidate : candidates) {
    if(IsReadable(candidate + "/co60.sou"))
      return candidate;
  }
  return candidates.back();
}

void GDecayNucleus::SetDecayDataDirectory(const std::string& directory) {
  DecayDataDirectory() = directory;
}

std::string GDecayNucleus::GetDecayDataDirectory() {
  if(DecayDataDirectory().empty())
    DecayDataDirectory() = DefaultDecayDataDirectory();
  return DecayDataDirectory();
}

GDecayNucleus::GDecayNucleus() {
  fDecayRadiation.SetOwner(kTRUE);
}

GDecayNucleus::GDecayNucleus(const char* parent)
  : GNucleus(parent) {
  fDecayRadiation.SetOwner(kTRUE);
  if(IsValid())
    LoadDecayRadiation();
}

GDecayNucleus::~GDecayNucleus() {
  ClearDecayRadiation();
}

void GDecayNucleus::ClearDecayRadiation() {
  fDecayRadiation.Delete();
  fDecayDataFile.clear();
  fDecayRadiationLoaded = false;
}

bool GDecayNucleus::LoadDecayRadiation(const char* filename) {
  ClearDecayRadiation();

  if(!IsValid())
    return false;

  std::string sourceFile;
  if(filename && *filename) {
    sourceFile = filename;
  } else {
    sourceFile = GetDecayDataDirectory() + "/" +
                 Lower(std::string(GetSymbol())) + std::to_string(GetA()) + ".sou";
  }

  std::ifstream input(sourceFile);
  if(!input) {
    std::cerr << "Could not open decay radiation data: " << sourceFile << '\n';
    return false;
  }

  std::string line;
  while(std::getline(input, line)) {
    const auto comment = line.find("//");
    if(comment != std::string::npos)
      line.erase(comment);

    std::stringstream stream(line);
    double energy = 0.0;
    if(!(stream >> energy))
      continue;

    const double unknown = std::numeric_limits<double>::quiet_NaN();
    double energyUncertainty = unknown;
    double intensity = unknown;
    double intensityUncertainty = unknown;
    stream >> energyUncertainty >> intensity >> intensityUncertainty;
    AddTransition(energy, intensity, energyUncertainty, intensityUncertainty);
  }

  fDecayDataFile = sourceFile;
  fDecayRadiationLoaded = GetNTransitions() > 0;
  if(!fDecayRadiationLoaded)
    std::cerr << "No decay gamma lines were found in: " << sourceFile << '\n';
  return fDecayRadiationLoaded;
}

void GDecayNucleus::AddTransition(double energy,
                                  double intensity,
                                  double energyUncertainty,
                                  double intensityUncertainty) {
  auto* transition = new GTransition;
  transition->SetEnergy(energy);
  transition->SetEnergyUncertainty(energyUncertainty);
  transition->SetIntensity(intensity);
  transition->SetIntensityUncertainty(intensityUncertainty);
  AddTransition(transition);
}

void GDecayNucleus::AddTransition(GTransition* transition) {
  if(transition)
    fDecayRadiation.Add(transition);
}

GTransition* GDecayNucleus::GetTransition(int index) {
  return static_cast<GTransition*>(fDecayRadiation.At(index));
}

const GTransition* GDecayNucleus::GetTransition(int index) const {
  return static_cast<const GTransition*>(fDecayRadiation.At(index));
}

void GDecayNucleus::Print(Option_t* opt) const {
  GNucleus::Print(opt);
  printf("Decay gamma lines: %d",GetNTransitions());
  if(!fDecayDataFile.empty())
    printf(" (%s)",fDecayDataFile.c_str());
  printf("\n");

  TIter next(&fDecayRadiation);
  while(auto* transition = static_cast<GTransition*>(next()))
    transition->Print(opt);
}
