#include <GNucleus.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

ClassImp(GNucleus)

namespace {
constexpr double kAtomicMassUnitMeV = 931.494043;

std::string Trim(std::string value) {
  const auto first = std::find_if_not(value.begin(), value.end(),
                                      [](unsigned char c) { return std::isspace(c); });
  const auto last = std::find_if_not(value.rbegin(), value.rend(),
                                     [](unsigned char c) { return std::isspace(c); }).base();
  return first < last ? std::string(first, last) : std::string();
}

std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return value;
}

std::string ElementSymbol(std::string value) {
  value = Trim(value);
  value.erase(std::remove_if(value.begin(), value.end(),
                             [](unsigned char c) { return std::isdigit(c); }), value.end());
  if(value.empty())
    return value;
  value = Lower(value);
  value.front() = std::toupper(static_cast<unsigned char>(value.front()));
  return value;
}

bool IsReadable(const std::string& filename) {
  std::ifstream input(filename);
  return input.good();
}
}

std::string& GNucleus::MassFilePath() {
  static std::string path;
  return path;
}

std::string GNucleus::DefaultMassFile() {
  if(const char* configured = std::getenv("GROOT_MASS_DATA"))
    return configured;

  std::vector<std::string> candidates;
  if(const char* gsys = std::getenv("GSYS")) {
    const std::string root = gsys;
    candidates.push_back(root + "/data/mass.dat");
    candidates.push_back(root + "/build/data/mass.dat");
    candidates.push_back(root + "/libraries/GTools/mass.dat");
  }
  candidates.emplace_back("libraries/GTools/mass.dat");

  for(const auto& candidate : candidates) {
    if(IsReadable(candidate))
      return candidate;
  }
  return candidates.back();
}

void GNucleus::SetMassFile(const std::string& filename) {
  MassFilePath() = filename;
}

std::string GNucleus::GetMassFile() {
  if(MassFilePath().empty())
    MassFilePath() = DefaultMassFile();
  return MassFilePath();
}

std::string GNucleus::SortName(const char* name) {
  if(!name)
    return {};

  std::string value = Trim(name);
  value.erase(std::remove_if(value.begin(), value.end(),
                             [](unsigned char c) { return std::isspace(c); }), value.end());
  if(value == "p") value = "1H";
  if(value == "d") value = "2H";
  if(value == "t") value = "3H";
  if(value == "a") value = "4He";

  const auto firstDigit = std::find_if(value.begin(), value.end(),
                                       [](unsigned char c) { return std::isdigit(c); });
  if(firstDigit == value.end())
    return {};

  std::string mass;
  std::string symbol;
  if(firstDigit == value.begin()) {
    const auto firstLetter = std::find_if_not(value.begin(), value.end(),
                                              [](unsigned char c) { return std::isdigit(c); });
    mass.assign(value.begin(), firstLetter);
    symbol.assign(firstLetter, value.end());
  } else {
    symbol.assign(value.begin(), firstDigit);
    mass.assign(firstDigit, value.end());
  }
  if(mass.empty() || symbol.empty() ||
     !std::all_of(mass.begin(), mass.end(), [](unsigned char c) { return std::isdigit(c); }))
    return {};

  symbol = Lower(symbol);
  symbol.front() = std::toupper(static_cast<unsigned char>(symbol.front()));
  return mass + symbol;
}

GNucleus::GNucleus(const char* name) {
  const std::string canonical = SortName(name);
  if(canonical.empty())
    return;
  TNamed::SetName(canonical.c_str());
  LoadMassData(nullptr);
}

GNucleus::GNucleus(int Z, int N, double mass, const char* symbol)
  : fN(N), fZ(Z), fMass(mass), fSymbol(ElementSymbol(symbol ? symbol : "")), fValid(Z >= 0 && N >= 0) {
  UpdateName();
}

GNucleus::GNucleus(int Z, int N, const char* massFile)
  : fN(N), fZ(Z) {
  LoadMassData(massFile);
}

bool GNucleus::LoadMassData(const char* massFile) {
  const std::string filename = massFile ? massFile : GetMassFile();
  std::ifstream input(filename);
  if(!input) {
    std::cerr << "Could not open mass data file: " << filename << '\n';
    return false;
  }

  const std::string requested = GetName()[0] ? Lower(SortName(GetName())) : std::string();
  std::string line;
  while(std::getline(input, line)) {
    std::stringstream stream(line);
    int N = 0;
    int Z = 0;
    std::string symbol;
    double massExcessKeV = 0.0;
    if(!(stream >> N >> Z >> symbol >> massExcessKeV))
      continue;

    const std::string candidate = Lower(symbol);
    const bool matchesName = !requested.empty() && candidate == requested;
    const bool matchesNumbers = requested.empty() && N == fN && Z == fZ;
    if(!matchesName && !matchesNumbers)
      continue;

    fN = N;
    fZ = Z;
    fSymbol = ElementSymbol(symbol);
    fMassExcess = massExcessKeV / 1000.0;
    SetMass();
    fValid = true;
    UpdateName();
    return true;
  }

  std::cerr << "Nucleus was not found in mass data: "
            << (requested.empty() ? std::to_string(fN + fZ) : requested) << '\n';
  return false;
}

void GNucleus::SetMass() {
  fMass = kAtomicMassUnitMeV * GetA() + fMassExcess;
}

void GNucleus::SetSymbol(const char* symbol) {
  fSymbol = ElementSymbol(symbol ? symbol : "");
  UpdateName();
}

void GNucleus::UpdateName() {
  if(fSymbol.empty() || GetA() <= 0)
    return;
  SetName((fSymbol + std::to_string(GetA())).c_str());
}

double GNucleus::GetEnergyFromBeta(double beta) const {
  if(fMass <= 0.0 || std::abs(beta) >= 1.0)
    return 0.0;
  return fMass * (1.0 / std::sqrt(1.0 - beta * beta) - 1.0);
}

double GNucleus::GetBetaFromEnergy(double energyMeV) const {
  if(fMass <= 0.0 || energyMeV < 0.0)
    return 0.0;
  const double gamma = energyMeV / fMass + 1.0;
  return std::sqrt(1.0 - 1.0 / (gamma * gamma));
}

double GNucleus::GetRadius() const {
  const int A = GetA();
  if(A <= 0)
    return 0.0;
  return 1.12 * std::pow(A, 1.0 / 3.0) - 0.94 * std::pow(A, -1.0 / 3.0);
}

void GNucleus::Print(Option_t*) const {
  if(!fValid) {
    printf("Invalid nucleus\n");
    return;
  }
  printf("Nucleus: %s (Z=%d, N=%d, mass=%.6f MeV, mass excess=%.6f MeV)\n",
         GetName(), fZ, fN, fMass, fMassExcess);
}
