#ifndef G_NUCLEUS_H
#define G_NUCLEUS_H

#include <string>

#include <TNamed.h>

// Physical isotope identity and mass data. Decay radiation belongs in
// GDecayNucleus, which is intentionally a separate derived type.
class GNucleus : public TNamed {
  public:
    GNucleus() = default;
    explicit GNucleus(const char* name);
    GNucleus(int Z, int N, double mass, const char* symbol);
    GNucleus(int Z, int N, const char* massFile = nullptr);
    ~GNucleus() override = default;

    static std::string SortName(const char* name);
    static void SetMassFile(const std::string& filename);
    static std::string GetMassFile();

    void SetZ(int Z) { fZ = Z; }
    void SetN(int N) { fN = N; }
    void SetMassExcess(double massExcess) { fMassExcess = massExcess; }
    void SetMass(double mass) { fMass = mass; }
    void SetMass();
    void SetSymbol(const char* symbol);

    int GetZ() const { return fZ; }
    int GetN() const { return fN; }
    int GetA() const { return fN + fZ; }
    double GetMassExcess() const { return fMassExcess; }
    double GetMass() const { return fMass; }
    const char* GetSymbol() const { return fSymbol.c_str(); }
    bool IsValid() const { return fValid; }

    double GetEnergyFromBeta(double beta) const;
    double GetBetaFromEnergy(double energyMeV) const;
    double GetRadius() const;

    void Print(Option_t* opt = "") const override;

  protected:
    void UpdateName();
    bool LoadMassData(const char* massFile);

  private:
    static std::string& MassFilePath();
    static std::string DefaultMassFile();

    int fN = 0;
    int fZ = 0;
    double fMass = 0.0;
    double fMassExcess = 0.0;
    std::string fSymbol;
    bool fValid = false;

  ClassDefOverride(GNucleus, 2)
};

#endif
