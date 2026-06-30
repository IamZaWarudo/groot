#ifndef G_DECAY_NUCLEUS_H
#define G_DECAY_NUCLEUS_H

#include <string>

#include <GNucleus.h>
#include <TList.h>

class GTransition;

// A parent nucleus together with the gamma radiation observed following its
// decay. Daughter branches remain an implementation detail of the dataset.
class GDecayNucleus : public GNucleus {
  public:
    GDecayNucleus();
    explicit GDecayNucleus(const char* parent);
    ~GDecayNucleus() override;

    static void SetDecayDataDirectory(const std::string& directory);
    static std::string GetDecayDataDirectory();

    bool LoadDecayRadiation(const char* filename = nullptr);
    void ClearDecayRadiation();
    bool HasDecayRadiation() const { return fDecayRadiationLoaded; }
    const char* GetDecayDataFile() const { return fDecayDataFile.c_str(); }

    void AddTransition(double energy,
                       double intensity = 0.0,
                       double energyUncertainty = 0.0,
                       double intensityUncertainty = 0.0);
    void AddTransition(GTransition* transition);
    GTransition* GetTransition(int index);
    const GTransition* GetTransition(int index) const;
    int GetNTransitions() const { return fDecayRadiation.GetSize(); }
    TList* GetTransitionList() { return &fDecayRadiation; }
    const TList* GetTransitionList() const { return &fDecayRadiation; }

    void Print(Option_t* opt = "") const override;

  private:
    static std::string& DecayDataDirectory();
    static std::string DefaultDecayDataDirectory();

    TList fDecayRadiation;
    std::string fDecayDataFile;
    bool fDecayRadiationLoaded = false;

  ClassDefOverride(GDecayNucleus, 1)
};

#endif
