#ifndef G_CALIBRATOR_H
#define G_CALIBRATOR_H

#include <string>
#include <vector>

#include <TNamed.h>
#include <TGraphErrors.h>

class GDecayNucleus;
class TH1;
class TH1D;
class TF1;

class GCalibrator : public TNamed {
  public:
    // A peak found in a spectrum. In automatic mode these are matched to
    // source lines and promoted to accepted calibration peaks.
    struct Candidate {
      double channel = 0.0;
      double channelError = 0.0;
      double counts = 0.0;
    };

    // An analyst-accepted channel-to-energy correspondence.
    struct Peak {
      double centroid = 0.0;
      double centroidError = 0.0;
      double energy = 0.0;
      double energyError = 0.0;
      double area = 0.0;
      double areaError = 0.0;
      double intensity = 0.0;
      double intensityError = 0.0;
      double resolution = 0.0;
      bool enabled = true;
      double residual = 0.0;
      std::string source;
    };

    GCalibrator();
    GCalibrator(TH1* spectrum,
                const char* source,
                bool autoFit = true,
                double sigma = 2.0,
                double threshold = 0.05,
                bool removeBackground = false);
    GCalibrator(TH1* spectrum,
                const char* source,
                double sigma,
                double threshold = 0.05,
                bool removeBackground = false);
    GCalibrator(TH1* (*histogramGetter)(int),
                const char* source,
                bool autoFit = true,
                double sigma = 2.0,
                double threshold = 0.05,
                bool removeBackground = false);
    GCalibrator(TH1* (*histogramGetter)(int),
                const char* source,
                double sigma,
                double threshold = 0.05,
                bool removeBackground = false);
    ~GCalibrator() override;

    void Clear(Option_t* opt = "") override;
    void Print(Option_t* opt = "") const override;
    void Draw(Option_t* opt = "") override;

    int AddData(TH1* spectrum,
                const GDecayNucleus& source,
                bool autoFit = true,
                double sigma = 2.0,
                double threshold = 0.05,
                bool removeBackground = false);
    int AddData(TH1* spectrum,
                const GDecayNucleus& source,
                double sigma,
                double threshold = 0.05,
                bool removeBackground = false);
    int AddData(TH1* spectrum,
                const char* source,
                bool autoFit = true,
                double sigma = 2.0,
                double threshold = 0.05,
                bool removeBackground = false);
    int AddData(TH1* spectrum,
                const char* source,
                double sigma,
                double threshold = 0.05,
                bool removeBackground = false);

    void AddPeak(double centroid,
                 double energy,
                 const std::string& source = "",
                 double centroidError = 0.0,
                 double energyError = 0.0,
                 double area = 0.0,
                 double areaError = 0.0,
                 double intensity = 0.0,
                 double intensityError = 0.0,
                 double resolution = 0.0);
    bool SetPeakEnabled(size_t index, bool enabled);

    bool Fit(int order = 1, bool includeOrigin = false);
    bool Check() const;
    bool HasFit() const { return fHasFit; }
    double Eval(double channel) const;
    double GetParameter(int index) const;
    double GetParameterError(int index) const;
    double GetChi2() const { return fChi2; }
    int GetNdf() const { return fNdf; }
    int GetFitOrder() const { return fFitOrder; }
    double GetEffParameter(int index) const;
    double GetEfficiencyParameter(int index) const { return GetEffParameter(index); }

    TGraphErrors& MakeEffGraph(double seconds = 3600.0,
                               double activity = 100000.0,
                               Option_t* opt = "draw");
    TGraphErrors* MakeEfficiencyGraph(double seconds = 3600.0,
                                      double activity = 100000.0,
                                      Option_t* opt = "draw") {
      return &MakeEffGraph(seconds, activity, opt);
    }
    TGraphErrors* MakeEffGraphPtr(double seconds = 3600.0,
                                  double activity = 100000.0,
                                  Option_t* opt = "draw") {
      return &MakeEffGraph(seconds, activity, opt);
    }
    TGraphErrors* MakeEfficiencyGraphPtr(double seconds = 3600.0,
                                         double activity = 100000.0,
                                         Option_t* opt = "draw") {
      return MakeEfficiencyGraph(seconds, activity, opt);
    }
    bool SaveEffGraph(std::string datafile = "", std::string fitfile = "");
    bool SaveEfficiencyGraph(std::string datafile = "", std::string fitfile = "") {
      return SaveEffGraph(datafile, fitfile);
    }
    std::string PrintEfficency(const char* filename = "");
    std::string PrintEfficiency(const char* filename = "") { return PrintEfficency(filename); }

    TH1D* ApplyCalibration(const TH1* source,
                           const char* name = nullptr,
                           const char* title = nullptr) const;

    bool Save(const char* filename) const;
    bool Load(const char* filename);

    size_t GetNPeaks() const { return fPeaks.size(); }
    const Peak& GetPeak(size_t index) const { return fPeaks.at(index); }
    const std::vector<Peak>& GetPeaks() const { return fPeaks; }

    size_t GetNCandidates() const { return fCandidates.size(); }
    const Candidate& GetCandidate(size_t index) const { return fCandidates.at(index); }
    const std::vector<Candidate>& GetCandidates() const { return fCandidates; }
    const char* GetReferenceSource() const { return fReferenceSource.c_str(); }

    const TGraphErrors* GetCalibrationGraph() const { return &fCalibrationGraph; }
    const TF1* GetFitFunction() const { return fFitFunction; }
    TGraphErrors* EffGraph() { return &fEfficiencyGraph; }
    const TGraphErrors* EffGraph() const { return &fEfficiencyGraph; }
    TGraphErrors* GetEfficiencyGraph() { return &fEfficiencyGraph; }
    const TGraphErrors* GetEfficiencyGraph() const { return &fEfficiencyGraph; }
    TF1* EffFit() { return fEfficiencyFit; }
    const TF1* EffFit() const { return fEfficiencyFit; }
    TF1* GetEfficiencyFit() { return fEfficiencyFit; }
    const TF1* GetEfficiencyFit() const { return fEfficiencyFit; }

  private:
    void ClearFit();
    void ClearEfficiencyFit();
    void UpdateResiduals();
    int FindCandidates(TH1* spectrum, double sigma, double threshold);
    int AutoMatchAndFit(TH1* spectrum, const GDecayNucleus& source);

    std::vector<Peak> fPeaks;
    std::vector<Candidate> fCandidates;
    TGraphErrors fCalibrationGraph;
    TGraphErrors fEfficiencyGraph;
    TF1* fFitFunction = nullptr; //! transient display/fit helper
    TF1* fEfficiencyFit = nullptr; //! transient efficiency fit helper
    std::string fReferenceSource;
    int fFitOrder = 1;
    bool fHasFit = false;
    double fParameters[3] = {0.0, 1.0, 0.0};
    double fParameterErrors[3] = {0.0, 0.0, 0.0};
    double fChi2 = 0.0;
    int fNdf = 0;
    double fEfficiencyParameters[4] = {0.0, 0.0, 0.0, 0.0};

  ClassDefOverride(GCalibrator, 0)
};

#endif
