#ifndef __GLEARNING_H__
#define __GLEARNING_H__

#include <TObject.h>

#include <GCalibration.h>

#include <string>
#include <vector>

class GF1;
class TH1;
class TH1D;
class TF1;
class TGraphErrors;

struct GLearningLine {
  double energy = 0;
  double energyError = 0;

  // Stored internally as probability from 0 to 1.
  double intensity = 0;
  double intensityError = 0;

  bool recommended = true;
};

struct GLearningPeak {
  double channel = 0;
  double channelError = 0;

  double energy = 0;
  double energyError = 0;

  double intensity = 0;
  double intensityError = 0;

  double area = 0;
  double areaError = 0;

  double efficiency = 0;
  double efficiencyError = 0;

  double chi2Ndf = 0;

  bool enabled = true;
};

class GLearning : public TObject {
  public:
    GLearning(TH1* spectrum=nullptr);
    virtual ~GLearning();

    void Clear(Option_t* opt="") override;

    void SetSpectrum(TH1* spectrum);
    TH1* GetSpectrum() const { return fSpectrum; }

    bool SetSource(const char* source);
    const char* GetSource() const { return fSource.c_str(); }

    void SetActivity(double activityBq) {
      fReferenceActivity = activityBq;
    }

    void SetHalfLifeYears(double years);
    void SetElapsedDays(double days);

    void SetLiveTime(double seconds) {
      fLiveTime = seconds;
    }

    double GetLiveTime() const {
      return fLiveTime;
    }

    double GetCurrentActivity() const;

    void Start() const;
    void PrintSourceLines() const;
    void PrintPeaks() const;

    bool AcceptFit(GF1* fit, double knownEnergy);
    bool AcceptLastFit(double knownEnergy);

    bool FitEnergyCalibration(int order=1);
    TGraphErrors* DrawEnergyCalibration();
    TH1D* MakeCalibratedSpectrum(const char* name="calibrated_spectrum");

    bool BuildEfficiencyPoints();
    TGraphErrors* DrawEfficiency();

    TF1* FitEfficiency();
    TF1* FitEfficiencyPolynomial(int order=3);

    bool EnablePoint(size_t index, bool enabled=true);
    bool Save(const char* filename) const;

    const std::vector<GLearningPeak>& GetPeaks() const {
      return fPeaks;
    }

    const GCalibration& GetCalibration() const {
      return fCalibration;
    }

  private:
    const GLearningLine* FindLine(double energy) const;
    void ClearGraphs();

    TH1* fSpectrum;

    std::string fSource;
    std::vector<GLearningLine> fLines;
    std::vector<GLearningPeak> fPeaks;

    GCalibration fCalibration;

    double fReferenceActivity;
    double fHalfLifeSeconds;
    double fElapsedSeconds;
    double fLiveTime;

    TGraphErrors* fEfficiencyGraph; //!
    TF1* fEfficiencyFit;            //!

  ClassDefOverride(GLearning, 1)
};

#endif
