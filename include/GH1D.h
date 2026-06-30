#ifndef __GH1D_H__
#define __GH1D_H__

#include<string>
#include <cmath>

#include<TH1.h>
#include<TH2.h>
//#include<TQObject.h>

class TF1;
class TVirtualPad;

class GH1D : public TH1D { //, public TQObject {
  public:
    GH1D();
    GH1D(std::string name,int nbinsx,const double *xbins);
    GH1D(std::string name,int nbinsx,const float  *xbins);
    GH1D(std::string name,int nbinsx,double xlow, double xhigh);

    GH1D(const char *name,const char *title,int nbinsx,const double *xbins);
    GH1D(const char *name,const char *title,int nbinsx,const float  *xbins);
    GH1D(const char *name,const char *title,int nbinsx,double xlow, double xhigh);

    //GH1D(const GH1D &h1d);
    GH1D(const TH1D &h1d);
    GH1D(const TH1F &h1f);
    //GH1D(const TH1 *h1);

    GH1D(const TVectorD &v);
    
    virtual ~GH1D(); 

    void Draw(Option_t *opt="") override;
    TH1* DrawCopy(Option_t *opt="", const char *name_postfix="_copy") const override; 
    TH1* DrawNormalized(Option_t *opt="", double norm=1) const override; 

    void Paint(Option_t *opt="") override;
    int DistancetoPrimitive(int px, int py) override;

    void SetShowResiduals(bool flag=true);
    void ToggleResiduals();
    bool ShowResiduals() const { return fShowResiduals; }
    bool IsShowingResiduals() const { return ShowResiduals(); }
    void SetResidualFit(TF1 *fit,double xlow,double xhigh);
    void ClearResidual();
    void UpdateResidualDisplay(TVirtualPad *pad=nullptr);

    TH1* Rebin(int ngroup=2,const char *newname="",const double *xbins=nullptr) override;
    void Unbin(int ngroup=-1);

    

    void SetBackground(int niter=12,Option_t* opt="compton");
    void ShowBackground();
    void ToggleBackground();
    TH1* GetBackground() const { return fBg; }

    //TH1  *ShowBackground(int niter=12,Option_t* opt="same") override;
    //void Background();
    

    void PeakSearch(double threshold=0.05,double sigma=1,Option_t *opt="") { }  

    int GetNbinsOriginal() const { return fOriginalBins; }

    bool  Add(const TH1 *h1, Double_t c1=1) override;
    //GH1D* AddNormalized(const TH1 *h1, Double_t c1=1);
    enum EstatusBits {
      kBackgroundRemoved = BIT(22),
      kProjectionX       = BIT(23)
    };

    void Print(Option_t *opt="") const override;

  public:
    void SetParent(TH2 *h)    { fParent = h;    }
		TH2* GetParent() const    { return fParent;  }

		void SetScale(double scale) { fScale = scale; } 

    bool IsNormalized() const { return fIsNormalized; }
    void Normalize();

    void UpdateFunctions();
    int  ShowPeaks(double thresh=.10,double sigma=1.0);
    bool RemovePeaks();

  private:
    void Init();
		void SetOriginal();
		void ResetToOriginal();
    TH1D *MakeResidualHist() const;

  private:
    //owned pointers
    TH1D *fOriginal;    
    TH1D *fBg;


		//TH1D *fSubtract;
		double fScale;

    int  fOriginalBins;
    bool fIsNormalized;
    bool fShowResiduals;
    TF1  *fResidualFit;  //! non-owning; fit is owned by the caller/histogram
    TH1D *fResidualHist; //!
    double fResidualXLow;
    double fResidualXHigh;
    
    //borrowed pointer
    TH2  *fParent;


  ClassDefOverride(GH1D,102)
};

#endif
