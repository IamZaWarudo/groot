#include<GH1D.h>

#include<GCommands.h>

#include<algorithm>
#include<cstdint>
#include<cstdio>
#include<cstring>
#include<limits>

#include<TSpectrum.h>
#include<TPolyMarker.h>
#include<TText.h>
#include<TString.h>
#include<TVirtualPad.h>
#include<TCanvas.h>
#include<TPad.h>
#include<TROOT.h>
#include<TBox.h>
#include<TF1.h>
#include<TGraph.h>
#include<TLine.h>
#include<TList.h>
#include<TArrayD.h>
#include<vector>

#include "GROI.h"
#include "TMath.h"

namespace {

TString ResidualPadName(const GH1D *hist,const char *suffix) {
  return Form("__gh1d_residual_%zx_%s",
              reinterpret_cast<std::uintptr_t>(hist),suffix);
}

TVirtualPad *FindPadContaining(TVirtualPad *pad,const TObject *object) {
  if(!pad || !object || !pad->GetListOfPrimitives())
    return nullptr;

  if(pad->GetListOfPrimitives()->FindObject(object))
    return pad;

  TIter iter(pad->GetListOfPrimitives());
  while(TObject *obj = iter.Next()) {
    if(!obj->InheritsFrom(TVirtualPad::Class()))
      continue;
    if(TVirtualPad *found = FindPadContaining(static_cast<TVirtualPad*>(obj),object))
      return found;
  }

  return nullptr;
}

TVirtualPad *FindPadByName(TVirtualPad *pad,const char *name) {
  if(!pad || !name || !pad->GetListOfPrimitives())
    return nullptr;

  if(strcmp(pad->GetName(),name) == 0)
    return pad;

  TIter iter(pad->GetListOfPrimitives());
  while(TObject *obj = iter.Next()) {
    if(!obj->InheritsFrom(TVirtualPad::Class()))
      continue;
    auto *child = static_cast<TVirtualPad*>(obj);
    if(strcmp(child->GetName(),name) == 0)
      return child;
    if(TVirtualPad *found = FindPadByName(child,name))
      return found;
  }

  return nullptr;
}

void LockFrameIfPresent(TVirtualPad *pad) {
  if(!pad || !pad->GetListOfPrimitives())
    return;

  TBox *frame = static_cast<TBox*>(pad->GetListOfPrimitives()->FindObject("TFrame"));
  if(frame)
    frame->SetBit(TBox::kCannotMove);
}

} // namespace


GH1D::GH1D() : TH1D(), fOriginal(0) { Init(); }

GH1D::GH1D(std::string name,int nbinsx,const double *xbins) :
  GH1D(name.c_str(),name.c_str(),nbinsx,xbins)  { }

GH1D::GH1D(std::string name,int nbinsx,const float  *xbins) :
  GH1D(name.c_str(),name.c_str(),nbinsx,xbins)  { }

GH1D::GH1D(std::string name,int nbinsx,double xlow, double xhigh) : 
  GH1D(name.c_str(),name.c_str(),nbinsx,xlow,xhigh)  { }

GH1D::GH1D(const char *name,const char *title,int nbinsx,const double *xbins) :
  TH1D(name,title,nbinsx,xbins), fOriginal(0), fParent(0)   { 
  Init();
}
 
GH1D::GH1D(const char *name,const char *title,int nbinsx,const float  *xbins) :
  TH1D(name,title,nbinsx,xbins), fOriginal(0), fParent(0)    { 
  Init();
} 

GH1D::GH1D(const char *name,const char *title,int nbinsx,double xlow, double xhigh) :
  TH1D(name,title,nbinsx,xlow,xhigh), fOriginal(0), fParent(0)     { 
  Init();
}

//GH1D::GH1D(const GH1D &h1d);
GH1D::GH1D(const TH1D &h1d) :
  TH1D(h1d), fOriginal(0), fParent(0)    { 
  Init();
} 

GH1D::GH1D(const TH1F &h1f) :
   fOriginal(0), fParent(0)     {  
  h1f.Copy(*this); 
  Init();
} 

/*
GH1D::GH1D(const TH1 *hist) {
  if(!hist) return;
  if(hist->InheritsFrom(TH1D::Class())) {
    TH1D *h1d = (TH1D*)hist;
    h1d->Copy(*this);
  } else if(hist->InheritsFrom(TH1F::Class())) {
    TH1F *h1f = (TH1F*)hist;
    h1f->Copy(*this);
  } else {
  
  }
  Init();
}
*/

GH1D::GH1D(const TVectorD &v) :
  TH1D(v), fOriginal(0), fParent(0)     { 
  Init();
}

GH1D::~GH1D() { 
  //printf("GH1D deleted\n"); fflush(stdout);  
  if(fOriginal) delete fOriginal;
  if(fBg)       delete fBg;
  if(fResidualHist) delete fResidualHist;

  //TH1D::~TH1D();
} 

void GH1D::Init() { 
  //this->SetBit(kNoTitle);
	SetOriginal();
  fParent = 0;
  fBg = 0;
  fIsNormalized = false;
  fShowResiduals = false;
  fResidualFit = nullptr;
  fResidualHist = nullptr;
  fResidualXLow = std::numeric_limits<double>::quiet_NaN();
  fResidualXHigh = std::numeric_limits<double>::quiet_NaN();
}


//void GH1D::Copy(TObject &newHist) const { 
//  TH1D::Copy(newHist);
//  (*(static_cast<GH1D*>(&newHist))).Init();
//}



void GH1D::Paint(Option_t *opt) {
  TH1D::Paint(opt);

  return;
  //printf("\t-in gh1d paint.\n");
  //fflush(stdout);
  TString sopt(opt);
  bool drawFunctions=false;
  if(!sopt.Length()) {
    if(this->GetSumw2()->GetSize()) {
      sopt.Append("hist");
      drawFunctions=true;
    }
  }
  TH1D::Paint(sopt.Data());
  if(drawFunctions) {
    TIter iter(this->GetListOfFunctions());
    while(TObject *obj = iter.Next()) {
      //if(obj->InheritsFrom(TF1::Class() || TLine::Class())
      if((obj->InheritsFrom("TF1")) || (obj->InheritsFrom("TLine")))
        obj->Paint("lsame");
    }
  }
}

int GH1D::DistancetoPrimitive(int px, int py) {

  /*
  if(gROOT->GetSelectedPad() && gROOT->GetSelectedPad()->GetListOfPrimitives()->FindObject(this)) {
    TObject *frame = gROOT->GetSelectedPad()->GetListOfPrimitives()->FindObject("TFrame");
    if(frame) 
      return frame->DistancetoPrimitive(px,py);
  }
  */
  return TH1D::DistancetoPrimitive(px,py);
}


void GH1D::Print(Option_t *opt) const {
  printf("this: %p\n",this->GetDirectory());
  printf("org:  %p\n",fOriginal->GetDirectory());
  return;
}


void GH1D::Draw(Option_t *opt) {
  //printf("GH1D Draw!\n"); fflush(stdout);
  TH1D::Draw(opt);
  if(gPad) {
    //gPad->Modified(); 
    gPad->Update();
    TBox *frame = (TBox*)gPad->GetListOfPrimitives()->FindObject("TFrame");
    if(frame) {
      frame->SetBit(TBox::kCannotMove); 
      //frame->SetBit(TObject::EStatusBits::kCannotPick);
    }
  }
  return; 
}

void GH1D::SetShowResiduals(bool flag) {
  fShowResiduals = flag;
  if(gPad)
    UpdateResidualDisplay(gPad);
}

void GH1D::ToggleResiduals() {
  SetShowResiduals(!ShowResiduals());
  printf("Residual display %s for %s\n",ShowResiduals() ? "on" : "off",GetName());
}

void GH1D::SetResidualFit(TF1 *fit,double xlow,double xhigh) {
  if(xlow>xhigh)
    std::swap(xlow,xhigh);

  fResidualFit = fit;
  fResidualXLow  = xlow;
  fResidualXHigh = xhigh;

  if(ShowResiduals())
    UpdateResidualDisplay(gPad);
}

void GH1D::ClearResidual() {
  fResidualFit = nullptr;
  if(fResidualHist) {
    delete fResidualHist;
    fResidualHist = nullptr;
  }
  fResidualXLow = std::numeric_limits<double>::quiet_NaN();
  fResidualXHigh = std::numeric_limits<double>::quiet_NaN();
  if(ShowResiduals())
    UpdateResidualDisplay(gPad);
}

TH1D *GH1D::MakeResidualHist() const {
  double viewLow = GetXaxis()->GetBinLowEdge(GetXaxis()->GetFirst());
  double viewHigh = GetXaxis()->GetBinUpEdge(GetXaxis()->GetLast());
  if(viewLow>viewHigh)
    std::swap(viewLow,viewHigh);

  TString name = Form("__%s_residuals",GetName());
  TH1D *residual = nullptr;
  const TArrayD *xbins = GetXaxis()->GetXbins();
  if(xbins && xbins->GetSize() > 0) {
    residual = new TH1D(name,"",
                        GetNbinsX(),xbins->GetArray());
  } else {
    residual = new TH1D(name,"",
                        GetNbinsX(),GetXaxis()->GetXmin(),GetXaxis()->GetXmax());
  }

  residual->SetDirectory(nullptr);
  residual->SetStats(false);
  residual->SetLineColor(kBlue + 2);
  residual->SetMarkerColor(kBlue + 2);
  residual->SetMarkerStyle(kFullCircle);
  residual->SetMarkerSize(0.35);

  double maxAbs = 0.0;
  for(int bin = 1; bin <= GetNbinsX(); ++bin) {
    const double x = GetXaxis()->GetBinCenter(bin);
    if(x < viewLow || x > viewHigh)
      continue;

    double fitLow = fResidualXLow;
    double fitHigh = fResidualXHigh;
    if(fitLow>fitHigh)
      std::swap(fitLow,fitHigh);

    const bool inFitRegion = fResidualFit && std::isfinite(fitLow) && std::isfinite(fitHigh) &&
                             x >= fitLow && x <= fitHigh;
    const double value = inFitRegion ? GetBinContent(bin) - fResidualFit->Eval(x) : 0.0;
    residual->SetBinContent(bin,value);
    residual->SetBinError(bin,GetBinError(bin));
    maxAbs = std::max(maxAbs,std::abs(value));
  }

  if(maxAbs <= 0.0 || !std::isfinite(maxAbs))
    maxAbs = 1.0;
  residual->SetMinimum(-1.15*maxAbs);
  residual->SetMaximum( 1.15*maxAbs);
  residual->GetXaxis()->SetRangeUser(viewLow,viewHigh);
  residual->GetXaxis()->SetTitle("");
  residual->GetYaxis()->SetTitle("");
  residual->GetYaxis()->SetTitleSize(0.10);
  residual->GetYaxis()->SetLabelSize(0.08);
  residual->GetXaxis()->SetTitleSize(0.11);
  residual->GetXaxis()->SetLabelSize(0.10);

  return residual;
}

void GH1D::UpdateResidualDisplay(TVirtualPad *pad) {
  TVirtualPad *rootPad = pad ? pad : gPad;
  if(!rootPad)
    return;
  if(rootPad->GetCanvas())
    rootPad = rootPad->GetCanvas();

  const TString histPadName = ResidualPadName(this,"hist");
  const TString residualPadName = ResidualPadName(this,"residual");

  TVirtualPad *drawingPad = FindPadContaining(rootPad,this);
  TVirtualPad *histPad = FindPadByName(rootPad,histPadName.Data());
  TVirtualPad *residualPad = FindPadByName(rootPad,residualPadName.Data());
  TVirtualPad *container = drawingPad ? drawingPad : pad;

  if(histPad && residualPad) {
    if(auto *tp = dynamic_cast<TPad*>(histPad)) {
      if(tp->GetMother())
        container = tp->GetMother();
    }
  }

  if(!container)
    return;

  const bool haveSplitPads = histPad && residualPad;
  TString drawOption = GetDrawOption();

  if(!ShowResiduals()) {
    if(haveSplitPads) {
      container->cd();
      container->Clear();
      TH1D::Draw(drawOption.Data());
      container->Modified();
      container->Update();
      LockFrameIfPresent(container);
      RequestCurrentPad(container);
    }
    return;
  }

  if(!haveSplitPads) {
    container->cd();
    container->Clear();

    TPad *newHistPad = new TPad(histPadName,histPadName,0.0,0.26,1.0,1.0);
    TPad *newResidualPad = new TPad(residualPadName,residualPadName,0.0,0.0,1.0,0.26);

    newHistPad->SetBottomMargin(0.0);
    newHistPad->SetFillStyle(4000);
    newHistPad->AddExec("groot_interact","GRootInteract()");
    newResidualPad->SetTopMargin(0.0);
    newResidualPad->SetBottomMargin(0.32);
    newResidualPad->SetFillStyle(4000);
    newResidualPad->AddExec("groot_interact","GRootInteract()");

    newHistPad->Draw();
    newResidualPad->Draw();
    histPad = newHistPad;
    residualPad = newResidualPad;
  }

  histPad->cd();
  histPad->Clear();
  TH1D::Draw(drawOption.Data());
  histPad->Modified();
  histPad->Update();
  LockFrameIfPresent(histPad);

  residualPad->cd();
  residualPad->Clear();
  if(fResidualHist) {
    delete fResidualHist;
    fResidualHist = nullptr;
  }
  fResidualHist = MakeResidualHist();
  if(fResidualHist) {
    fResidualHist->Draw("hist");

    double viewLow = GetXaxis()->GetBinLowEdge(GetXaxis()->GetFirst());
    double viewHigh = GetXaxis()->GetBinUpEdge(GetXaxis()->GetLast());
    if(viewLow>viewHigh)
      std::swap(viewLow,viewHigh);

    double fitLow = fResidualXLow;
    double fitHigh = fResidualXHigh;
    if(fitLow>fitHigh)
      std::swap(fitLow,fitHigh);

    if(fResidualFit && std::isfinite(fitLow) && std::isfinite(fitHigh)) {
      const double fillLow = std::max(viewLow,fitLow);
      const double fillHigh = std::min(viewHigh,fitHigh);
      if(fillLow < fillHigh) {
        std::vector<double> x;
        std::vector<double> y;
        x.push_back(fillLow);
        y.push_back(0.0);
        for(int bin = 1; bin <= fResidualHist->GetNbinsX(); ++bin) {
          const double center = fResidualHist->GetXaxis()->GetBinCenter(bin);
          if(center < fillLow || center > fillHigh)
            continue;
          x.push_back(std::max(fillLow,fResidualHist->GetXaxis()->GetBinLowEdge(bin)));
          y.push_back(fResidualHist->GetBinContent(bin));
          x.push_back(std::min(fillHigh,fResidualHist->GetXaxis()->GetBinUpEdge(bin)));
          y.push_back(fResidualHist->GetBinContent(bin));
        }
        x.push_back(fillHigh);
        y.push_back(0.0);
        if(x.size() > 2) {
          TGraph *fill = new TGraph(static_cast<int>(x.size()),x.data(),y.data());
          fill->SetFillColorAlpha(kAzure + 2,0.30);
          fill->SetLineColorAlpha(kAzure + 2,0.0);
          fill->SetBit(TObject::kCanDelete);
          fill->Draw("f same");
          fResidualHist->Draw("hist same");
        }
      }
    }

    TLine *zero = new TLine(viewLow,0.0,viewHigh,0.0);
    zero->SetLineColor(kGray + 2);
    zero->SetLineStyle(kDashed);
    zero->SetBit(TObject::kCanDelete);
    zero->Draw();
  }
  residualPad->Modified();
  residualPad->Update();

  histPad->cd();
  RequestCurrentPad(histPad);
}

TH1* GH1D::DrawCopy(Option_t *opt, const char *name_postfix) const {
  //printf("GH1D Draw!\n"); fflush(stdout);
  return TH1D::DrawCopy(opt,name_postfix);
}

TH1* GH1D::DrawNormalized(Option_t *opt, double norm) const {
  //printf("GH1D Draw!\n"); fflush(stdout);
  return TH1D::DrawNormalized(opt,norm);
} 

void GH1D::SetOriginal()   {
  if(fOriginal) {
    delete fOriginal;
    fOriginal = 0;
  }
	if(!fOriginal) {
  	fOriginal = new TH1D();  
    fOriginalBins = GetNbinsX();
    dynamic_cast<TH1D*>(this)->Copy(*fOriginal);
    fOriginal->SetName(Form("_%s_copy",this->GetName()));
  }
  fOriginal->SetDirectory(nullptr);
  //Print();
}

void GH1D::ResetToOriginal() { 

  if(!fOriginal) return;
  //double xlow = GetXaxis()->GetBinLowEdge(GetXaxis()->GetFirst());
  //double xup  = GetXaxis()->GetBinUpEdge(GetXaxis()->GetLast());
  std::string fname = this->GetName();
  //printf("Unbinng:\n");
  //printf("name: %s\n",fname.c_str());
  //printf("current  bins = %i\n",GetNbinsX());  
  //printf("original bins = %i\n\n",GetNbinsOriginal());  

  TDirectory *current = this->GetDirectory();
  fOriginal->Copy(*(dynamic_cast<TH1D*>(this)));
  this->SetName(fname.c_str());
  this->SetDirectory(current);
  //Print();
} 


TH1* GH1D::Rebin(int ngroup,const char *newname,const double *xbins) {
  TString sname(newname);
  //find the current viewing range
  //double xlow = GetXaxis()->GetBinLowEdge(GetXaxis()->GetFirst());
  //double xup  = GetXaxis()->GetBinUpEdge(GetXaxis()->GetLast());
  /*
  if(sname.Length()==0) {
    // we are not going to return a new histogram
    // so lets copy the current histogram to allow us to unbin latter...
    if(!fOriginal) {
      fOriginal = new TH1D();  
      fOriginal->SetDirectory(0);
      fOriginalBins = GetNbinsX();
      dynamic_cast<TH1D*>(this)->Copy(*fOriginal);
      fOriginal->SetName(Form("_%s_copy",this->GetName()));
      fOriginal->SetDirectory(0);
    }
  }
	*/
  TH1 *temp = TH1D::Rebin(ngroup,newname,xbins);
  //printf("Rebinng:\n");
  //printf("current  bins = %i\n",GetNbinsX());  
  //printf("original bins = %i\n\n",GetNbinsOriginal());  
  //if(sname.Length()==0) 
  //  GetXaxis()->SetRangeUser(xlow,xup);
  UpdateFunctions();
  return temp;
}

void GH1D::Unbin(int ngroup) {
  //printf("unbin called!\n");

  if(!fOriginal) return;
  //double xlow = GetXaxis()->GetBinLowEdge(GetXaxis()->GetFirst());
  //double xup  = GetXaxis()->GetBinUpEdge(GetXaxis()->GetLast());
  //std::string fname = this->GetName();
  //printf("Unbinng:\n");
  //printf("name: %s\n",fname.c_str());
  //printf("current  bins = %i\n",GetNbinsX());  
  //printf("original bins = %i\n\n",GetNbinsOriginal());  

  int oldBins = GetNbinsX();
	
	/*
  TDirectory *current = this->GetDirectory();
  fOriginal->Copy(*(dynamic_cast<TH1D*>(this)));
  this->SetName(fname.c_str());
  //printf("name2 = %s\n",this->GetName());
  this->SetDirectory(current);
  //printf("name3 = %s\n",this->GetName());
	*/

	ResetToOriginal();

  if(ngroup>0) {
    int aim = oldBins*ngroup;
    //printf("aim = %i\n",aim);
    if(aim<GetNbinsOriginal()) {
      while(GetNbinsX()>(aim+2)) {
        //printf("\tcurrent = %i\n",GetNbinsX());
        Rebin(ngroup);
      }
    }
  }
  UpdateFunctions();
}

/*
void GH1D::SetSubtract(TH1D *h, double scale) {
	
  if(!fOriginal) {
    fOriginal = new TH1D();  
    fOriginal->SetDirectory(0);
    fOriginalBins = GetNbinsX();
    dynamic_cast<TH1D*>(this)->Copy(*fOriginal);
    fOriginal->SetName(Form("_%s_copy",this->GetName()));
    fOriginal->SetDirectory(0);
  }


	if(fSubtract) {
		printf("DELETING SUBTRACT\n");
		fflush(stdout);
		delete fSubtract;
	}
	fSubtract = new TH1D();
	fSubtract->SetDirectory(0);
	h->Copy(*fSubtract);
	fSubtract->SetName(Form("%s_%s",this->GetName(),h->GetName()));
	fSubtract->SetDirectory(0);

	SetScale(scale);

	DoSubtract();
}

void GH1D::DoSubtract() {
	if(!fOriginal || !fSubtract)
		return;
	
  std::string fname = this->GetName();
  TDirectory *current = this->GetDirectory();
  fOriginal->Copy(*(dynamic_cast<TH1D*>(this)));
  this->SetName(fname.c_str());
  //printf("name2 = %s\n",this->GetName());
  this->SetDirectory(current);
  //printf("name3 = %s\n",this->GetName());

	this->Add(fSubtract,fScale);

}
*/


void GH1D::SetBackground(int niter,Option_t* opt) {
  if(!fBg) {
    fBg = static_cast<TH1D*>(TSpectrum::StaticBackground(this,niter,opt));
    if(fBg) fBg->SetDirectory(nullptr);
  }
  if(fBg->GetNbinsX() != this->GetNbinsX()) {
    delete fBg;
    fBg = static_cast<TH1D*>(TSpectrum::StaticBackground(this,niter,opt));
    if(fBg) fBg->SetDirectory(nullptr);
  }
}

void GH1D::ShowBackground() {
  SetBackground();
  if(gPad) {
    if(gPad->GetListOfPrimitives()->FindObject(fBg)) // the bg is already drawn...
      gPad->GetListOfPrimitives()->Remove(fBg);
    else
      fBg->Draw("same");
  }
}

void GH1D::ToggleBackground() {
  double x1=sqrt(-1);
  double x2=sqrt(-1);
  if(this->GetXaxis()) {
    x1 = this->GetXaxis()->GetBinLowEdge(this->GetXaxis()->GetFirst());
    x2 = this->GetXaxis()->GetBinUpEdge(this->GetXaxis()->GetLast());
  }

  if(this->TestBit(GH1D::kBackgroundRemoved)) {
    //put the background back...
    std::string fname = this->GetName();
    TDirectory *current = this->GetDirectory();
    fOriginal->Copy(*(dynamic_cast<TH1D*>(this)));
    this->SetName(fname.c_str());
    //printf("name2 = %s\n",this->GetName());
    this->SetDirectory(current);
    this->SetBit(GH1D::kBackgroundRemoved,0);

  } else {
    //remove the background...
    if(!fOriginal) {
      fOriginal = new TH1D();  
      fOriginal->SetDirectory(0);
      fOriginalBins = GetNbinsX();
      dynamic_cast<TH1D*>(this)->Copy(*fOriginal);
      fOriginal->SetName(Form("_%s_copy",this->GetName()));
      fOriginal->SetDirectory(0);
    }
    this->GetXaxis()->UnZoom();
    //TH1 *bg = TSpectrum::StaticBackground(this,12,"COMPTON");
    SetBackground();
    TH1D::Add(fBg,-1); 
    this->SetBit(GH1D::kBackgroundRemoved);
  }
  if(x1==x1) {
    //printf("%s I AM HERE!\n",GetName());
    this->GetXaxis()->SetRangeUser(x1,x2);
  }
}


bool GH1D::Add(const TH1 *h1, Double_t c1) {
  ResetToOriginal();
  bool result = TH1::Add(h1,c1);
  SetOriginal();
  return result;
}
/*
GH1D *GH1D::AddNormalized(const TH1 *h1, Double_t c1) {
  ResetToOriginal();
  GH1D *temp = new GH1D(this);
  if(!h1->InheritsFrom(GH1D::Class())) {
    h1 = new GH1D(h1);
  }
  printf("here 1\n"); fflush(stdout);
  temp->Normalize();
  printf("here 2\n"); fflush(stdout);
  ((GH1D*)h1)->Normalize();
  printf("here 3\n"); fflush(stdout);
  temp->Add(h1,c1);
  printf("here 4\n"); fflush(stdout);
  delete h1;
  return temp;
}
*/






void GH1D::Normalize() {
  
  if(!fIsNormalized) {
    double sum = GetSumOfWeights();
    if(sum==0) {
      printf("GH1D::Normalize sum of weights is zero\n");
      return;
    }
    double max = GetMaximum();
    double min = GetMinimum();
    this->Scale(1.0/sum);
    if (TMath::Abs(max+1111) > 1e-3) SetMaximum(max*1.0/sum);
    if (TMath::Abs(min+1111) > 1e-3) SetMinimum(min*1.0/sum);
    Sumw2(false);  
    fIsNormalized = true;
  } else if(fIsNormalized) {
    std::string fname = this->GetName();
    //fOriginal->Copy(*(dynamic_cast<TH1D*>(this)));
    //this->SetName(fname.c_str());
    ResetToOriginal();
    fIsNormalized = false;
  }

}

void GH1D::UpdateFunctions() {
  TIter iter(this->GetListOfFunctions());
  while(TObject *obj = iter.Next()) {
    if(obj->InheritsFrom(GROI::Class())) {
      static_cast<GROI*>(obj)->Update();
    }
  }

}


bool GH1D::RemovePeaks() {
  bool flag = false;
  if(TObject *obj = GetListOfFunctions()->FindObject("PeakLabels")) {
    GetListOfFunctions()->Remove(obj);
    ((TObjArray*)obj)->Delete();
    flag = true;
  }
  return flag;
}



int GH1D::ShowPeaks(double thresh,double sigma) {
  RemovePeaks();

  TSpectrum::StaticSearch(this,sigma,"Qnodraw",thresh);
  TPolyMarker *pm = (TPolyMarker*)GetListOfFunctions()->FindObject("TPolyMarker");
  if(!pm) {
    //something has gone wrong....
    return 0;
  }
  TObjArray *array = new TObjArray();
  array->SetName("PeakLabels");
  int n = pm->GetN();
  if(n==0)
    return n;
  TText *text;
  double *x = pm->GetX();
  //  double *y = pm->GetY();
  for(int i=0;i<n;i++) {
    //y[i] += y[i]*0.15;
    double y = 0;
    for(int i_x = x[i]-3;i_x<x[i]+3;i_x++){
      if((GetBinContent(GetXaxis()->FindBin(i_x)))>y){
	      y = GetBinContent(GetXaxis()->FindBin(i_x));
      }
    }
    y+=y*0.1;
    text = new TText(x[i],y,Form("%.1f",x[i]));
    text->SetTextSize(0.025);
    text->SetTextAngle(90);
    text->SetTextAlign(12);
    text->SetTextFont(42);
    text->SetTextColor(GetLineColor());
    array->Add(text);
  }
  GetListOfFunctions()->Remove(pm);
  pm->Delete();
  GetListOfFunctions()->Add(array);
  return n;
}
