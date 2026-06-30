
#include<GMarker.h>

#include<algorithm>
#include<iostream>

#include<TH1.h>
#include<TList.h>
#include<TVirtualPad.h>
#include<TPad.h>
//#include<TTimer.h>
#include<TEnv.h>
#include<TCutG.h>
#include<TMarker.h>
#include<TGraph.h>

#include<vector>

namespace {

const char* MarkerTypeName(GMarkerType type) {
  switch(type) {
    case GMarkerType::kPrimary:
      return "Primary";
    case GMarkerType::kBackground:
      return "Background";
    case GMarkerType::kZoom:
      return "Zoom";
    case GMarkerType::kFit:
      return "Fit";
    case GMarkerType::kCut:
      return "Cut";
    case GMarkerType::kProjection:
      return "Projection";
    case GMarkerType::kAll:
    default:
      return "Primary";
  }
}

Color_t DefaultMarkerColor(GMarkerType type) {
  switch(type) {
    case GMarkerType::kBackground:
      return kAzure + 2;
    case GMarkerType::kZoom:
      return kGreen + 2;
    case GMarkerType::kFit:
      return kViolet + 1;
    case GMarkerType::kCut:
      return kOrange + 7;
    case GMarkerType::kProjection:
      return kMagenta + 1;
    case GMarkerType::kAll:
    case GMarkerType::kPrimary:
    default:
      return kRed;
  }
}

Style_t DefaultLineStyle(GMarkerType type) {
  switch(type) {
    case GMarkerType::kBackground:
    case GMarkerType::kZoom:
      return kDashed;
    case GMarkerType::kProjection:
      return kDotted;
    case GMarkerType::kAll:
    case GMarkerType::kPrimary:
    case GMarkerType::kFit:
    case GMarkerType::kCut:
    default:
      return kSolid;
  }
}

double DefaultMarkerSize(GMarkerType type) {
  return type == GMarkerType::kCut ? 0.65 : 1.1;
}

int MarkerEnvInt(GMarkerType type,const char* field,int defaultValue) {
  return gEnv->GetValue(Form("GMarker.%s.%s",MarkerTypeName(type),field),defaultValue);
}

double MarkerEnvDouble(GMarkerType type,const char* field,double defaultValue) {
  return gEnv->GetValue(Form("GMarker.%s.%s",MarkerTypeName(type),field),defaultValue);
}

} // namespace

int GMarker::GetMaxMarkers(GMarkerType type) { 
  switch(type) {
    case GMarkerType::kPrimary:
      return gEnv->GetValue("GMarker.Primary.Max",2);
    case GMarkerType::kBackground:
      return gEnv->GetValue("GMarker.Background.Max",2);
    case GMarkerType::kZoom:
      return gEnv->GetValue("GMarker.Zoom.Max",2);
    case GMarkerType::kFit:
      return gEnv->GetValue("GMarker.Fit.Max",-1);
    case GMarkerType::kCut:
      // Polygon-gate vertices should not silently displace earlier vertices.
      return gEnv->GetValue("GMarker.Cut.Max",-1);
    case GMarkerType::kProjection:
      return gEnv->GetValue("GMarker.Projection.Max",2);
    default: 
      return 2;
  };
  return 2;
}

void GMarker::SetMaxMarkers(GMarkerType type,int value) { 
  switch(type) {
    case GMarkerType::kPrimary:
      gEnv->SetValue("GMarker.Primary.Max",value);
      break;
    case GMarkerType::kBackground:
      gEnv->SetValue("GMarker.Background.Max",value);
      break;
    case GMarkerType::kZoom:
      gEnv->SetValue("GMarker.Zoom.Max",value);
      break;
    case GMarkerType::kFit:
      gEnv->SetValue("GMarker.Fit.Max",value);
      break;
    case GMarkerType::kCut:
      gEnv->SetValue("GMarker.Cut.Max",value);
      break;
    case GMarkerType::kProjection:
      gEnv->SetValue("GMarker.Projection.Max",value);
      break;
    default: 
      return;
  };
  return;
}

GMarker::GMarker() :fHist(0), fLineX(0), fLineY(0), fPoint(0), fX(sqrt(-1)), fY(sqrt(-1)) {
  //SetLineWidth(2);
  //SetLineColor(kRed);
  //SetName("GMarker");


}

GMarker::~GMarker() {
  //printf("gmarker deleted\n");
  if(fLineX) delete fLineX;
  if(fLineY) delete fLineY;
  if(fPoint) delete fPoint;
}

void GMarker::SetType(GMarkerType type) {
  fType = type;
  if(fType == GMarkerType::kCut && fHist && fHist->GetDimension() == 2 && !fPoint)
    fPoint = new TMarker;
  UpdateStyle();
}

void GMarker::UpdateStyle() {
  Color_t color = MarkerEnvInt(fType,"Color",DefaultMarkerColor(fType));
  Style_t style = MarkerEnvInt(fType,"Style",DefaultLineStyle(fType));
  Width_t width = MarkerEnvInt(fType,"Width",2);
  Style_t markerStyle = MarkerEnvInt(fType,"MarkerStyle",kFullCircle);
  Size_t markerSize = MarkerEnvDouble(fType,"MarkerSize",DefaultMarkerSize(fType));

  if(fLineX) {
    fLineX->SetLineColor(color);
    fLineX->SetLineStyle(style);
    fLineX->SetLineWidth(width);
  }
  if(fLineY) {
    fLineY->SetLineColor(color);
    fLineY->SetLineStyle(style);
    fLineY->SetLineWidth(width);
  }
  if(fPoint) {
    fPoint->SetMarkerColor(color);
    fPoint->SetMarkerStyle(markerStyle);
    fPoint->SetMarkerSize(markerSize);
  }
}

void GMarker::PaintBackgroundRegion() const {
  if(!fHist || fHist->GetDimension() != 1 || fType != GMarkerType::kBackground || !gPad)
    return;

  const auto markers = Get(fHist,GMarkerType::kBackground);
  // A background is one interval. Draw it once, from the first marker only.
  if(markers.size() != 2 || markers.front() != this)
    return;

  double low = markers[0]->X();
  double high = markers[1]->X();
  if(low > high) std::swap(low,high);

  const int firstBin = fHist->GetXaxis()->FindBin(low);
  const int lastBin = fHist->GetXaxis()->FindBin(high);
  if(firstBin > lastBin)
    return;

  const double baseline = gPad->GetUymin();
  std::vector<double> x;
  std::vector<double> y;
  x.reserve(2*(lastBin-firstBin+1)+3);
  y.reserve(2*(lastBin-firstBin+1)+3);

  x.push_back(fHist->GetXaxis()->GetBinLowEdge(firstBin));
  y.push_back(baseline);
  for(int bin = firstBin; bin <= lastBin; ++bin) {
    const double content = fHist->GetBinContent(bin);
    x.push_back(fHist->GetXaxis()->GetBinLowEdge(bin));
    y.push_back(content);
    x.push_back(fHist->GetXaxis()->GetBinUpEdge(bin));
    y.push_back(content);
  }
  x.push_back(fHist->GetXaxis()->GetBinUpEdge(lastBin));
  y.push_back(baseline);

  TGraph region(static_cast<int>(x.size()),x.data(),y.data());
  const Color_t fillColor = MarkerEnvInt(GMarkerType::kBackground,"Color",
                                         DefaultMarkerColor(GMarkerType::kBackground));
  const double fillAlpha = MarkerEnvDouble(GMarkerType::kBackground,"FillAlpha",0.25);
  region.SetFillColorAlpha(fillColor,fillAlpha);
  region.SetLineColorAlpha(fillColor,0.0);
  region.Paint("F");
}

void GMarker::AddTo(TH1 *h, double x, double y,bool ignoreMax,Option_t *opt) {
  if(!h) return;
  fHist = h;
  //if(h && h->GetDimension() == 1) {
  x = fHist->GetXaxis()->GetBinLowEdge(fHist->GetXaxis()->FindBin(x));
  SetX(x);

  if(!fLineX) fLineX = new TLine;

  if(h->GetDimension() == 2) {
    y = fHist->GetYaxis()->GetBinLowEdge(fHist->GetYaxis()->FindBin(y));
    SetY(y);

    if(!fLineY) fLineY = new TLine;
    if(fType == GMarkerType::kCut && !fPoint)
      fPoint = new TMarker;
  }

  UpdateStyle();

  if(!ignoreMax) {
    const int maxMarkers = GetMaxMarkers(fType);
    auto markers = Get(fHist,fType);

    while(maxMarkers >= 0 && static_cast<int>(markers.size()) >= maxMarkers) {
      markers.front()->Remove();
      markers.erase(markers.begin());
    }
  }
  fHist->GetListOfFunctions()->Add(this);
} 

void GMarker::Remove() {
  TH1* hist = fHist;  
  fHist = nullptr;
  if(hist && hist->GetListOfFunctions()) 
    hist->GetListOfFunctions()->Remove(this);
  delete this;

}


void GMarker::Paint(Option_t *opt) {
  TString sopt(opt);
  sopt.ToLower();
  //if(sopt.Length() == 0) sopt = this->GetDrawOption();

  sopt.Append("ndc");

  if(!gPad || !fHist) return;


  double lm = gPad->GetLeftMargin();
  double rm = 1.-gPad->GetRightMargin();
  double tm = 1.-gPad->GetTopMargin();
  double bm = gPad->GetBottomMargin();

  double xndc  = (rm-lm)*((fX-gPad->GetUxmin())/(abs(gPad->GetUxmax()-gPad->GetUxmin())))+lm;
  double yndc  = (tm-bm)*((fY-gPad->GetUymin())/(abs(gPad->GetUymax()-gPad->GetUymin())))+bm;

  /*
     printf("lm:  %.04f\n",lm);
     printf("rm:  %.04f\n",rm);
     printf("tm:  %.04f\n",tm);
     printf("bm:  %.04f\n",bm);

     printf("fX:  %.04f\n",fX);
     printf("fY:  %.04f\n",fY);
     printf("Uxmin: %.04f\n",gPad->GetUxmin());
     printf("Uxmax: %.04f\n",gPad->GetUxmax());
     printf("Uymin: %.04f\n",gPad->GetUymin());
     printf("Uymax: %.04f\n",gPad->GetUymax());



     printf("xndc:  %.04f\n",xndc);
     printf("yndc:  %.04f\n",yndc);
   */



  if(fType == GMarkerType::kCut && fHist->GetDimension() == 2) {
    if(fPoint) {
      fPoint->SetX(fX);
      fPoint->SetY(fY);
      fPoint->Paint();
    }

    auto markers = Get(fHist,GMarkerType::kCut);
    auto found = std::find(markers.begin(),markers.end(),this);
    if(found != markers.end() && markers.size() > 1) {
      const size_t index = static_cast<size_t>(std::distance(markers.begin(),found));
      if(index > 0) {
        TLine segment(markers[index-1]->X(),markers[index-1]->Y(),fX,fY);
        segment.SetLineColor(MarkerEnvInt(GMarkerType::kCut,"Color",
                                          DefaultMarkerColor(GMarkerType::kCut)));
        segment.SetLineWidth(MarkerEnvInt(GMarkerType::kCut,"Width",2));
        segment.SetLineStyle(MarkerEnvInt(GMarkerType::kCut,"Style",
                                          DefaultLineStyle(GMarkerType::kCut)));
        segment.Paint();
      }
      if(index == markers.size()-1 && markers.size() > 2) {
        TLine closing(fX,fY,markers.front()->X(),markers.front()->Y());
        closing.SetLineColor(MarkerEnvInt(GMarkerType::kCut,"Color",
                                          DefaultMarkerColor(GMarkerType::kCut)));
        closing.SetLineWidth(MarkerEnvInt(GMarkerType::kCut,"Width",2));
        closing.SetLineStyle(MarkerEnvInt(GMarkerType::kCut,"ClosingStyle",kDashed));
        closing.Paint();
      }
    }
    return;
  }

  PaintBackgroundRegion();

  if(fLineX) {
    if(!fLineX->TestBit(TLine::kLineNDC))
      fLineX->SetBit(TLine::kLineNDC,true);

    fLineX->SetX1(xndc);
    fLineX->SetX2(xndc);
    /*
       double xuser = ((fLineX->GetX1()-lm)/(rm-lm))*(gPad->GetUxmax()-gPad->GetUxmin())+gPad->GetUxmin();
       if(fX !=xuser) {
       fX= xuser;
       }
     */
    fLineX->SetY1(bm);    
    fLineX->SetY2(tm);    
    fLineX->SetBit(TLine::kLineNDC,true);
  }
  if(fLineY) {
    if(!fLineY->TestBit(TLine::kLineNDC))
      fLineY->SetBit(TLine::kLineNDC,true);

    fLineY->SetY1(yndc);
    fLineY->SetY2(yndc);
    /*
       double yuser = ((fLineY->GetY1()-tm)/(bm-tm))*(gPad->GetUymax()-gPad->GetUymin())+gPad->GetUymin();
       if(fY !=yuser) {
       fY= yuser;
       }
     */
    fLineY->SetX1(lm);    
    fLineY->SetX2(rm);    
    fLineY->SetBit(TLine::kLineNDC,true);
  }

  if(fLineX) {  fLineX->Paint(sopt.Data()); }
  if(fLineY) {  fLineY->Paint(sopt.Data()); }

}

void GMarker::RemoveAll(TH1 *h,bool removeBGMarkers) {
  //remove all markers from h
  TIter iter(h->GetListOfFunctions(),kIterBackward);
  while(TObject *obj = iter.Next()) {
    if(obj->InheritsFrom(GMarker::Class())) {
      GMarker *marker = ((GMarker*)obj);
      if(marker->GetType()!=GMarkerType::kBackground || removeBGMarkers)
        marker->Remove();
    }
  }
}


std::vector<GMarker*> GMarker::Get(TH1 *h,GMarkerType type) {
  //type:
  //  - 0: all
  //  - 1: primary
  //  - 2: bg
  //return all markers in h
  std::vector<GMarker*> toReturn;
  if(!h) return toReturn;

  TIter iter(h->GetListOfFunctions());
  while(TObject *obj = iter.Next()) {
    if(!obj->InheritsFrom(GMarker::Class())) continue;
    auto *marker = static_cast<GMarker*>(obj);
    if(type == GMarkerType::kAll || marker->GetType() == type)
      toReturn.push_back(marker);

  }
  return toReturn;
}

std::vector<GMarker*> GMarker::GetBG(TH1 *h) {
  return Get(h,GMarkerType::kBackground);
}

void GMarker::SetLineColor(Color_t color) {
  if(fLineX) fLineX->SetLineColor(color);
  if(fLineY) fLineY->SetLineColor(color);
  if(fPoint) fPoint->SetMarkerColor(color);
}

void GMarker::ExecuteEvent(int event, int px, int py) { 
  if(!gPad || !fHist) return;

  switch(event) {
    case kButton1Down:
      SetBit(kCannotPick, false);
      break;

    case kButton1Motion:
    case kButton1Up: {
      double x = gPad->PadtoX(gPad->AbsPixeltoX(px));
      x = fHist->GetXaxis()->GetBinLowEdge(fHist->GetXaxis()->FindBin(x));

      SetX(x);
      if(fHist->GetDimension() == 2) {
        double y = gPad->PadtoY(gPad->AbsPixeltoY(py));
        y = fHist->GetYaxis()->GetBinLowEdge(fHist->GetYaxis()->FindBin(y));
        SetY(y);
      }

      gPad->Modified();
      gPad->Update();
      break;
    }

    default:
      break;
  }
 return;









 //printf("GMarker::ExecuteEvent(%i,%i,%i)\n",event,px,py);
  if(!fLineX || fLineY) return;
  //printf("here\n");

  //if(fLineX)
  //  d1 = fLineX->DistancetoPrimitive(px,py);
  //if(fLineY)
  //  d2 = fLineY->DistancetoPrimitive(px,py);
  //printf("d1 = %i, d2 = %i\n",d1,d2);
  //if(d1<d2) { 
  if(fLineX) {
    //printf("here!"); 
    fLineX->ExecuteEvent(event,px,py);
  }
  // } else {
  //   if(fLineY)
  //     fLineY->ExecuteEvent(event,px,py);
  // } 

  return;
} 


//TODO: Currently events are not being sent to the lines...
int GMarker::DistancetoPrimitive(int px, int py) { 
  int d1 = 9999;
  int d2 = 9999;
  if(fType == GMarkerType::kCut && fPoint)
    d1 = fPoint->DistancetoPrimitive(px,py);
  if(fLineX)
    d1 = std::min(d1,fLineX->DistancetoPrimitive(px,py));
  //if(fLineY)
  //  d2 = fLineY->DistancetoPrimitive(px,py);
  //return (d1 < d2) ? d1 : d2;
  return d1;
}


/****************
Current issue:
- Primary markers in 2D draw as full vertical/horizontal lines.
- This becomes visually cluttered while building polygon gates.

Desired direction:
- Marker appearance should depend on GMarkerType and/or interaction mode.
- Keep current functionality for 1D range markers.
- Improve temporary visualization for 2D gate construction.

Ideas:
- Different drawing styles by marker type:
    kPrimary      -> current crosshair
    kCut          -> point/small cross
    kProjection   -> line
    kZoom         -> dashed line
    kBackground   -> alternate color/style

- Possibly add:
    GMarker::UpdateStyle()
    GMarker::PaintAsPoint()
    GMarker::PaintAsCross()
    GMarker::PaintAsLine()

- While constructing polygon cuts:
    - draw connecting polyline between markers
    - avoid full-screen crosshair clutter

- Consider:
    marker rendering strategy depending on histogram dimension
    and current interaction context.
****************/




TCutG* GMarker::MakeTCutG(TH1* h,GMarkerType type) {
  if(!h || h->GetDimension() != 2) 
    return nullptr;

  auto markers = GMarker::Get(h,type);
  std::vector<double> x;
  std::vector<double> y;


  if(markers.size() < 2) {
    return nullptr;
  } else if(markers.size()==2) {
    double x1 = markers[0]->X();
    double y1 = markers[0]->Y();
    double x2 = markers[1]->X();
    double y2 = markers[1]->Y();

    if(x1 > x2) std::swap(x1, x2);
    if(y1 > y2) std::swap(y1, y2);

    //midpoints
    const double xm = 0.5 * (x1 + x2);
    const double ym = 0.5 * (y1 + y2);
    x = {x1, xm, x2, x2, x2, xm, x1, x1, x1};
    y = {y2, y2, y2, ym, y1, y1, y1, ym, y2};
  } else { // more than 2
    const int nMarkers = static_cast<int>(markers.size());

    for(int i = 0; i < nMarkers; ++i) {
      const int next = (i + 1) % nMarkers;
      const double x1 = markers[i]->X();
      const double y1 = markers[i]->Y();
      const double x2 = markers[next]->X();
      const double y2 = markers[next]->Y();

      x.push_back(x1);
      y.push_back(y1);

      //midpoints
      x.push_back(0.5 * (x1 + x2));
      y.push_back(0.5 * (y1 + y2));
    }
    x.push_back(x.front());
    y.push_back(y.front());
  }
  TCutG *cut = new TCutG("CUT",static_cast<int>(x.size()),x.data(),y.data());
  cut->SetLineWidth(MarkerEnvInt(GMarkerType::kCut,"Width",2));
  cut->SetLineColor(MarkerEnvInt(GMarkerType::kCut,"Color",
                                 DefaultMarkerColor(GMarkerType::kCut)));
  cut->SetLineStyle(MarkerEnvInt(GMarkerType::kCut,"Style",
                                 DefaultLineStyle(GMarkerType::kCut)));
  return cut;
}
