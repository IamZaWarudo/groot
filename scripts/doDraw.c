

void doDraw(TH2D* mat) {

  GH2D gmat(*mat);
  
  gmat->RebinX(2);

  GH1D *p1 = gmat->ProjectionX(192,192);
  GH1D *p1 = gmat->ProjectionX(193,193);
  GH1D *p1 = gmat->ProjectionX(194,194);
  GH1D *p1 = gmat->ProjectionX(195,195);
  GH1D *p1 = gmat->ProjectionX(196,196);
  GH1D *p1 = gmat->ProjectionX(197,197);
  GH1D *p1 = gmat->ProjectionX(198,198);
  GH1D *p1 = gmat->ProjectionX(199,199);
  GH1D *p1 = gmat->ProjectionX(201,201);
  GH1D *p1 = gmat->ProjectionX(202,202);
  GH1D *p1 = gmat->ProjectionX(203,203);
  GH1D *p1 = gmat->ProjectionX(204,204);
  GH1D *p1 = gmat->ProjectionX(205,205);
  GH1D *p1 = gmat->ProjectionX(206,206);
}
