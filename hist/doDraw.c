

void doDraw() {

TCanvas *c = new TCanvas; c->Divide(2,2); c->cd(1); Forced_205->Draw(); c->cd(2); Forced_204->Draw(); c->cd(3); Forced_203->Draw(); c->cd(4); Forced_202->Draw();

TH2D *qdc = (TH2D*)_file0->Get("QDCsummary");


TH1D *p1 = qdc->ProjectionX("p1",204,204);
TH1D *p21 = qdc->ProjectionX("p1",204,204);
TH1D *p13 = qdc->ProjectionX("p1",204,204);
TH1D *p31 = qdc->ProjectionX("p1",204,204);


}
