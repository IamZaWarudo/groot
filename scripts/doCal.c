





void doCal(double p1,double p2,double p3, double p4, double p5) {

 //1332
 //1173
 // 356 
 // 302
 //  80

 double e[5] = {80,302,356,1173,1332};
 double c[5] = {p1,p2,p3,p4,p5}; 
 //double c[5] = 


 TGraph *gre = new TGraph(5,c,e);

 TVirtualPad *current = gPad;

 new GCanvas;

 gre->Draw("A*"); 


 gre->Fit("pol1");

 if(current) 
   current->cd();

}




