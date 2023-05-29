TCutG* get_cut_dEdx_m2_proton(){
//========= Macro generated from object:"cut_dEdx_m2_proton/Graph
//========= by ROOT version6.26/06
   
   auto cutg = new TCutG("cut_dEdx_m2_proton",13);
   cutg->SetVarX("trdEdx");
   cutg->SetVarY("trTofM2");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,0.923524,1.48188);
   cutg->SetPoint(1,0.847547,0.893508);
   cutg->SetPoint(2,0.928426,0.474399);
   cutg->SetPoint(3,1.16371,0.480981);
   cutg->SetPoint(4,1.19312,0.695252);
   cutg->SetPoint(5,1.9823,0.774128);
   cutg->SetPoint(6,3.32294,0.767733);
   cutg->SetPoint(7,3.38666,1.00436);
   cutg->SetPoint(8,3.25186,1.09176);
   cutg->SetPoint(9,1.34753,1.0939);
   cutg->SetPoint(10,1.20292,1.4968);
   cutg->SetPoint(11,0.928426,1.50107);
   cutg->SetPoint(12,0.923524,1.48188);
   return cutg;
}
