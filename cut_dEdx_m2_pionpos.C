TCutG* get_cut_dEdx_m2_pionpos(){
//========= Macro generated from object:"cut_dEdx_m2_pionpos/Graph
//========= by ROOT version6.26/06
   
   auto cutg = new TCutG("cut_dEdx_m2_pionpos",16);
   cutg->SetVarX("trdEdx");
   cutg->SetVarY("trTofM2");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,0.835319,0.0308953);
   cutg->SetPoint(1,0.920345,0.0617093);
   cutg->SetPoint(2,1.01607,0.098686);
   cutg->SetPoint(3,1.12548,0.154151);
   cutg->SetPoint(4,1.19564,0.194826);
   cutg->SetPoint(5,1.29256,0.208384);
   cutg->SetPoint(6,1.34548,0.130733);
   cutg->SetPoint(7,1.36213,0.0210349);
   cutg->SetPoint(8,1.33299,-0.071407);
   cutg->SetPoint(9,1.26283,-0.142895);
   cutg->SetPoint(10,1.16175,-0.13057);
   cutg->SetPoint(11,0.987534,-0.0270349);
   cutg->SetPoint(12,0.887048,-0.00854651);
   cutg->SetPoint(13,0.837103,0.00254651);
   cutg->SetPoint(14,0.837103,0.00254651);
   cutg->SetPoint(15,0.835319,0.0308953);
   return cutg;
}
