TCutG* get_cut_dEdx_m2_pionneg(){
//========= Macro generated from object:"cut_dEdx_m2_pionneg/Graph
//========= by ROOT version6.26/06
   
   auto cutg = new TCutG("cut_dEdx_m2_pionneg",16);
   cutg->SetVarX("trdEdx");
   cutg->SetVarY("trTofM2");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,0.83303,0.0318915);
   cutg->SetPoint(1,0.918434,0.0624031);
   cutg->SetPoint(2,1.01289,0.110713);
   cutg->SetPoint(3,1.11243,0.202248);
   cutg->SetPoint(4,1.19105,0.231488);
   cutg->SetPoint(5,1.29059,0.208605);
   cutg->SetPoint(6,1.34319,0.131054);
   cutg->SetPoint(7,1.36016,0.0217209);
   cutg->SetPoint(8,1.33075,-0.0710853);
   cutg->SetPoint(9,1.26061,-0.142279);
   cutg->SetPoint(10,1.15937,-0.129566);
   cutg->SetPoint(11,0.985173,-0.0265891);
   cutg->SetPoint(12,0.883933,0.00265116);
   cutg->SetPoint(13,0.834727,0.00265116);
   cutg->SetPoint(14,0.834727,0.00265116);
   cutg->SetPoint(15,0.83303,0.0318915);
   return cutg;
}
