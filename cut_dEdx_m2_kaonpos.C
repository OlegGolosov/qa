TCutG* get_cut_dEdx_m2_kaonpos(){
//========= Macro generated from object:"cut_dEdx_m2_kaonpos/Graph
//========= by ROOT version6.26/06
   
   auto cutg = new TCutG("cut_dEdx_m2_kaonpos",15);
   cutg->SetVarX("trdEdx");
   cutg->SetVarY("trTofM2");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,0.968998,0.380221);
   cutg->SetPoint(1,1.03574,0.419074);
   cutg->SetPoint(2,1.0923,0.410748);
   cutg->SetPoint(3,1.1251,0.371895);
   cutg->SetPoint(4,1.15338,0.29974);
   cutg->SetPoint(5,1.16469,0.234523);
   cutg->SetPoint(6,1.14885,0.185957);
   cutg->SetPoint(7,1.07533,0.148492);
   cutg->SetPoint(8,1.01877,0.144329);
   cutg->SetPoint(9,0.970129,0.154043);
   cutg->SetPoint(10,0.936193,0.170694);
   cutg->SetPoint(11,0.912439,0.256725);
   cutg->SetPoint(12,0.912439,0.315004);
   cutg->SetPoint(13,0.953161,0.374671);
   cutg->SetPoint(14,0.968998,0.380221);
   return cutg;
}
