TCutG* get_cut_dEdx_m2_kaonneg(){
//========= Macro generated from object:"cut_dEdx_m2_kaonneg/Graph
//========= by ROOT version6.26/06
   
   auto cutg = new TCutG("cut_dEdx_m2_kaonneg",15);
   cutg->SetVarX("trdEdx");
   cutg->SetVarY("trTofM2");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,0.943642,0.414818);
   cutg->SetPoint(1,1.00977,0.452957);
   cutg->SetPoint(2,1.06708,0.445012);
   cutg->SetPoint(3,1.10015,0.405283);
   cutg->SetPoint(4,1.12826,0.333771);
   cutg->SetPoint(5,1.13928,0.268616);
   cutg->SetPoint(6,1.12385,0.219353);
   cutg->SetPoint(7,1.05772,0.152609);
   cutg->SetPoint(8,0.994893,0.125593);
   cutg->SetPoint(9,0.948602,0.151019);
   cutg->SetPoint(10,0.910577,0.20505);
   cutg->SetPoint(11,0.88688,0.290864);
   cutg->SetPoint(12,0.88688,0.349663);
   cutg->SetPoint(13,0.928211,0.408461);
   cutg->SetPoint(14,0.943642,0.414818);
   return cutg;
}
