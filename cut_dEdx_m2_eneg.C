TCutG* get_cut_dEdx_m2_eneg(){
//========= Macro generated from object:"cut_dEdx_m2_eneg/Graph
//========= by ROOT version6.26/06
   
   auto cutg = new TCutG("cut_dEdx_m2_eneg",11);
   cutg->SetVarX("trdEdx");
   cutg->SetVarY("trTofM2");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,1.40111,0.00581395);
   cutg->SetPoint(1,1.41355,-0.0406977);
   cutg->SetPoint(2,1.52637,-0.0833333);
   cutg->SetPoint(3,1.65085,-0.0988372);
   cutg->SetPoint(4,1.78571,-0.0639535);
   cutg->SetPoint(5,1.82363,0.0639535);
   cutg->SetPoint(6,1.69271,0.153101);
   cutg->SetPoint(7,1.57251,0.141473);
   cutg->SetPoint(8,1.45446,0.114341);
   cutg->SetPoint(9,1.40789,0.0523256);
   cutg->SetPoint(10,1.40111,0.00581395);
   return cutg;
}
