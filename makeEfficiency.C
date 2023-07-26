void makeEfficiency(const char *nameIn="qa_mc.root", const char *nameOut="pbpb13.eff.root")
{
  TFile fIn(nameIn), fOut(nameOut, "recreate");

  vector<pair<int,string>> particles=
  {
    {-211, "pionneg"},
    {211, "pionpos"},
    {2212, "proton"},
  };

  for(auto &part:particles)
  {
    auto eff=(TH2*)fIn.Get(Form("h2_trY_trPt_%s", part.second.c_str()))->Clone(Form("h2_eff_trY_trPt_%i", part.first));
    auto sim=(TH2*)fIn.Get(Form("h2_simY_simPt_sim_%s_prim", part.second.c_str()));
    eff->Rebin2D(5,5);
    sim->Rebin2D(5,5);
    eff->Divide(sim);
    fOut.cd();
    eff->Write();
  }

  fOut.Close(); 
}
