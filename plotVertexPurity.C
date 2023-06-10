void plotVertexPurity(string targetInFileName="~/desktop/analysis/na61_qa/pbpb/13agev/16_043/pbpb13_16_043_goodEventBeforeVtxZcut_t4.qa.root")
{
  const char* histName="h2_vtxZ_vtxNtrFit";
  int maxNtracks=100;
  float vtxZnormMin=-678, vtxZnormMax=-660;
  string targetOutFileName=regex_replace(targetInFileName, regex(".qa.root"), "_targetOut.qa.root");
  string vtxZcutFileName=regex_replace(targetInFileName, regex(".qa.root"), ".vtxPurity.root");
  cout << targetInFileName << endl;
  cout << targetOutFileName << endl;
  cout << vtxZcutFileName << endl;
  TFile targetInFile(targetInFileName.c_str());
  TFile targetOutFile(targetOutFileName.c_str());
  TFile vtxZcutFile(vtxZcutFileName.c_str(), "recreate");
  auto h2In=(TH2F*)targetInFile.Get(histName);
  auto h2Out=(TH2F*)targetOutFile.Get(histName);
//  h2In->RebinX(2);
//  h2Out->RebinX(2);
  auto h2vtxPurity=(TH2F*)h2In->Clone("h2vtxPurity");
  h2vtxPurity->Reset("ICES");
  auto unity=new TF1("unity","1", h2In->GetXaxis()->GetXmin(), h2In->GetXaxis()->GetXmax());
  for (int i=2;i<=h2In->GetNbinsY();i++)
  {
    THStack hs(Form("hs_%s_%i", histName, i), Form("hs_%s_%i;Z_{vtx} (cm);N_{events}", histName, i));
    auto hIn=h2In->ProjectionX(Form("%s_in_%i", histName, i), i, i);
    auto hOut=h2Out->ProjectionX(Form("%s_out_%i", histName, i), i, i);
    hIn->SetLineColor(kBlack);
    hOut->SetLineColor(kRed);
    hIn->SetTitle("target-in");
    hOut->SetTitle("target-out");
    float integralIn=hIn->Integral(hIn->GetXaxis()->FindBin(vtxZnormMin), hIn->GetXaxis()->FindBin(vtxZnormMax), "w");
    float integralOut=hOut->Integral(hOut->GetXaxis()->FindBin(vtxZnormMin), hOut->GetXaxis()->FindBin(vtxZnormMax), "w");
    hOut->Scale(integralIn/integralOut);
    hs.Add(hIn);
    hs.Add(hOut,"hist");
    auto hVtxPurity=(TH1*)hOut->Clone(Form("hVtxPurity_%i", i));
    hVtxPurity->Divide(hIn);
    for(int j=1;j<=hVtxPurity->GetNbinsX();j++)
    {
      auto val=hVtxPurity->GetBinContent(j);
      auto err=hVtxPurity->GetBinError(j);
      if(!isnan(val) && !isinf(val))
      {
        val=hVtxPurity->GetBinContent(j);
        err=hVtxPurity->GetBinError(j);
      }
      else
      {
        val=0.5*(hVtxPurity->GetBinContent(j-1)+hVtxPurity->GetBinContent(j+1));
        err=0.5*(hVtxPurity->GetBinError(j-1)+hVtxPurity->GetBinError(j+1));
        cout << val << "\t" << err << endl;
        if(isnan(val) || isinf(val))
        {
          val=0;
          err=0;
        }
      }
      h2vtxPurity->SetBinContent(j, i, val);
      h2vtxPurity->SetBinError(j, i, err);
    }
    hVtxPurity->Write();
    hs.Write();
  }
  h2In->Write(Form("%s_In", h2In->GetName()));
  h2Out->Write(Form("%s_Out", h2Out->GetName()));
  h2vtxPurity->Write();
  vtxZcutFile.Close(); 
}
