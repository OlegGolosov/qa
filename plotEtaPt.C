void plotEtaPt(TString fInName="qa13.root", float pBeam=13)
{
  int cbmCorner=2;
  float cbmXshift = 0.;
  float cbmYshift = 0.;
  vector<vector<float>> cbmPosition = 
  {
    {0.13, 0.85},
    {0.58, 0.85},
    {0.13, 0.15},
    {0.58, 0.15},
  }; 
  int kinCorner=1;
  float kinXshift = 0.07;
  float kinYshift = 0.;
  vector<vector<float>> kinPosition = 
  {
    {0.13, 0.85},
    {0.59, 0.85},
    {0.13, 0.15},
    {0.59, 0.15},
  }; 
  vector <TString> partNames={"pionneg","pionpos","proton"};
  vector <int> pids={-211,211,2212};
  vector <TString> partTitles={"#scale[1.2]{#pi^{-}}","#scale[1.2]{#pi^{+}}","proton"};
  vector <vector<float>> kinCuts;
  if(pBeam==13)
    kinCuts={{0.8,1.2,0.2,2.2},{0.8,1.2,0.2,2.2},{0.8,1.2,0.4,2.2}};
  if(pBeam==30)
    kinCuts={{0.6,1.0,0.2,2.2},{0.2,1.2,0.2,2.2},{0.2,0.8,0.8,3.0}};
  vector <vector<float>> kinPlotRanges={{1.01,5.49,0.01,2.49},{1.01,5.49,0.01,2.49},{1.01,5.49,0.0,2.99}};
  vector<float> psdRadii={6.7,20,40,60};
  vector<int> psdColors={kBlue,kRed,kGreen+2};
  vector<TString> psdTitles={"PSD1","PSD2","PSD3"};

  float eBeam=sqrt(pBeam*pBeam+0.938*0.938);
  float yMid = 0.25 * log ((eBeam + pBeam) / (eBeam - pBeam));

  vector <TBox*> psdBoxes;
  for (uint i=0;i<psdRadii.size()-1;i++)
  {
    float psdRmin=psdRadii.at(i), psdRmax=psdRadii.at(i+1), psdz=800.;
    float psdPtMin=0.01, psdPtMax=1.8, etaShift=0.02;
    float psdEtaMin=-log(0.5*psdRmax/psdz)+etaShift;
    float psdEtaMax=-log(0.5*psdRmin/psdz);
    psdBoxes.push_back(new TBox(psdEtaMin, psdPtMin, psdEtaMax, psdPtMax));
    psdBoxes.at(i)->SetLineColor(psdColors.at(i));
    psdBoxes.at(i)->SetLineWidth(4);
    psdBoxes.at(i)->SetFillStyle(0);
  }
  bool logz=1;
  gStyle->SetOptStat(0);
  //c->Divide(partNames.size(),1/*,0,0*/);
  auto f=new TFile(fInName,"read");
  for (int i=0;i<partNames.size();i++)
  { 
    auto c=new TCanvas(Form("c_%s",partNames.at(i).Data()),"c",800,600);
    float m=TDatabasePDG::Instance()->GetParticle(pids.at(i))->Mass();
    TH2* h = (TH2*)f->Get(Form("h2_trEta_trPt_%s",partNames.at(i).Data()));
    c->cd(i+1); 
    c->SetLogz(logz);
    c->SetRightMargin(0.12);
    
    h->GetXaxis()->SetTitle("#eta");
    h->GetXaxis()->SetNdivisions(505);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetRangeUser(kinPlotRanges.at(i).at(0),kinPlotRanges.at(i).at(1));
    h->GetYaxis()->SetTitle("p_{T} (GeV/#it{c})");
    h->GetYaxis()->SetNdivisions(505);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetRangeUser(kinPlotRanges.at(i).at(2),kinPlotRanges.at(i).at(3));
    h->GetZaxis()->SetNdivisions(505);
    h->GetZaxis()->SetLabelSize(0.05);
    
    auto *track_box=new TBox(kinCuts.at(i).at(0), kinCuts.at(i).at(2), kinCuts.at(i).at(1), kinCuts.at(i).at(3));
 
    TAxis *ax=h->GetXaxis();
    TAxis *ay=h->GetYaxis();
    for (uint i=1;i<=ax->GetNbins();i++)
    {
      double eta=ax->GetBinCenter(i);
      for (uint j=1;j<=ay->GetNbins();j++)
      {
        double pt=ay->GetBinCenter(j);
        TLorentzVector v;
      	v.SetPtEtaPhiM(pt,eta,0.,m);
      	float y=v.Rapidity()-yMid;
        if (!track_box->IsInside(y,pt))
          h->SetBinContent(i,j,0);
      } 
    }
    h->Draw("colz");
    for (uint b=0;b<psdBoxes.size();b++)
    {
      auto box=psdBoxes.at(b); 
      box->Draw();
      auto tex = new TLatex(box->GetX2()-0.05,box->GetY2()-0.47," "+psdTitles.at(b));
      tex->SetTextSize(0.05);
      tex->SetTextFont(42);
      tex->SetTextAngle(90);
      tex->SetTextColor(psdColors.at(b));
      tex->Draw();
    }
    float kinx=kinPosition.at(kinCorner-1).at(0)+kinXshift;
    float kiny=kinPosition.at(kinCorner-1).at(1)+kinYshift;
    TLatex *tex = new TLatex(kinx,kiny,Form("%s, #it{y}#in [%.1f, %.1f]", partTitles.at(i).Data(), kinCuts.at(i).at(0), kinCuts.at(i).at(1)));
    tex->SetNDC();
    tex->SetTextSize(0.045);
    tex->SetTextFont(42);
    tex->Draw();
    tex = new TLatex(kinx,kiny-0.05,Form("p_{T}#in [%.1f, %.1f] GeV/#it{c}", kinCuts.at(i).at(2), kinCuts.at(i).at(3)));
    tex->SetNDC();
    tex->SetTextSize(0.045);
    tex->SetTextFont(42);
    tex->Draw();
    float cbmx=cbmPosition.at(cbmCorner-1).at(0)+cbmXshift;
    float cbmy=cbmPosition.at(cbmCorner-1).at(1)+cbmYshift;
    tex = new TLatex(cbmx,cbmy,"NA61/SHINE performance");
    tex->SetNDC();
    tex->SetTextSize(0.045);
    tex->SetTextFont(42);
    tex->Draw();
    tex = new TLatex(cbmx,cbmy-0.05,Form("Pb+Pb @ %i#it{A} GeV/#it{c}", int(pBeam)));
    tex->SetNDC();
    tex->SetTextSize(0.045);
    tex->SetTextFont(42);
    tex->Draw();
//    tex = new TLatex(cbmx,cbmy-0.1,"DCM-QGSM-SMM");
//    tex->SetNDC();
//    tex->SetTextSize(0.045);
//    tex->SetTextFont(42);
//    tex->Draw();
    TString saveName=partNames.at(i)+"_psdEtaPt";
//    c->Print(saveName+".pdf");
  }
}
