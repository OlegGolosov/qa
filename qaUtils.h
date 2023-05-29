using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;
using fourVector=LorentzVector<PtEtaPhiE4D<double>>;
using filteredDF=ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>;
using definedDF=ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void>;

TChain* makeChain(string& filename, const char* treename) {
  cout << "Adding files to chain:" << endl;
  TChain *chain = new TChain(treename);
  if (filename.rfind(".root") < filename.size())
    chain->Add(filename.data());
  else {
    TFileCollection fc("fc", "", filename.c_str());
    chain->AddFileInfoList((TCollection*)fc.GetList());
  }
  chain->ls();
  return chain;
}

void saveHists(filteredDF &dd, vector <pair <vector<string>, TH1DModel>> &models, TDirectory &dir)
{
  map <string, RResultPtr<::TH1D >> hists;
  for (auto& model:models)
  {
    auto m=model.second;
    const char *v=model.first.at(0).c_str(), *w=model.first.at(1).c_str();
    if (!w[0])
    {
      if(m.fName=="") 
        m.fName=Form("h_%s", v);
      if(m.fTitle=="") 
        m.fTitle=Form("%s;%s;nEntries", v, v);
      if(m.fNbinsX==0)
        hists.emplace(m.fName, dd.Histo1D(v));
      else
        hists.emplace(m.fName, dd.Histo1D(m, v));
    }
    else
    {
      if(m.fName=="") 
        m.fName=Form("h_%s_%s", v, w);
      if(m.fTitle=="") 
        m.fTitle=Form("%s {%s};%s;nEntries", v, w, v);
      if(m.fNbinsX==0)
        hists.emplace(m.fName, dd.Histo1D(v, w));
      else
        hists.emplace(m.fName, dd.Histo1D(m, v, w)); 
    }
  }

  dir.cd();
  for (auto& hist:hists)
  {
    cout << hist.first << endl;
    hist.second->Write(hist.first.c_str());
  }
}

void saveHists(filteredDF &dd, vector <pair <vector<string>, TH2DModel>> &models, TDirectory &dir)
{
  map <string, RResultPtr<::TH2D >> hists;
  for (auto& model:models)
  {
    auto m=model.second;
    const char *v1=model.first.at(0).c_str(), *v2=model.first.at(1).c_str(), *w=model.first.at(2).c_str();
    if (!w[0])
    {
      if (m.fName=="") 
        m.fName=Form("h2_%s_%s", v1, v2);
      if (m.fTitle=="") 
        m.fTitle=Form("%s:%s;%s;%s;nEntries", v2, v1, v1, v2);
      hists.emplace(m.fName, dd.Histo2D(m, v1, v2));
    }
    else
    {
      if (m.fName=="") 
        m.fName=Form("h2_%s_%s_%s", v1, v2, w);
      if (m.fTitle=="") 
        m.fTitle=Form("%s:%s {%s};%s;%s;nEntries", v2, v1, w, v1, v2);
      hists.emplace(m.fName, dd.Histo2D(m, v1, v2, w));
    }
  }

  dir.cd();
  for (auto& hist:hists)
  {
    cout << hist.first << endl;
    hist.second->Write();
  }
}

void saveHists(filteredDF &dd, vector <pair <vector<string>, TH3DModel>> &models, TDirectory &dir)
{
  map <string, RResultPtr<::TH3D >> hists;
  for (auto& model:models)
  {
    auto m=model.second;
    const char *v1=model.first.at(0).c_str(), *v2=model.first.at(1).c_str(), *v3=model.first.at(2).c_str(), *w=model.first.at(3).c_str();
    if (!w[0])
    {
      if (m.fName=="") 
        m.fName=Form("h3_%s_%s_%s", v1, v2, v3);
      if (m.fTitle=="") 
        m.fTitle=Form("%s:%s:%s;%s;%s;%s;nEntries", v3,v2, v1, v1, v2, v3);
      hists.emplace(m.fName, dd.Histo3D(m, v1, v2, v3));
    }
    else
    {
      if (m.fName=="") 
        m.fName=Form("h3_%s_%s_%s_%s", v1, v2, v3, w);
      if (m.fTitle=="") 
        m.fTitle=Form("%s:%s:%s {%s};%s;%s;%s;nEntries", v3, v2, v1, w, v1, v2, v3);
      hists.emplace(m.fName, dd.Histo3D(m, v1, v2, v3, w));
    }
  }

  dir.cd();
  for (auto& hist:hists)
  {
    cout << hist.first << endl;
    hist.second->Write();
  }
}

void saveHists(filteredDF &dd, vector <pair <vector<string>, TProfile1DModel>> &models, TDirectory &dir)
{
  map <string, RResultPtr<::TProfile >> hists;
  for (auto& model:models)
  {
    auto m=model.second;
    const char *v1=model.first.at(0).c_str(), *v2=model.first.at(1).c_str(), *w=model.first.at(2).c_str();
    if (!w[0])
    {
      if (m.fName=="") 
        m.fName=Form("p_%s_%s", v1, v2);
      if (m.fTitle=="") 
        m.fTitle=Form("#LT%s#GT:%s;%s;#LT%s#GT", v2, v1, v1, v2);
      hists.emplace(m.fName, dd.Profile1D(m, v1, v2));
    }
    else
    {
      if (m.fName=="") 
        m.fName=Form("p_%s_%s_%s", v1, v2, w);
      if (m.fTitle=="") 
        m.fTitle=Form("#LT%s#GT:%s {%s};%s;#LT%s#GT", v2, v1, w, v1, v2);
      hists.emplace(m.fName, dd.Profile1D(m, v1, v2, w));
    }
  }

  dir.cd();
  for (auto& hist:hists)
  {
    cout << hist.first << endl;
    hist.second->Write();
  }
}

void saveHists(filteredDF &dd, vector <pair <vector<string>, TProfile2DModel>> &models, TDirectory &dir)
{
  map <string, RResultPtr<::TProfile2D >> hists;
  for (auto& model:models)
  {
    auto m=model.second;
    const char *v1=model.first.at(0).c_str(), *v2=model.first.at(1).c_str(), *v3=model.first.at(2).c_str(), *w=model.first.at(3).c_str();
    if (!w[0])
    {
      if (m.fName=="") 
        m.fName=Form("p2_%s_%s_%s", v1, v2, v3);
      if (m.fTitle=="") 
        m.fTitle=Form("#LT%s#GT:%s:%s;%s;%s;#LT%s#GT", v3, v2, v1, v1, v2, v3);
      hists.emplace(m.fName, dd.Profile2D(m, v1, v2, v3));
    }
    else
    {
      if (m.fName=="") 
        m.fName=Form("p2_%s_%s_%s_%s", v1, v2, v3, w);
      if (m.fTitle=="") 
        m.fTitle=Form("#LT%s#GT:%s:%s {%s};%s;%s;#LT%s#GT", v3, v2, v1, w, v1, v2, v3);
      hists.emplace(m.fName, dd.Profile2D(m, v1, v2, v3, w));
    }
  }

  dir.cd();
  for (auto& hist:hists)
  {
    cout << hist.first << endl;
    hist.second->Write();
  }
}

void unfoldTH2 (RResultPtr<::TH2D> _h2, string projectionAxis="x", int rebinX=1, int rebinY=1, int firstBin=1, int lastBin=-1, bool log=true, const char *format="\%.0f")
{
  auto h2=(TH2*)_h2->Clone();
  h2->RebinX(rebinX);
  h2->RebinY(rebinY);
  TAxis *axis;
  const char* axisTitle;
  if (projectionAxis=="x")
  {
    axis=h2->GetYaxis();
    axisTitle=h2->GetXaxis()->GetTitle();
  }
  else if (projectionAxis=="y")
  {
    axis=h2->GetXaxis();
    axisTitle=h2->GetYaxis()->GetTitle();
  }
  else 
  {
    cout << "Not valid projection axis: " << projectionAxis << endl;
    return;
  }
  if (lastBin<0) 
    lastBin=axis->GetNbins();
  auto hs=new THStack(Form("hs_%s_log%i", h2->GetName(), log), h2->GetTitle());
  auto c=new TCanvas(Form("c_%s_log%i", h2->GetName(), log), h2->GetTitle(), 800, 600);
  int npads=lastBin-firstBin+1;
  int npadsx = int (ceil (sqrt (npads)));
  int npadsy = int (ceil (1. * npads / npadsx));
  c->Divide(npadsx, npadsy);
  const char* namePattern=Form("%s-%s", format, format);
  for (int i=firstBin;i<=lastBin;i++)
  {
    TH1* h;
    float binMin=axis->GetBinLowEdge(i), binMax=axis->GetBinUpEdge(i);
    if (projectionAxis=="x")
      h=h2->ProjectionX(Form(namePattern, binMin, binMax), i, i);
    else if (projectionAxis=="y")
      h=h2->ProjectionY(Form(namePattern, binMin, binMax), i, i);
    h->SetTitle(Form("%s;%s;N_{entries}", h->GetName(), axisTitle));
    h->SetLineColor(i);
    c->cd(i);
    if(log) gPad->SetLogy();
    h->Draw();
    hs->Add(h);
  }
  c->Write();
  hs->Write();
  delete h2;
}
