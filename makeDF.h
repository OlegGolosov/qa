#include "cut_dEdx_m2_eneg.C"       
#include "cut_dEdx_m2_epos.C"     
#include "cut_dEdx_m2_kaonneg.C"  
#include "cut_dEdx_m2_kaonpos.C"
#include "cut_dEdx_m2_pionneg.C"
#include "cut_dEdx_m2_pionpos.C"
#include "cut_dEdx_m2_proton.C"

float nSigma(float x, const TF1 &f)
{
  auto p=f.GetParameters();
  return fabs(x-p[1])/p[2];
};

float nSigma(float x, float y, const TF2 &f2)
{
  auto p=f2.GetParameters();
  return sqrt(pow((x-p[1])/p[2],2)+pow((y-p[3])/p[4],2));
};

float r(float x, float y)
{
  return sqrt(x*x+y*y);
}

RVec<int> matchedPdg(RVec<int> trSimIndex, RVec<int> simPdg)
{
  RVec<int> pdgs;
  for (auto& ind:trSimIndex)
  {
    if(ind>-1)
      pdgs.push_back(simPdg.at(ind));
    else
      pdgs.push_back(0);
  }
  return pdgs;
}

RVec<TLorentzVector> trSimMom4(RVec<int> trSimIndex, RVec<TLorentzVector> simMom4)
{
  RVec<TLorentzVector> mom4;
  for (auto& ind:trSimIndex)
  {
    if(ind>-1)
      mom4.push_back(simMom4.at(ind));
    else
      mom4.push_back({0,0,0,0});
  }
  return mom4;
}

RVec<int> trSimMotherId(RVec<int> trSimIndex, RVec<int> simMotherId)
{
  RVec<int> motherId;
  for (auto& ind:trSimIndex)
  {
    if(ind>-1)
      motherId.push_back(simMotherId.at(ind));
    else
      motherId.push_back(-999);
  }
  return motherId;
}

auto getMass(RVec<int> pdgs)
{
  RVec<float> masses; 
  for(auto& pdg:pdgs)
  {
    float mass;
    auto part=TDatabasePDG::Instance()->GetParticle(pdg);
    if(!part) mass=0;
    else if(abs(pdg)>10000) mass=(pdg%10000/10)*0.931;
    else mass=part->Mass();
    masses.push_back(mass); 
  }
  return masses;
}

bool isSimulation=false;

filteredDF makeDataFrame (definedDF &d, TDirectory *dir=nullptr)
{
  int pBeam=round(*d.Range(0,100).Mean("beamPz"));
  isSimulation=*d.Range(0,100).Mean("isSimulation");
  auto mp=0.938, mNucl=0.931;
  auto eBeam=sqrt(pBeam*pBeam+mp*mp);
  auto yBeam=0.25*log((eBeam+pBeam)/(eBeam-pBeam));
  auto sqrt2MoverT=sqrt(2*mNucl/eBeam);
  printf("pBeam = %i \t yBeam = %f \n", pBeam, yBeam);

  auto rapidityCMS = [yBeam] (RVec<TLorentzVector> mom4){RVec<float> y; for(auto &mom:mom4) y.push_back(mom.Rapidity()-yBeam); return y;};

  auto fPid=new TFile(Form("pbpb%i.pid.root", pBeam));
  auto pidGetterNeg=(Pid::Getter*)fPid->Get("pidGetterNeg");
  auto pidGetterPos=(Pid::Getter*)fPid->Get("pidGetterPos");
  float pidPurity=0.9;
  
  auto getPid_dEdx=[pidGetterNeg, pidGetterPos, pidPurity](RVec<float> log20p, RVec<float> dEdx, RVec<int> charge)
  {
    RVec<int> pids(log20p.size(), 0);
    Pid::Getter *getter;
    for(int i=0;i<log20p.size();i++)
    {
      if(charge.at(i)>0) getter=pidGetterPos;
      else getter=pidGetterNeg;
      if(!getter) continue;
      pids.at(i)=getter->GetPid(log20p.at(i), dEdx.at(i), pidPurity);
    }
    return pids; 
  };

  double log20pmin=0, log20pmax=10;
  int nBinsLog20p=200;
  map <int, TH1F*> hPidEffCorrs;
  for (auto pid:{-211,211,2212})
  {
    Pid::Getter *getter;
    auto charge=abs(pid)/pid;
    if(charge>0) getter=pidGetterPos;
    else getter=pidGetterNeg;
    if(!getter) continue;
    hPidEffCorrs.emplace(pid, new TH1F(Form("pidEffCorr_%i", pid), Form("pidEffCorr_%i", pid), nBinsLog20p, log20pmin, log20pmax));
    hPidEffCorrs.at(pid)->SetDirectory(0);

    for(int i=1;i<=nBinsLog20p;i++)
    {
      auto h=hPidEffCorrs.at(pid);
      auto eff=1./getter->GetEfficiency(h->GetXaxis()->GetBinCenter(i), pid, pidPurity);
      if (isnan(eff) || isinf(eff) || eff>10) eff=1.;
      h->SetBinContent(i, eff);
    }

    if(dir)
    {
      dir->cd();
      hPidEffCorrs.at(pid)->Clone()->Write();
    }
  }

  auto pidEfficiencyCorrection=[hPidEffCorrs](RVec<float> log20p, RVec<int> pid)
  {
    RVec <float> eff(pid.size(), 1.);
    for (int i=0;i<pid.size();i++)
    {
      if(hPidEffCorrs.find(pid.at(i))==hPidEffCorrs.end())
        continue;
      auto h=hPidEffCorrs.at(pid.at(i));
      eff.at(i)=h->GetBinContent(h->GetXaxis()->FindBin(log20p.at(i)));
    }
    return eff; 
  };

  map <int, TCutG*> dEdx_m2_cuts =
  {
    {-11, get_cut_dEdx_m2_eneg()},
    {11, get_cut_dEdx_m2_epos()},
    {-321, get_cut_dEdx_m2_kaonneg()},
    {321, get_cut_dEdx_m2_kaonpos()},
    {-211, get_cut_dEdx_m2_pionneg()},
    {211, get_cut_dEdx_m2_pionpos()},
    {2212, get_cut_dEdx_m2_proton()},
//    {1000010020, get_cut_dEdx_m2_deuteron()},
  };

  auto getPid_dEdx_m2=[dEdx_m2_cuts](RVec<float> dEdx, RVec<float> m2, RVec<int> charge)
  {
    RVec<int> pids(dEdx.size());
    for(int i=0;i<pids.size();i++)
      for(auto &p:dEdx_m2_cuts)
        if(p.second->IsInside(dEdx.at(i), m2.at(i)) && charge.at(i)==(p.first/abs(p.first)))
        {
          pids.at(i)=p.first;
          continue;
        }
    return pids; 
  };

  TDatabasePDG::Instance(); // for some reason needed at lxplus
/*                          
  vector <float> xproton={4.63796, 6.15895, 7.91560, 5.10925, 3.73822, 3.45973, 4.12382, 4.63796};
  vector <float> yproton={1.14662, 1.37218, 1.11654, 0.657895, 0.620301, 0.883459, 1.07143, 1.14662};
  TCutG proton("proton", xproton.size(), &xproton[0], &yproton[0]);
  vector <float> xpionneg={-5.94473, -5.55913, -4.68081, -3.63111, -2.64567, -1.89589, -1.59597, -1.21037, -1.78877, -2.75278, -3.7168 , -5.04499, -5.92331, -6.11611, -6.11611, -5.94473};
  vector <float> ypionneg={1.54511, 1.57519, 1.5, 1.40226, 1.30451, 1.28947, 1.38722, 1.17669, 0.800752, 0.733083, 0.800752, 1.01128, 1.32707, 1.5, 1.5, 1.54511};
  TCutG pionneg("pionneg", xpionneg.size(), &xpionneg[0], &ypionneg[0]);

  auto getPid_dEdx_cut=[proton, pionneg](RVec<float> Qlog20p, RVec<float> dEdx)
  {
    RVec<int> pids;
    for(int i=0;i<Qlog20p.size();i++)
    {
      int pid=0;
      if(proton.IsInside(Qlog20p.at(i), dEdx.at(i))) pid=2212;
      else if(pionneg.IsInside(Qlog20p.at(i), dEdx.at(i))) pid=-211;
      pids.push_back(pid);
    }
    return pids; 
  };
*/
  TFile fCentr(Form("pbpb%i.centrality_Epsd.root", pBeam));
  auto centrGetter=(Centrality::Getter*)fCentr.Get("centr_getter_1D");
  if(!centrGetter) centrGetter=(Centrality::Getter*)fCentr.Get("centr_getter_1d");

  auto getCentrality=[centrGetter](float E)
  {
    if(!centrGetter) return 0.f;
    return centrGetter->GetCentrality(E+21.); 
  };

  TFile fVtxPurity(Form("pbpb%i.vtxPurity.root", pBeam));
  auto h2VtxPurity=(TH2*)fVtxPurity.Get("h2vtxPurity");
  if(h2VtxPurity) h2VtxPurity->SetDirectory(0);

  auto getVtxPurity=[h2VtxPurity](float vtxZ, int vtxNtrFit)
  {
    if(!h2VtxPurity) return 1.f;
    int binX=h2VtxPurity->GetXaxis()->FindBin(vtxZ);
    int binY=h2VtxPurity->GetYaxis()->FindBin(vtxNtrFit);
    float purity=1.-h2VtxPurity->GetBinContent(binX, binY);
//    if(purity<0.01 || purity>1 || isnan(purity) || isinf(purity)) cout << vtxZ << "\t" << vtxNtrFit << "\t" << purity << endl;
    return purity; 
  };

  TFile fEff(Form("pbpb%i.eff.root", pBeam));
  map <int, TH2*> h2Effs;
  cout << "Tracking efficiency histograms:\n";
  for(auto &pdg:{-211, 211, 2212})
  {
    cout << endl << Form("h2_eff_trY_trPt_%i", pdg);
    auto h=(TH2*)fEff.Get(Form("h2_eff_trY_trPt_%i", pdg));
    if(!h) continue;
    h->SetDirectory(0);
    h2Effs.emplace(pdg, h);
    cout << "\tFOUND!!!" << endl;
  }

  auto trackingEfficiencyCorrection=[h2Effs](RVec<float> y, RVec<float> pt, RVec<int> pdg)
  {
    RVec <float> eff(pdg.size(), 1.);
    for (int i=0;i<pdg.size();i++)
    {
      if(h2Effs.find(pdg.at(i))==h2Effs.end())
        continue; 
      auto h2Eff=h2Effs.at(pdg.at(i));
      int binX=h2Eff->GetXaxis()->FindBin(y.at(i));
      int binY=h2Eff->GetYaxis()->FindBin(pt.at(i));
      auto val=1./h2Eff->GetBinContent(binX, binY);
      if(!isinf(val)) eff.at(i)=val;
    }
    return eff; 
  };

  string pileUpCut, goodWfa="Sum(abs(timesS1_1)<4000)==1";
  string goodS1S2="fabs(adcS1-150)<35 && fabs(adcS2-195)<55 && adcS1*1.258+adcS2<390";
  float vtxZnominal=-591.9, vtxZdelta=2;
  string psdModSub1="psdModId<=16 || psdModId==45";
  string psdModSub2="17<=psdModId && psdModId<=28";
  string psdModSub3="29<=psdModId && psdModId<=44";
  string psdModSub0="psdModId==6 || psdModId==7 || psdModId==10 || psdModId==11 || psdModId==45";
  string psdModSub0wo45="psdModId==6 || psdModId==7 || psdModId==10 || psdModId==11";
  string psdModSub1wo45="psdModId<=16";
  if(pBeam==13)
  {
    pileUpCut="psdE<3800*(1-Mgood/230.)";
  }
  else if(pBeam==30)
  {
    pileUpCut="psdE<8000*(1-pow(Mgood/400., 2))";
  }
  else if(pBeam==41)
  {
    pileUpCut="true";
    goodWfa="true";
    goodS1S2="true";
    vtxZnominal=-581.1;
    psdModSub1="psdModId<=4";
    psdModSub2="5<=psdModId && psdModId<5+1*24";      //psd-like
    psdModSub3="5+1*24<=psdModId && psdModId<5+6*24"; //psd-like
//    psdModSub2="5<=psdModId && psdModId<5+3*24";      //psd-like 3x3 rings
//    psdModSub3="5+3*24<=psdModId && psdModId<5+6*24"; //psd-like 3x3 rings
//    psdModSub2="5<=psdModId && psdModId<5+5*24";       //half-rings 
//    psdModSub3="5+5*24<=psdModId && psdModId<5+10*24"; //half-rings
    psdModSub0="5+6*24<=psdModId && psdModId<5+10*24"; 
    psdModSub0wo45="psdModSub0";
    psdModSub1wo45="psdModSub1";
  }

  
//  cout << "Number of events: " << *(d.Count()) << endl;
  auto dd=d
    //.Range(0,1000)
    .Filter("eventId>0")
    .Define("noDeltaEvents", "Sum(trCharge>0)>0")
    .Define("vtxChi2ndf", "vtxChi2/vtxNdf")
    .Define("psdModPhi", "atan2(psdModY,psdModX)")
    .Define("psdModSub1", psdModSub1)
    .Define("psdModSub2", psdModSub2)
    .Define("psdModSub3", psdModSub3)
    .Define("psdModSub0", psdModSub0)
    .Define("psdModSub0wo45", psdModSub0wo45)
    .Define("psdModSub1wo45", psdModSub1wo45)
    .Define("psdE", "Sum(psdModE)")
    .Define("psd0E", "Sum(psdModE*psdModSub0)")
    .Define("psd1E", "Sum(psdModE*psdModSub1)")
    .Define("psd2E", "Sum(psdModE*psdModSub2)")
    .Define("psd12E", "psd1E+psd2E")
    .Define("psd3E", "Sum(psdModE*psdModSub3)")
    .Define("trId", "RVec<int> id; for(int i=0;i<trPt.size();i++)id.push_back(i); return id;")
    .Define("trMom", "RVec<TVector3> mom; for(int i=0;i<nTracks;i++){TVector3 v;v.SetPtEtaPhi(trPt[i], trEta[i], trPhi[i]);mom.push_back(v);} return mom;")
    .Define("trPx", "RVec<float> px; for(auto &mom:trMom)px.push_back(mom.Px()); return px;")
    .Define("trPy", "RVec<float> py; for(auto &mom:trMom)py.push_back(mom.Py()); return py;")
    .Define("trPz", "RVec<float> pz; for(auto &mom:trMom)pz.push_back(mom.Pz()); return pz;")
    .Define("trP", "RVec<float> trP; for(auto &mom:trMom)trP.push_back(mom.Mag()); return trP;")
    .Define("trQp", "trCharge*trP")
    .Define("trLog10p", "log10(trP)")
    .Define("trLog20p", "log(20*trP)")
    .Define("trQoverP", "trCharge/trP")
    .Define("trQlog20p", "trCharge*log(20*trP)")
    .Define("trQdEdx", "trCharge*trdEdx")
    .Define("trNclustVTPC", "trNclustVTPC1+trNclustVTPC2")
    .Define("trNclustTPC", "trNclustVTPC1+trNclustVTPC2+trNclustMTPC")
    .Define("trNclustPotTPC", "trNclustPotVTPC1+trNclustPotVTPC2+trNclustPotMTPC")
    .Define("trNclustToPot", "RVec<float> ratio(nTracks, -10); for(int i=0;i<nTracks;i++)if(trNclustPotTPC.at(i)>0)ratio.at(i)=(float)trNclustTPC.at(i)/trNclustPotTPC.at(i); return ratio;")
    .Define("trChi2ndf", "trChi2/trNdf")
    .Define("trNclustCut", "trNclustVTPC>=15 && trNclustTPC>=30")
    .Define("trNclustToPotCut", "trNclustToPot>0.55 && trNclustToPot<1.1")
    .Define("trDcaCut", "abs(trDcaZ)<.5 && abs(trDcaX-0.02)<2 && abs(trDcaY-0.02)<1")
    .Define("trGoodNclust", "trNclustCut && trNclustToPotCut")
    .Define("trChi2ndfCut", "trChi2ndf<10")
    .Define("trChi2Cut", "trChi2<100000")
    .Define("trMid", "trPt>0.02 && trEta<3.5")
    .Define("trGood", "trGoodNclust && trDcaCut && trChi2Cut")
    .Define("vtxPurity", getVtxPurity, {"vtxZ", "vtxNtrFit"})
    .Define("trGoodPos", "trGood && trCharge>0")
    .Define("trGoodNeg", "trGood && trCharge<0")
    .Define("Mgood", "Sum(trGood)")
    .Define("MgoodMid", "Sum(trGood*trMid)")
  ;

  if(pBeam==13)
    dd=dd.Define("centrality", getCentrality, {"psdE"});
  else if(pBeam==30)
    dd=dd.Define("centrality", getCentrality, {"psd1E"});
  else if(pBeam==41)
  {
    dd=dd
      .Define("centrality", "centrality49")
      .Redefine("psdE", "vetoEcal")
    ;
  }

  if (!isSimulation)
  {
    dd=dd
      .Filter("trigT4||trigT2")
//      .Filter("trigT4")
      .Define("goodWfa", goodWfa)
      .Define("goodS1S2_16_011", "fabs(adcS1-145)<30 && fabs(adcS2-325)<175 && adcS1*1.258+adcS2<390")
      .Define("goodS1S2", goodS1S2)
      .Define("goodBpd", "fabs(bpd1x-0.1)<0.4 && fabs(bpd1y-0.2)<0.9 && fabs(bpd2x+0.2)<0.4 && fabs(bpd2y+0.1)<0.4 && fabs(bpd3x+0.5)<0.5 && fabs(bpd3y+0.25)<0.5")
      .Define("pileUpCut", pileUpCut)
      .Define("vtxXwrtBpd3", "vtxX-bpd3x")
      .Define("vtxYwrtBpd3", "vtxY-bpd3y")
      .Define("goodEventBeforeVtxCut", "goodS1S2 && goodWfa && pileUpCut && vtxFitPerfect && noDeltaEvents")
      .Define("goodEventBeforeVtxCut_16_011", "goodS1S2_16_011 && goodWfa && vtxFitPerfect")
      .Define("goodVtxPos_16_011", "fabs(vtxX+0.3045)<0.4355 && fabs(vtxY+0.321)<0.343 && fabs(vtxZ+592)<2")
      .Define("goodEvent_16_011", "goodEventBeforeVtxCut_16_011 && goodVtxPos_16_011")
//      .Define("trTofM2", "RVec<float> trM2(trId.size(), -999); for(int i=0;i<tofHitM2.size();i++) if(tofHitTrackId.at(i)>=0) trM2.at(tofHitTrackId.at(i))=tofHitM2.at(i); return trM2;")
//      .Define("trPid", getPid_dEdx_cut, {"trQp", "trdEdx"})
      .Define("trPid", getPid_dEdx, {"trLog20p", "trdEdx", "trCharge"})
      .Define("trPidEffCorr", pidEfficiencyCorrection, {"trLog20p", "trPid"})
      .Define("trPid_dEdx_m2", getPid_dEdx_m2, {"trdEdx", "trTofM2", "trCharge"})
      .Define("trGoodPosWeightNdEdx", "(trGood && trCharge>0)*trNdEdxClust")
      .Define("trGoodNegWeightNdEdx", "(trGood && trCharge<0)*trNdEdxClust")
      .Define("pionneg_dEdx_m2", "trGood && trPid_dEdx_m2==-211")
      .Define("kaonneg_dEdx_m2", "trGood && trPid_dEdx_m2==-321")
      .Define("eneg_dEdx_m2", "trGood && trPid_dEdx_m2==-11")
      .Define("pionpos_dEdx_m2", "trGood && trPid_dEdx_m2==211")
      .Define("kaonpos_dEdx_m2", "trGood && trPid_dEdx_m2==321")
      .Define("epos_dEdx_m2", "trGood && trPid_dEdx_m2==11")
      .Define("proton_dEdx_m2", "trGood && trPid_dEdx_m2==2212")
    ;
/*
    vector<int> pids{11,211,321,2212,1000010020};
    map <int, TH3*> purityHistMap;
    TFile fPurities("pbpb13.pid_kit.root");
    for(auto &charge:{-1,1})
      for(auto &pid:pids)
        for (auto &wst:{0,1})
        {
          const char *name=Form("%i_%i", pid*charge, wst);
          int key=charge*pid+wst;
          auto h=(TH3*)fPurities.Get(name);
          h->SetDirectory(0);
          purityHistMap.emplace(key, h); 
        }
    
    auto getPidPurity_kit=[purityHistMap, pids](RVec<float> p, RVec<float> pt, RVec<float> dEdx, RVec<int> charge,  RVec<float> px)
    {
      RVec<map<int,float>> purities;
      auto h=purityHistMap.begin()->second;
      for(int i=0;i<p.size();i++)
      {
        map<int,float> pur;
        auto wst=(charge.at(i)*px.at(i)<0);
        auto binx = h->GetXaxis()->FindBin(p.at(i)); 
        auto biny = h->GetYaxis()->FindBin(pt.at(i)); 
        auto binz = h->GetZaxis()->FindBin(log(dEdx.at(i)));
        for(auto &pid:pids)
        {
          int key=charge.at(i)*pid+wst;
          pur.emplace(pid, purityHistMap.at(key)->GetBinContent(binx, biny, binz));
        }
        cout << endl;
        purities.push_back(pur);
      }
      return purities; 
    };

    dd=dd
      .Define("trPidPurity", getPidPurity_kit, {"trP", "trPt", "trdEdx", "trCharge", "trPx"})
      .Define("purityElectron", "RVec<float> pur; for(auto& p:trPidPurity) pur.push_back(p.at(11)); return pur;")
      .Define("purityPion", "RVec<float> pur; for(auto& p:trPidPurity) pur.push_back(p.at(211)); return pur;")
      .Define("purityKaon", "RVec<float> pur; for(auto& p:trPidPurity) pur.push_back(p.at(321)); return pur;")
      .Define("purityProton", "RVec<float> pur; for(auto& p:trPidPurity) pur.push_back(p.at(2212)); return pur;")
      .Define("purityDeuteron", "RVec<float> pur; for(auto& p:trPidPurity) pur.push_back(p.at(1000010020)); return pur;")
      .Define("pionneg_kit", "trGood && trCharge<0 && purityPion>0.9")
      .Define("kaonneg_kit", "trGood && trCharge<0 && purityKaon>0.9")
      .Define("eneg_kit", "trGood && trCharge<0 && purityElectron>0.9")
      .Define("pionpos_kit", "trGood && trCharge>0 && purityPion>0.9")
      .Define("kaonpos_kit", "trGood && trCharge>0 && purityKaon>0.9")
      .Define("epos_kit", "trGood && trCharge>0 && purityElectron>0.9")
      .Define("proton_kit", "trGood && trCharge>0 && purityProton>0.9")
      .Define("deuteron_kit", "trGood && trCharge>0 && purityDeuteron>0.9")
    ;
*/
  }

  if (isSimulation)
  {
    dd=dd
      .Define("trPid", matchedPdg, {"trSimIndex", "simPdg"})
      .Define("simM", getMass, {"simPdg"})
      .Define("simMom4", "RVec<TLorentzVector> mom4; for(int i=0;i<nSim;i++){TLorentzVector v;v.SetPtEtaPhiM(simPt.at(i), simEta.at(i), simPhi.at(i), simM.at(i));mom4.push_back(v);} return mom4;")
      .Define("simY", rapidityCMS, {"simMom4"})
      .Define("trSimMom4", trSimMom4, {"trSimIndex", "simMom4"})
      .Define("trSimMotherId", trSimMotherId, {"trSimIndex", "simMotherId"})
      .Define("sim_pionneg_prim", "simMotherId==-1 && simPdg==-211")
      .Define("sim_pionpos_prim", "simMotherId==-1 && simPdg==211")
      .Define("sim_proton_prim", "simMotherId==-1 && simPdg==2212")
      .Define("trPrimary", "trSimMotherId==-1")
      .Redefine("trGood", "trGood && trPrimary")
      .Define("goodMcEventBeforeVtxCut", "vtxFitPerfect")
    ;
  }

  dd=dd
    .Define("trM", getMass, {"trPid"})
    .Define("trMom4", "RVec<TLorentzVector> mom4; for(int i=0;i<nTracks;i++){TLorentzVector v;v.SetPtEtaPhiM(trPt[i], trEta[i], trPhi[i], trM[i]);mom4.push_back(v);} return mom4;")
    .Define("trY", rapidityCMS, {"trMom4"})
    .Define("trY0", Form("trY/%f", yBeam))
    .Define("trUt0", Form("trPt/trM*%f", sqrt2MoverT))
//    .Define("trYproton", [yBeam](RVec<TLorentzVector> mom4){RVec<float> y; for(auto &mom:mom4){mom.SetPtEtaPhiM(mom.Pt(), mom.Eta(), mom.Phi(), TDatabasePDG::Instance()->GetParticle(2212)->Mass()); y.push_back(mom.Rapidity()-yBeam);} return y;}, {"trMom4"})
//    .Define("trYpion", [yBeam](RVec<TLorentzVector> mom4){RVec<float> y; for(auto &mom:mom4){mom.SetPtEtaPhiM(mom.Pt(), mom.Eta(), mom.Phi(), TDatabasePDG::Instance()->GetParticle(211)->Mass()); y.push_back(mom.Rapidity()-yBeam);} return y;}, {"trMom4"})
    .Define("pionneg", "trGood && trPid==-211")
    .Define("kaonneg", "trGood && trPid==-321")
    .Define("eneg", "trGood && trPid==-11")
    .Define("pionpos", "trGood && trPid==211")
    .Define("kaonpos", "trGood && trPid==321")
    .Define("epos", "trGood && trPid==11")
    .Define("proton", "trGood && trPid==2212")
    .Define("deuteron", "trGood && trPid==1000010020")
    .Define("trTrackingEffCorr", trackingEfficiencyCorrection, {"trY", "trPt", "trPid"})
    .Define("pionneg_effTr", "pionneg*trTrackingEffCorr")
    .Define("pionpos_effTr", "pionpos*trTrackingEffCorr")
    .Define("proton_effTr", "proton*trTrackingEffCorr")
    .Define("Mproton", "Sum(proton)")
    .Define("Mpionneg", "Sum(pionneg)")
    .Define("MprotonMid", "Sum(proton*trMid)")
    .Define("MpionnegMid", "Sum(pionneg*trMid)")
  ;

  if (!isSimulation)
  {
    dd=dd
      .Define("trTrackingPidEffCorr", "trTrackingEffCorr*trPidEffCorr")
      .Define("pionneg_effPid", "pionneg*trPidEffCorr")
      .Define("pionpos_effPid", "pionpos*trPidEffCorr")
      .Define("proton_effPid", "proton*trPidEffCorr")
      .Define("pionneg_effTrPid", "pionneg*trTrackingPidEffCorr")
      .Define("pionpos_effTrPid", "pionpos*trTrackingPidEffCorr")
      .Define("proton_effTrPid", "proton*trTrackingPidEffCorr")
    ;
  }

  auto temp=dd;
  if (isSimulation) dd=dd.Filter("goodMcEventBeforeVtxCut");
  if (!isSimulation) temp=dd.Filter("goodEventBeforeVtxCut");

  float vtxZfitOffset=0.2;
  TF1 f_vtxXfit("f_vtxXfit","gaus", -2, 2);
  TF1 f_vtxYfit("f_vtxYfit","gaus", -2, 2);
  TF1 f_vtxZfit("f_vtxZfit","gaus",vtxZnominal-vtxZfitOffset,vtxZnominal+vtxZfitOffset);
  auto h_vtxXfit=temp.Histo1D({"h_vtxXfit", "vtxXfit", 500, -2, 2}, "vtxX");
  h_vtxXfit->Fit(&f_vtxXfit, "w", "",  -2, 2);
  auto h_vtxYfit=temp.Histo1D({"h_vtxYfit", "vtxYfit", 500, -2, 2}, "vtxY");
  h_vtxYfit->Fit(&f_vtxYfit, "w", "",  -2, 2);
  auto h_vtxZfit=temp.Histo1D({"h_vtxZfit", "vtxZfit", 500, vtxZnominal-5,vtxZnominal+5}, "vtxZ");
  h_vtxZfit->Fit(&f_vtxZfit, "w", "", vtxZnominal-vtxZfitOffset, vtxZnominal+vtxZfitOffset);

  if(dir)
  {
    dir->cd();
    h_vtxXfit->Write();
    h_vtxYfit->Write();
    h_vtxZfit->Write();
  }

  dd=dd
    .Define("nSigmaVtxX", [f_vtxXfit](float z){return nSigma(z, f_vtxXfit);}, {"vtxX"})
    .Define("nSigmaVtxY", [f_vtxYfit](float z){return nSigma(z, f_vtxYfit);}, {"vtxY"})
    .Define("nSigmaVtxXY", r, {"nSigmaVtxX", "nSigmaVtxY"})
    .Define("goodVtxXY", "nSigmaVtxXY<=3")
    .Define("goodVtxZ", [vtxZnominal, vtxZdelta](float vtxZ){return fabs(vtxZ-vtxZnominal)<vtxZdelta;}, {"vtxZ"})
    .Define("goodVtxPos", "goodVtxZ && goodVtxXY")
    .Define("goodVtx", "vtxFitPerfect && goodVtxPos")
  ;

  if (isSimulation)
  {
    dd=dd
      .Define("goodMcEvent", "goodMcEventBeforeVtxCut && goodVtx")
    ;
  }

  if (!isSimulation)
  {
    TF1 f_vtxXwrtBpd3fit("f_vtxXwrtBpd3fit","gaus", -.2, .2);
    TF1 f_vtxYwrtBpd3fit("f_vtxYwrtBpd3fit","gaus", -.2, .2);
    auto h_vtxXwrtBpd3fit=temp.Histo1D({"h_vtxXwrtBpd3fit", "vtxXwrtBpd3fit", 500, -.2, .2}, "vtxXwrtBpd3");
    h_vtxXwrtBpd3fit->Fit(&f_vtxXwrtBpd3fit, "w", "",  -.2, .2);
    auto h_vtxYwrtBpd3fit=temp.Histo1D({"h_vtxYwrtBpd3fit", "vtxYwrtBpd3fit", 500, -.2, .2}, "vtxYwrtBpd3");
    h_vtxYwrtBpd3fit->Fit(&f_vtxYwrtBpd3fit, "w", "",  -.2, .2);
    dd=dd
      .Define("nSigmaVtxXwrtBpd3", [f_vtxXwrtBpd3fit](float z){return nSigma(z, f_vtxXwrtBpd3fit);}, {"vtxXwrtBpd3"})
      .Define("nSigmaVtxYwrtBpd3", [f_vtxYwrtBpd3fit](float z){return nSigma(z, f_vtxYwrtBpd3fit);}, {"vtxYwrtBpd3"})
      .Define("nSigmaVtxXYwrtBpd3", r, {"nSigmaVtxXwrtBpd3", "nSigmaVtxYwrtBpd3"})
      .Define("goodVtxXYwrtBpd3", "nSigmaVtxXYwrtBpd3<=3")
      .Define("goodEvent", "goodEventBeforeVtxCut && goodVtx")
      .Define("goodEvent_before_vtxZ", "goodEventBeforeVtxCut && goodVtxXY")
    ;
    if(dir)
    {
      dir->cd();
      h_vtxXwrtBpd3fit->Write();
      h_vtxYwrtBpd3fit->Write();
    }
  }
  return dd; 
}
