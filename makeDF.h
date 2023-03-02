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

filteredDF makeDataFrame (definedDF &d, TDirectory *dir=nullptr)
{
  int pBeam=round(*d.Range(0,0).Mean("beamPz"));
  float mp=0.938;
  float eBeam=sqrt(pBeam*pBeam+mp*mp);
  float yBeam=0.25*log((eBeam+pBeam)/(eBeam-pBeam));
  printf("pBeam = %i \t yBeam = %f \n", pBeam, yBeam);

  TFile fPid(Form("pid_pbpb%i_graphical.root", pBeam));
  auto pidGetter=(Pid::Getter*)fPid.Get("pid_getter");
  float pidPurity=0.9;
  TDatabasePDG::Instance(); // for some reason needed at lxplus

  auto getPid=[pidGetter, pidPurity](RVec<float> p, RVec<float> dEdx)
  {
    RVec<int> pids;
    for(int i=0;i<p.size();i++)
      pids.push_back(pidGetter->GetPid(p.at(i), dEdx.at(i), pidPurity)); 
    return pids; 
  };

  TFile fCentr(Form("centrality_pbpb%i_hE.root", pBeam));
  auto centrGetter=(Centrality::Getter*)fCentr.Get("centr_getter_1D");

  auto getCentrality=[centrGetter](float E)
  {
    return centrGetter->GetCentrality(E); 
  };

  string pileUpCut;
  if(pBeam==13)
    pileUpCut="psdE<3800*(1-Mgood/230.)";
  else if(pBeam==30)
    pileUpCut="psdE<8000*(1-pow(Mgood/400., 2))";
  
  cout << "Number of events: " << *(d.Count()) << endl;
  auto dd=d
    .Filter("trigT4||trigT2")
//    .Filter("trigT4")
    .Define("goodWfa","Sum(abs(timesS1_1)<4000)==1")
    .Define("goodS1S2", "fabs(adcS1-150)<35 && fabs(adcS2-195)<55 && adcS1*1.258+adcS2<390")
    .Define("psdModPhi", "atan2(psdModY,psdModX)")
    .Define("psdModSub0wo45", "psdModId==6 || psdModId==7 || psdModId==10 || psdModId==11")
    .Define("psdModSub1wo45", "psdModId<=16")
    .Define("psdModSub0", "psdModId==6 || psdModId==7 || psdModId==10 || psdModId==11 || psdModId==45")
    .Define("psdModSub1", "psdModId<=16 || psdModId==45")
    .Define("psdModSub2", "17<=psdModId && psdModId<=28")
    .Define("psdModSub3", "29<=psdModId && psdModId<=44")
    .Define("psdE", "Sum(psdModE)")
    .Define("psd0E", "Sum(psdModE*psdModSub0)")
    .Define("psd1E", "Sum(psdModE*psdModSub1)")
    .Define("psd2E", "Sum(psdModE*psdModSub2)")
    .Define("psd12E", "psd1E+psd2E")
    .Define("psd3E", "Sum(psdModE*psdModSub3)")
    .Define("centrality", getCentrality, {"psdE"})
    .Define("trId", "RVec<int> id; for(int i=0;i<trPt.size();i++)id.push_back(i); return id;")
    .Define("trMom", "RVec<TVector3> mom; for(int i=0;i<nTracks;i++){TVector3 v;v.SetPtEtaPhi(trPt[i], trEta[i], trPhi[i]);mom.push_back(v);} return mom;")
    .Define("trPz", "RVec<float> pz; for(auto &mom:trMom)pz.push_back(mom.Pz()); return pz;")
    .Define("trP", "RVec<float> trP; for(auto &mom:trMom)trP.push_back(mom.Mag()); return trP;")
    .Define("trQp", "trCharge*trP")
    .Define("trQlog20p", "trCharge*log(20*trP)")
    .Define("trNhitsToPot", "RVec<float> ratio(nTracks, -10); for(int i=0;i<nTracks;i++)if(trNhitsPot.at(i)>0)ratio.at(i)=(float)trNhits.at(i)/trNhitsPot.at(i); return ratio;")
    .Define("trNhitsVTPC", "trNhitsVTPC1+trNhitsVTPC2")
    .Define("trNhitsCut", "trNhitsVTPC>15 && trNhits>30")
    .Define("trNhitsToPotCut", "trNhitsToPot>0.55 && trNhitsToPot<1.1")
    .Define("trDcaCut", "abs(trDcaX)<2 && abs(trDcaY)<1")
    .Define("trGoodNhits", "trNhitsCut && trNhitsToPotCut")
    .Define("trGood", "trNhitsCut && trNhitsToPotCut && trDcaCut")
    .Define("trTofM2", "RVec<float> trM2(trId.size(), -999); for(int i=0;i<tofHitM2.size();i++) if(tofHitTrackId.at(i)>=0) trM2.at(tofHitTrackId.at(i))=tofHitM2.at(i); return trM2;")
    .Define("trPid", getPid, {"trQp", "trdEdx"})
    .Define("trM", "RVec<float> m; for(auto& pid:trPid)m.push_back(TDatabasePDG::Instance()->GetParticle(pid)->Mass()); return m;")
    .Define("trMom4", "RVec<TLorentzVector> mom4; for(int i=0;i<nTracks;i++){TLorentzVector v;v.SetPtEtaPhiM(trPt[i], trEta[i], trPhi[i], trM[i]);mom4.push_back(v);} return mom4;")
    .Define("trY", [yBeam](RVec<TLorentzVector> mom4){RVec<float> y; for(auto &mom:mom4)y.push_back(mom.Rapidity()-yBeam); return y;}, {"trMom4"})
    .Define("pionneg", "trGood && trPid==-211")
    .Define("proton", "trGood && trPid==2212")
    .Define("trMid", "trPt>0.1 && trEta<3.")
    .Define("Mgood", "Sum(trGood)")
    .Define("Mproton", "Sum(proton)")
    .Define("Mpionneg", "Sum(pionneg)")
    .Define("MgoodMid", "Sum(trGood*trMid)")
    .Define("MprotonMid", "Sum(proton*trMid)")
    .Define("MpionnegMid", "Sum(pionneg*trMid)")
    .Define("vtxChi2ndf", "vtxChi2/vtxNdf")
    .Define("vtxXwrtBpd3", "vtxX-bpd3x")
    .Define("vtxYwrtBpd3", "vtxY-bpd3y")
    .Define("pileUpCut", pileUpCut)
    .Define("goodEventBeforeVtxCut", "goodS1S2 && goodWfa && pileUpCut && vtxFitPerfect")
    .Define("goodVtxOld", "-594.0<vtxZ && vtxZ<-590.0 && -0.74<vtxX && vtxX<0.131 && -0.66<vtxY && vtxY<0.022")
    .Define("goodEventOld", "goodEventBeforeVtxCut && goodVtxOld")
  ;
  
  auto temp=dd.Filter("goodEventBeforeVtxCut");

  float vtxZnominal=-591.9, vtxZfitOffset=0.2;
  TF1 f_vtxXfit("f_vtxXfit","gaus", -2, 2);
  TF1 f_vtxYfit("f_vtxYfit","gaus", -2, 2);
  TF1 f_vtxZfit("f_vtxZfit","gaus",vtxZnominal-vtxZfitOffset,vtxZnominal+vtxZfitOffset);
  TF1 f_vtxXwrtBpd3fit("f_vtxXwrtBpd3fit","gaus", -.2, .2);
  TF1 f_vtxYwrtBpd3fit("f_vtxYwrtBpd3fit","gaus", -.2, .2);
  auto h_vtxXfit=temp.Histo1D({"h_vtxXfit", "vtxXfit", 500, -2, 2}, "vtxX");
  h_vtxXfit->Fit(&f_vtxXfit, "w", "",  -2, 2);
  auto h_vtxYfit=temp.Histo1D({"h_vtxYfit", "vtxYfit", 500, -2, 2}, "vtxY");
  h_vtxYfit->Fit(&f_vtxYfit, "w", "",  -2, 2);
  auto h_vtxZfit=temp.Histo1D({"h_vtxZfit", "vtxZfit", 500, vtxZnominal-2,vtxZnominal+2}, "vtxZ");
  h_vtxZfit->Fit(&f_vtxZfit, "w", "", vtxZnominal-vtxZfitOffset, vtxZnominal+vtxZfitOffset);
  auto h_vtxXwrtBpd3fit=temp.Histo1D({"h_vtxXwrtBpd3fit", "vtxXwrtBpd3fit", 500, -.2, .2}, "vtxXwrtBpd3");
  h_vtxXwrtBpd3fit->Fit(&f_vtxXwrtBpd3fit, "w", "",  -.2, .2);
  auto h_vtxYwrtBpd3fit=temp.Histo1D({"h_vtxYwrtBpd3fit", "vtxYwrtBpd3fit", 500, -.2, .2}, "vtxYwrtBpd3");
  h_vtxYwrtBpd3fit->Fit(&f_vtxYwrtBpd3fit, "w", "",  -.2, .2);
  if(dir)
  {
    dir->cd();
    h_vtxXfit->Write();
    h_vtxYfit->Write();
    h_vtxZfit->Write();
    h_vtxXwrtBpd3fit->Write();
    h_vtxYwrtBpd3fit->Write();
  }
  
  dd=dd
    .Define("nSigmaVtxX", [f_vtxXfit](float z){return nSigma(z, f_vtxXfit);}, {"vtxX"})
    .Define("nSigmaVtxY", [f_vtxYfit](float z){return nSigma(z, f_vtxYfit);}, {"vtxY"})
    .Define("nSigmaVtxZ", [f_vtxZfit](float z){return nSigma(z, f_vtxZfit);}, {"vtxZ"})
    .Define("nSigmaVtxXwrtBpd3", [f_vtxXwrtBpd3fit](float z){return nSigma(z, f_vtxXwrtBpd3fit);}, {"vtxXwrtBpd3"})
    .Define("nSigmaVtxYwrtBpd3", [f_vtxYwrtBpd3fit](float z){return nSigma(z, f_vtxYwrtBpd3fit);}, {"vtxYwrtBpd3"})
    .Define("nSigmaVtxXY", r, {"nSigmaVtxX", "nSigmaVtxY"})
    .Define("nSigmaVtxXYwrtBpd3", r, {"nSigmaVtxXwrtBpd3", "nSigmaVtxYwrtBpd3"})
    .Define("goodVtxZ", "nSigmaVtxZ<=11.2")
    .Define("goodVtxXY", "nSigmaVtxXY<=3")
    .Define("goodVtxXYwrtBpd3", "nSigmaVtxXYwrtBpd3<=3")
    .Define("goodVtxPos", "goodVtxZ && goodVtxXY")
    .Define("goodVtx", "vtxFitPerfect && goodVtxPos")
    .Define("goodEvent", "goodEventBeforeVtxCut && goodVtx")
  ;
  return dd; 
}
