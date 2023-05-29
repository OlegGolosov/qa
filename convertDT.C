#include "qaUtils.h"

class DataTreeEvent;

void convertDT(string in="dt.root", const char* out="treeDT.root") {
  auto c=makeChain(in, "DataTree");
  RDataFrame d(*c, {"DTEvent"});
  //RDataFrame d("DataTree", in);
  auto dd=d
//    .Range(0,100)
    .Define("runId", "fRunId")
    .Define("eventId", "fEventId")
    .Define("vtxNdf", "return 1.;")
    .Define("vtxChi2", "return 1.;")
    .Define("vtxX", "(float)DTEvent.GetVertexPositionComponent(0,0)")
    .Define("vtxY", "(float)DTEvent.GetVertexPositionComponent(1,0)")
    .Define("vtxZ", "(float)DTEvent.GetVertexPositionComponent(2,0)")
    .Define("vtxFitPerfect", "DTEvent.GetVertexQuality(0)==1?true:false")
    .Define("nTracks", "fNumberOfVertexTracks")
    .Define("tofNhits", "fNumberOfTOFHits")
    .Define("psdNmod", "fNumberOfPSDModules")
    .Define("tracks", "RVec<DataTreeTrack>t;for(int i=0;i<nTracks;i++) t.push_back(*DTEvent.GetVertexTrack(i)); return t;")
    .Define("tofHits", "RVec<DataTreeTOFHit>t;for(int i=0;i<tofNhits;i++) t.push_back(*DTEvent.GetTOFHit(i)); return t;")
    .Define("psdModules", "RVec<DataTreePSDModule>t;for(int i=0;i<psdNmod;i++) t.push_back(*DTEvent.GetPSDModule(i)); return t;")
    .Define("trPt", "RVec<float> v; for (auto &t:tracks) v.push_back(t.GetPt()); return v;")
    .Define("trEta", "RVec<float> v; for (auto &t:tracks) v.push_back(t.GetEta()); return v;")
    .Define("trPhi", "RVec<float> v; for (auto &t:tracks) v.push_back(t.GetPhi()); return v;")
    .Define("trDcaX", "RVec<float> v; for (auto &t:tracks) v.push_back(t.GetDCAComponent(0)); return v;")
    .Define("trDcaY", "RVec<float> v; for (auto &t:tracks) v.push_back(t.GetDCAComponent(1)); return v;")
    .Define("trDcaZ", "RVec<float> v; for (auto &t:tracks) v.push_back(t.GetDCAComponent(2)); return v;")
    .Define("trChi2", "RVec<float> v; for (auto &t:tracks) v.push_back(t.GetChi2()); return v;")
    .Define("trNdf", "RVec<int> v; for (auto &t:tracks) v.push_back(t.GetNDF()); return v;")
    .Define("trCharge", "RVec<int> v; for (auto &t:tracks) v.push_back(t.GetCharge()); return v;")
    .Define("trNhitsVTPC1", "RVec<int> v; for (auto &t:tracks) v.push_back(t.GetNumberOfHits(EnumTPC::kVTPC1)); return v;")
    .Define("trNhitsVTPC2", "RVec<int> v; for (auto &t:tracks) v.push_back(t.GetNumberOfHits(EnumTPC::kVTPC2)); return v;")
    .Define("trNhitsMTPC", "RVec<int> v; for (auto &t:tracks) v.push_back(t.GetNumberOfHits(EnumTPC::kMTPC)); return v;")
    .Define("trNhits", "RVec<int> v; for (auto &t:tracks) v.push_back(t.GetNumberOfHits(EnumTPC::kTPCAll)); return v;")
    .Define("trNhitsPotVTPC1", "RVec<int> v; for (auto &t:tracks) v.push_back(t.GetNumberOfHitsPotential(EnumTPC::kVTPC1)); return v;")
    .Define("trNhitsPotVTPC2", "RVec<int> v; for (auto &t:tracks) v.push_back(t.GetNumberOfHitsPotential(EnumTPC::kVTPC2)); return v;")
    .Define("trNhitsPotMTPC", "RVec<int> v; for (auto &t:tracks) v.push_back(t.GetNumberOfHitsPotential(EnumTPC::kMTPC)); return v;")
    .Define("trNhitsPot", "RVec<int> v; for (auto &t:tracks) v.push_back(t.GetNumberOfHitsPotential(EnumTPC::kTPCAll)); return v;")
    .Define("trNdEdxClustVTPC1", "RVec<int> v; for (auto &t:tracks) v.push_back(t.GetNumberOfdEdxClusters(EnumTPC::kVTPC1)); return v;")
    .Define("trNdEdxClustVTPC2", "RVec<int> v; for (auto &t:tracks) v.push_back(t.GetNumberOfdEdxClusters(EnumTPC::kVTPC2)); return v;")
    .Define("trNdEdxClustMTPC", "RVec<int> v; for (auto &t:tracks) v.push_back(t.GetNumberOfdEdxClusters(EnumTPC::kMTPC)); return v;")
    .Define("trNdEdxClust", "RVec<int> v; for (auto &t:tracks) v.push_back(t.GetNumberOfdEdxClusters(EnumTPC::kTPCAll)); return v;")
    .Define("trdEdxVTPC1", "RVec<float> v; for (auto &t:tracks) v.push_back(t.GetdEdx(EnumTPC::kVTPC1)); return v;")
    .Define("trdEdxVTPC2", "RVec<float> v; for (auto &t:tracks) v.push_back(t.GetdEdx(EnumTPC::kVTPC2)); return v;")
    .Define("trdEdxMTPC", "RVec<float> v; for (auto &t:tracks) v.push_back(t.GetdEdx(EnumTPC::kMTPC)); return v;")
    .Define("trdEdx", "RVec<float> v; for (auto &t:tracks) v.push_back(t.GetdEdx(EnumTPC::kTPCAll)); return v;")
    .Define("tofHitStatus", "RVec<int> v; for (auto &h:tofHits) v.push_back(h.GetStatus()); return v;")
    .Define("tofHitCharge", "RVec<int> v; for (auto &h:tofHits) v.push_back(h.GetCharge()); return v;")
    .Define("tofHitTrackId", "RVec<int> v; for (auto &h:tofHits) v.push_back(h.GetFirstRecoTrackId()); return v;")
    .Define("tofHitX", "RVec<float> v; for (auto &h:tofHits) v.push_back(h.GetX()); return v;")
    .Define("tofHitY", "RVec<float> v; for (auto &h:tofHits) v.push_back(h.GetY()); return v;")
    .Define("tofHitZ", "RVec<float> v; for (auto &h:tofHits) v.push_back(h.GetZ()); return v;")
    .Define("tofHitTime", "RVec<float> v; for (auto &h:tofHits) v.push_back(h.GetTime()); return v;")
    .Define("tofHitLength", "RVec<float> v; for (auto &h:tofHits) v.push_back(h.GetPathLength()); return v;")
    .Define("tofHitM2", "RVec<float> v; for (auto &h:tofHits) v.push_back(h.GetSquaredMass()); return v;")
    .Define("tofHitM2err", "RVec<float> v; for (auto &h:tofHits) v.push_back(h.GetSquaredMassError()); return v;")
    .Define("psdPosX", "(float)DTEvent.GetPsdPositionX()")
    .Define("psdPosY", "(float)DTEvent.GetPsdPositionY()")
    .Define("psdPosZ", "(float)DTEvent.GetPsdPositionZ()")
    .Define("psdModId", "RVec<int> v; for (auto &m:psdModules) v.push_back(m.GetId()+1); return v;")
    .Define("psdModX", "RVec<float> v; for (auto &m:psdModules) v.push_back(m.GetPositionComponent(0)); return v;")
    .Define("psdModY", "RVec<float> v; for (auto &m:psdModules) v.push_back(m.GetPositionComponent(1)); return v;")
    .Define("psdModE", "RVec<float> v; for (auto &m:psdModules) v.push_back(m.GetEnergy()); return v;")
    .Define("beamStatus", "DTEvent.GetBeamStatus()")
    .Define("trigT1", "(float)DTEvent.GetTrigger(EnumTrigger::kT1)->GetIsFired()")
    .Define("trigT2", "(float)DTEvent.GetTrigger(EnumTrigger::kT2)->GetIsFired()")
    .Define("trigT4", "(float)DTEvent.GetTrigger(EnumTrigger::kT4)->GetIsFired()")
    .Define("adcS1", "(float)DTEvent.GetTrigger(EnumTrigger::kS1)->GetSignal()")
    .Define("adcS2", "(float)DTEvent.GetTrigger(EnumTrigger::kS2)->GetSignal()")
    .Define("adcS3", "(float)DTEvent.GetTrigger(EnumTrigger::kS3)->GetSignal()")
    .Define("adcV1", "(float)DTEvent.GetTrigger(EnumTrigger::kV1)->GetSignal()")
    .Define("adcV1p", "(float)DTEvent.GetTrigger(EnumTrigger::kV1p)->GetSignal()")
    .Define("adcPSD", "(float)DTEvent.GetTrigger(EnumTrigger::kPSD)->GetSignal()")
    .Define("bpd1x", "(float)DTEvent.GetBPD(EnumBPD::kBPD1)->GetPositionComponent(0)")
    .Define("bpd2x", "(float)DTEvent.GetBPD(EnumBPD::kBPD2)->GetPositionComponent(0)")
    .Define("bpd3x", "(float)DTEvent.GetBPD(EnumBPD::kBPD3)->GetPositionComponent(0)")
    .Define("bpd1y", "(float)DTEvent.GetBPD(EnumBPD::kBPD1)->GetPositionComponent(1)")
    .Define("bpd2y", "(float)DTEvent.GetBPD(EnumBPD::kBPD2)->GetPositionComponent(1)")
    .Define("bpd3y", "(float)DTEvent.GetBPD(EnumBPD::kBPD3)->GetPositionComponent(1)")
    .Define("beamPx", "(float)DTEvent.GetBeamMomentumX()")
    .Define("beamPy", "(float)DTEvent.GetBeamMomentumY()")
    .Define("beamPz", "(float)DTEvent.GetBeamMomentumZ()")
    .Define("nTimesS1_1", "(int)DTEvent.GetWFA(EnumWFA::kS1_1)->GetNHits()")
    .Define("nTimesT4", "(int)DTEvent.GetWFA(EnumWFA::kT4)->GetNHits()")
    .Define("timesS1_1", "RVec<float> v; for(auto& t:DTEvent.GetWFA(EnumWFA::kS1_1)->GetTimesWFA()) v.push_back(t); return v;")
    .Define("timesT4", "RVec<float> v; for(auto& t:DTEvent.GetWFA(EnumWFA::kT4)->GetTimesWFA()) v.push_back(t); return v;")
  ;
  
  dd.Foreach([](const ULong64_t entry){if (entry % 100 == 0) cout << "\r" << entry;}, {"rdfentry_"}); // progress display 
  cout << endl;

  vector<string> definedNames;
  vector<string> toExclude={"tracks", "tofHits", "psdModules"};
  for (auto& definedName:dd.GetDefinedColumnNames())
  {
    bool exclude=false;
    for (auto &nameToExclude:toExclude)
      if (definedName==nameToExclude)
        exclude=true;
    if (!exclude)
      definedNames.push_back(definedName);
  }
  dd.Snapshot("t",out, definedNames);
}
