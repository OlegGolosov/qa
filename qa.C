#include "qaUtils.h"
#include "makeDF.h"

//void qa(string in="tree13.root", const char* out="qa13.root", string eventSelection="goodEvent") {
void qa(string in="tree30.root", const char* out="qa30.root", string eventSelection="goodEvent") {
//void qa(string in="tree41.root", const char* out="qa41.root", string eventSelection="goodEvent") {
//void qa(string in="tree_mc.root", const char* out="qa_mc.root", string eventSelection="goodMcEvent") {
  auto c=makeChain(in, "t"); 
  RDataFrame d(*c);
  TFile outFile(out, "recreate");
  auto dd=makeDataFrame(d, &outFile);
  if (eventSelection.size()>0)
    dd=dd.Filter(eventSelection);
  cout << "Number of selected events: " << *(dd.Count()) << endl;

  int nBinsP=500, nBinsEta=500, nBinsY=500, nBinsPhi=500, Mmax=500, nBinsVtxXY=500, nBinsVtxZ=1000, nBinsdEdx=400, nBinsM2=1000, nBinsDca=200;
  float ptMax=5, etaMax=7, yMin=-1, yMax=3, psdEmax=8000, vtxXYmax=2, vtxZmin=-700, vtxZmax=-500, dEdxMax=10, M2Min=-3, M2Max=7, dcaMax=5;
  vector <pair <vector<string>, TH1DModel>> h1common = {
    {{"goodVtxXY", ""},             {"", "", 2, 0, 2}},
    {{"vtxPurity", ""},             {"", "", 110, 0, 1.1}},
    {{"nTracks", ""},               {"", "", Mmax, 0, double(Mmax)}},
    {{"vtxNtrFit", ""},             {"", "", Mmax, 0, double(Mmax)}},
    {{"Mgood", ""},                 {"", "", Mmax, 0, double(Mmax)}},
    {{"Mproton", ""},               {"", "", Mmax, 0, double(Mmax)}},
    {{"Mpionneg", ""},              {"", "", Mmax, 0, double(Mmax)}},
    {{"MgoodMid", ""},              {"", "", Mmax, 0, double(Mmax)}},
    {{"MprotonMid", ""},            {"", "", Mmax, 0, double(Mmax)}},
    {{"MpionnegMid", ""},           {"", "", Mmax, 0, double(Mmax)}},
    {{"psdE", ""},                  {"", "", int(psdEmax/10), 0, psdEmax}},
    {{"psd0E", ""},                 {"", "", int(psdEmax/10), 0, psdEmax}},
    {{"psd1E", ""},                 {"", "", int(psdEmax/10), 0, psdEmax}},
    {{"psd2E", ""},                 {"", "", int(psdEmax/10), 0, psdEmax}},
    {{"psd12E", ""},                {"", "", int(psdEmax/10), 0, psdEmax}},
    {{"trNclustTPC", ""},           {"", "", 250, 0, 250}},
    {{"trNclustTPC", "trGood"},     {"", "", 250, 0, 250}},
    {{"trNclustVTPC1", ""},         {"", "", 250, 0, 250}},
    {{"trNclustVTPC1", "trGood"},   {"", "", 250, 0, 250}},
    {{"trNclustVTPC2", ""},         {"", "", 250, 0, 250}},
    {{"trNclustVTPC2", "trGood"},   {"", "", 250, 0, 250}},
    {{"trNclustVTPC", ""},          {"", "", 250, 0, 250}},
    {{"trNclustVTPC", "trGood"},    {"", "", 250, 0, 250}},
    {{"trNclustMTPC", ""},          {"", "", 250, 0, 250}},
    {{"trNclustMTPC", "trGood"},    {"", "", 250, 0, 250}},
    {{"trNclustPotTPC", ""},        {"", "", 250, 0, 250}},
    {{"trNclustPotTPC", "trGood"},  {"", "", 250, 0, 250}},
    {{"trNclustToPot", ""},         {"", "", 100, 0, 2}},
    {{"trNclustToPot", "trGood"},   {"", "", 100, 0, 2}},
    {{"trNclustToPotCut", ""},      {"", "", 2, 0, 2}},
    {{"trNclustCut", ""},           {"", "", 2, 0, 2}},
    {{"trNdf", ""},                 {"", "", 500, 0, 500}},
    {{"trNdf", "trGood"},           {"", "", 500, 0, 500}},
    {{"trChi2ndf", ""},             {"", "", 1000, 0, 100}},
    {{"trChi2ndf", "trGood"},       {"", "", 1000, 0, 100}},
    {{"trChi2ndfCut", ""},          {"", "", 2, 0, 2}},
    {{"trQp", ""},                  {"", "", nBinsP, -10, 20}},
    {{"trQp", "trGood"},            {"", "", nBinsP, -10, 20}},
    {{"trQoverP", ""},              {"", "", nBinsP, -10, 10}},
    {{"trQoverP", "trGood"},        {"", "", nBinsP, -10, 10}},
    {{"trPx", ""},                  {"", "", nBinsP, -3, 3}},
    {{"trPx", "trGood"},            {"", "", nBinsP, -3, 3}},
    {{"trPy", ""},                  {"", "", nBinsP, -3, 3}},
    {{"trPy", "trGood"},            {"", "", nBinsP, -3, 3}},
    {{"trPz", ""},                  {"", "", nBinsP, 0, 30}},
    {{"trPz", "trGood"},            {"", "", nBinsP, 0, 30}},
    {{"trPt", ""},                 {"", "", nBinsP, 0, ptMax}},
    {{"trPt", "trGood"},           {"", "", nBinsP, 0, ptMax}},
    {{"trPt", "pionneg"},          {"", "", nBinsP, 0, ptMax}},
    {{"trPt", "pionpos"},          {"", "", nBinsP, 0, ptMax}},
    {{"trPt", "proton"},           {"", "", nBinsP, 0, ptMax}},
    {{"trEta", ""},                 {"", "", nBinsEta, 0, etaMax}},
    {{"trEta", "trGood"},           {"", "", nBinsEta, 0, etaMax}},
    {{"trEta", "pionneg"},          {"", "", nBinsEta, 0, etaMax}},
    {{"trEta", "pionpos"},          {"", "", nBinsEta, 0, etaMax}},
    {{"trEta", "proton"},           {"", "", nBinsEta, 0, etaMax}},
    {{"trPhi", ""},                 {"", "", nBinsPhi, -3.15, 3.15}},
    {{"trPhi", "trGood"},           {"", "", nBinsPhi, -3.15, 3.15}},
    {{"trPhi", "pionneg"},          {"", "", nBinsPhi, -3.15, 3.15}},
    {{"trPhi", "pionpos"},          {"", "", nBinsPhi, -3.15, 3.15}},
    {{"trPhi", "proton"},           {"", "", nBinsPhi, -3.15, 3.15}},
    {{"trY", ""},                 {"", "", nBinsY, yMin, yMax}},
    {{"trY", "trGood"},           {"", "", nBinsY, yMin, yMax}},
    {{"trY", "pionneg"},          {"", "", nBinsY, yMin, yMax}},
    {{"trY", "pionpos"},          {"", "", nBinsY, yMin, yMax}},
    {{"trY", "proton"},           {"", "", nBinsY, yMin, yMax}},
    {{"trDcaX", ""}, 	              {"", "", nBinsDca, -dcaMax, dcaMax}},
    {{"trDcaX", "trGood"},          {"", "", nBinsDca, -dcaMax, dcaMax}},
    {{"trDcaY", ""}, 	              {"", "", nBinsDca, -dcaMax, dcaMax}},
    {{"trDcaY", "trGood"},          {"", "", nBinsDca, -dcaMax, dcaMax}},
    {{"trDcaZ", ""}, 	              {"", "", nBinsDca, -0.5*dcaMax, 0.5*dcaMax}},
    {{"trDcaZ", "trGood"},          {"", "", nBinsDca, -0.5*dcaMax, 0.5*dcaMax}},
    {{"trDcaCut", ""},              {"", "", 2, 0, 2}},
    {{"trGood", ""},                {"", "", 2, 0, 2}},
    {{"trGoodNclust", ""},          {"", "", 2, 0, 2}},
    {{"vtxZ", ""},                  {"", "", nBinsVtxZ, vtxZmin, vtxZmax}},
    {{"trQlog20p", ""},             {"", "", 500, -10, 10}},
    {{"trQlog20p", "trGood"},       {"", "", 500, -10, 10}},
  };
  
  vector <pair <vector<string>, TH1DModel>> h1recOnly = {
    {{"adcV1", ""},                 {"", "", 1000, 0, 1000}},
    {{"adcV1p", ""},                {"", "", 100, 0, 100}},
    {{"adcPSD", ""},                {"", "", 1400, 0, 1400}},
//    {{"timesS1_1", ""},             {"", "", int(6e4), -3e4, 3e4}},
    {{"goodS1S2", ""},              {"", "", 2, 0, 2}},
    {{"goodWfa", ""},               {"", "", 2, 0, 2}},
    {{"goodEventBeforeVtxCut", ""}, {"", "", 2, 0, 2}},
    {{"goodBpd", ""},               {"", "", 2, 0, 2}},
    {{"goodVtxXYwrtBpd3", ""},      {"", "", 2, 0, 2}},
    {{"goodEvent", ""},             {"", "", 2, 0, 2}},
    {{"trQdEdx", ""},       	      {"", "", 400, -4, 4}},
    {{"trQdEdx", "trGood"}, 	      {"", "", 400, -4, 4}},
    {{"trQdEdx", "pionneg"},	      {"", "", 400, -4, 4}},
    {{"trQdEdx", "proton"}, 	      {"", "", 400, -4, 4}},
    {{"trTofM2", ""},               {"", "", nBinsM2, M2Min, M2Max}},
    {{"trTofM2", "trGood"},         {"", "", nBinsM2, M2Min, M2Max}},
    {{"trTofM2", "pionneg"},        {"", "", nBinsM2, M2Min, M2Max}},
    {{"trTofM2", "proton"},         {"", "", nBinsM2, M2Min, M2Max}},
  };
  
  vector <pair <vector<string>, TH1DModel>> h1simOnly = {
    {{"goodMcEventBeforeVtxCut", ""}, {"", "", 2, 0, 2}},
    {{"trPrimary", ""}, {"", "", 2, 0, 2}},
    {{"trPrimary", "trDcaCut"}, {"", "", 2, 0, 2}},
    {{"trPrimary", "trGood"}, {"", "", 2, 0, 2}},
  };

  vector <pair <vector<string>, TH2DModel>> h2common = {
    {{"vtxX", "vtxY", ""},                  {"", "", nBinsVtxXY, -vtxXYmax, vtxXYmax, nBinsVtxXY, -vtxXYmax, vtxXYmax}},
    {{"vtxX", "Mgood", ""},                 {"", "", nBinsVtxXY, -vtxXYmax, vtxXYmax, Mmax, 0, double(Mmax)}},
    {{"vtxY", "Mgood", ""},                 {"", "", nBinsVtxXY, -vtxXYmax, vtxXYmax, Mmax, 0, double(Mmax)}},
    {{"vtxZ", "Mgood", ""},                 {"", "", nBinsVtxZ, vtxZmin, vtxZmax, Mmax, 0, double(Mmax)}},
    {{"vtxX", "nTracks", ""},               {"", "", nBinsVtxXY, -vtxXYmax, vtxXYmax, Mmax, 0, double(Mmax)}},
    {{"vtxY", "nTracks", ""},               {"", "", nBinsVtxXY, -vtxXYmax, vtxXYmax, Mmax, 0, double(Mmax)}},
    {{"vtxZ", "nTracks", ""},               {"", "", nBinsVtxZ, vtxZmin, vtxZmax, Mmax, 0, double(Mmax)}},
    {{"vtxX", "vtxNtrFit", ""},             {"", "", nBinsVtxXY, -vtxXYmax, vtxXYmax, Mmax, 0, double(Mmax)}},
    {{"vtxY", "vtxNtrFit", ""},             {"", "", nBinsVtxXY, -vtxXYmax, vtxXYmax, Mmax, 0, double(Mmax)}},
    {{"vtxZ", "vtxNtrFit", ""},             {"", "", nBinsVtxZ, vtxZmin, vtxZmax, Mmax, 0, double(Mmax)}},
    {{"vtxChi2ndf", "Mgood", ""},           {"", "", 120, 0, 12, Mmax, 0, double(Mmax)}},
    {{"vtxZ", "vtxPurity", ""},             {"", "", 500, -597, -587, 110, 0, 1.1}},
    {{"Mgood", "vtxPurity", ""},            {"", "", Mmax, 0, double(Mmax), 100, 0, 1.1}},
    {{"psdE", "vtxPurity", ""},             {"", "", int(psdEmax/10), 0, psdEmax, 100, 0, 1.1}},
    {{"Mgood", "psdE", ""},                 {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mgood", "psd0E", ""},                {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mgood", "psd1E", ""},                {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mgood", "psd2E", ""},                {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mgood", "psd12E", ""},               {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mgood", "psd3E", ""},                {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mgood", "centrality", ""},           {"", "", Mmax, 0, double(Mmax), 20, 0, 100}},
    {{"Mproton", "psdE", ""},               {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mproton", "psd0E", ""},              {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mproton", "psd1E", ""},              {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mproton", "psd2E", ""},              {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mproton", "psd12E", ""},             {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mproton", "psd3E", ""},              {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mproton", "centrality", ""},         {"", "", Mmax, 0, double(Mmax), 20, 0, 100}},
    {{"Mpionneg", "psdE", ""},              {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mpionneg", "psd0E", ""},             {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mpionneg", "psd1E", ""},             {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mpionneg", "psd2E", ""},             {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mpionneg", "psd12E", ""},            {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mpionneg", "psd3E", ""},             {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"Mpionneg", "centrality", ""},        {"", "", Mmax, 0, double(Mmax), 20, 0, 100}},
    {{"MgoodMid", "psdE", ""},              {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MgoodMid", "psd0E", ""},             {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MgoodMid", "psd1E", ""},             {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MgoodMid", "psd2E", ""},             {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MgoodMid", "psd12E", ""},            {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MgoodMid", "psd3E", ""},             {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MgoodMid", "centrality", ""},        {"", "", Mmax, 0, double(Mmax), 20, 0, 100}},
    {{"MprotonMid", "psdE", ""},            {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MprotonMid", "psd0E", ""},           {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MprotonMid", "psd1E", ""},           {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MprotonMid", "psd2E", ""},           {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MprotonMid", "psd12E", ""},          {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MprotonMid", "psd3E", ""},           {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MprotonMid", "centrality", ""},      {"", "", Mmax, 0, double(Mmax), 20, 0, 100}},
    {{"MpionnegMid", "psdE", ""},           {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MpionnegMid", "psd0E", ""},          {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MpionnegMid", "psd1E", ""},          {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MpionnegMid", "psd2E", ""},          {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MpionnegMid", "psd12E", ""},         {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MpionnegMid", "psd3E", ""},          {"", "", Mmax, 0, double(Mmax), int(psdEmax/10), 0, psdEmax}},
    {{"MpionnegMid", "centrality", ""},     {"", "", Mmax, 0, double(Mmax), 20, 0, 100}},
    {{"psd0E", "psd1E", ""},                {"", "", int(psdEmax/10), 0, psdEmax, int(psdEmax/10), 0, psdEmax}},
    {{"psd0E", "psd2E", ""},                {"", "", int(psdEmax/10), 0, psdEmax, int(psdEmax/10), 0, psdEmax}},
    {{"psd0E", "psd3E", ""},                {"", "", int(psdEmax/10), 0, psdEmax, int(psdEmax/10), 0, psdEmax}},
    {{"psd1E", "psd2E", ""},                {"", "", int(psdEmax/10), 0, psdEmax, int(psdEmax/10), 0, psdEmax}},
    {{"psd1E", "psd3E", ""},                {"", "", int(psdEmax/10), 0, psdEmax, int(psdEmax/10), 0, psdEmax}},
    {{"psd2E", "psd3E", ""},                {"", "", int(psdEmax/10), 0, psdEmax, int(psdEmax/10), 0, psdEmax}},
    {{"psdE", "centrality", ""},            {"", "", int(psdEmax/10), 0, psdEmax, 20, 0, 100}},
    {{"psd0E", "centrality", ""},           {"", "", int(psdEmax/10), 0, psdEmax, 20, 0, 100}},
    {{"psd1E", "centrality", ""},           {"", "", int(psdEmax/10), 0, psdEmax, 20, 0, 100}},
    {{"psd12E", "centrality", ""},          {"", "", int(psdEmax/10), 0, psdEmax, 20, 0, 100}},
    {{"psdModId", "psdModE", ""},           {"", "", 244, 1, 244, int(psdEmax/20), 0, 0.5*psdEmax}},
    {{"trDcaX", "trDcaY", ""},              {"", "", nBinsDca, -dcaMax, dcaMax, nBinsDca, -dcaMax, dcaMax}},
    {{"trDcaX", "trDcaY", "trGood"},        {"", "", nBinsDca, -dcaMax, dcaMax, nBinsDca, -dcaMax, dcaMax}},
    {{"trNclust", "trNclustPot", ""},       {"", "", 250, 0, 250, 250, 0, 250}},
    {{"trNclust", "trNclustPot", "trGood"}, {"", "", 250, 0, 250, 250, 0, 250}},
    {{"trEta", "trPt", ""},                 {"", "", nBinsEta, 0, etaMax, nBinsP, 0, ptMax}},
    {{"trEta", "trPt", "trGood"},           {"", "", nBinsEta, 0, etaMax, nBinsP, 0, ptMax}},
    {{"trEta", "trPt", "trMid"},            {"", "", nBinsEta, 0, etaMax, nBinsP, 0, ptMax}},
    {{"trEta", "trPt", "eneg"},          {"", "", nBinsEta, 0, etaMax, nBinsP, 0, ptMax}},
    {{"trEta", "trPt", "pionneg"},          {"", "", nBinsEta, 0, etaMax, nBinsP, 0, ptMax}},
    {{"trEta", "trPt", "kaonneg"},          {"", "", nBinsEta, 0, etaMax, nBinsP, 0, ptMax}},
    {{"trEta", "trPt", "epos"},          {"", "", nBinsEta, 0, etaMax, nBinsP, 0, ptMax}},
    {{"trEta", "trPt", "pionpos"},          {"", "", nBinsEta, 0, etaMax, nBinsP, 0, ptMax}},
    {{"trEta", "trPt", "kaonpos"},          {"", "", nBinsEta, 0, etaMax, nBinsP, 0, ptMax}},
    {{"trEta", "trPt", "proton"},           {"", "", nBinsEta, 0, etaMax, nBinsP, 0, ptMax}},
    {{"trEta", "trPt", "deuteron"},           {"", "", nBinsEta, 0, etaMax, nBinsP, 0, ptMax}},
    {{"trPhi", "trPt", ""},                 {"", "", nBinsPhi, -3.15, 3.15, nBinsP, 0, ptMax}},
    {{"trPhi", "trPt", "trGood"},           {"", "", nBinsPhi, -3.15, 3.15, nBinsP, 0, ptMax}},
    {{"trPhi", "trPt", "trMid"},            {"", "", nBinsPhi, -3.15, 3.15, nBinsP, 0, ptMax}},
    {{"trPhi", "trPt", "pionneg"},          {"", "", nBinsPhi, -3.15, 3.15, nBinsP, 0, ptMax}},
    {{"trPhi", "trPt", "pionpos"},          {"", "", nBinsPhi, -3.15, 3.15, nBinsP, 0, ptMax}},
    {{"trPhi", "trPt", "proton"},          {"", "", nBinsPhi, -3.15, 3.15, nBinsP, 0, ptMax}},
    {{"trPhi", "trEta", ""},                {"", "", nBinsPhi, -3.15, 3.15, nBinsEta, 0, etaMax}},
    {{"trPhi", "trEta", "trGood"},          {"", "", nBinsPhi, -3.15, 3.15, nBinsEta, 0, etaMax}},
    {{"trPhi", "trEta", "trMid"},           {"", "", nBinsPhi, -3.15, 3.15, nBinsEta, 0, etaMax}},
    {{"trPhi", "trEta", "pionneg"},         {"", "", nBinsPhi, -3.15, 3.15, nBinsEta, 0, etaMax}},
    {{"trPhi", "trEta", "proton"},          {"", "", nBinsPhi, -3.15, 3.15, nBinsEta, 0, etaMax}},
    {{"trY", "trPt", "pionneg"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "pionpos"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "proton"},             {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "pionneg_effTr"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "pionpos_effTr"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "proton_effTr"},             {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trYproton", "trPt", "pionpos"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trYpion", "trPt", "proton"},             {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trYproton", "trPt", "deuteron"},             {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trYpion", "trPt", "deuteron"},             {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trPhi", "trY", "pionneg"},           {"", "", nBinsPhi, -3.15, 3.15, nBinsY, yMin, yMax}},
    {{"trPhi", "trY", "pionpos"},           {"", "", nBinsPhi, -3.15, 3.15, nBinsY, yMin, yMax}},
    {{"trPhi", "trY", "proton"},            {"", "", nBinsPhi, -3.15, 3.15, nBinsY, yMin, yMax}},
  };

  vector <pair <vector<string>, TH2DModel>> h2recOnly = {
    {{"adcS1", "adcS2", ""},                {"", "", 100, 100, 200, 120, 130, 250}},
    {{"vtxX", "bpd3x", ""},                 {"", "", nBinsVtxXY, -vtxXYmax, vtxXYmax, nBinsVtxXY, -vtxXYmax, vtxXYmax}},
    {{"vtxY", "bpd3y", ""},                 {"", "", nBinsVtxXY, -vtxXYmax, vtxXYmax, nBinsVtxXY, -vtxXYmax, vtxXYmax}},
    {{"bpd1x", "bpd1y", ""},                {"", "", nBinsVtxXY, -vtxXYmax, vtxXYmax, nBinsVtxXY, -vtxXYmax, vtxXYmax}},
    {{"bpd2x", "bpd2y", ""},                {"", "", nBinsVtxXY, -vtxXYmax, vtxXYmax, nBinsVtxXY, -vtxXYmax, vtxXYmax}},
    {{"bpd3x", "bpd3y", ""},                {"", "", nBinsVtxXY, -vtxXYmax, vtxXYmax, nBinsVtxXY, -vtxXYmax, vtxXYmax}},
    {{"vtxXwrtBpd3", "vtxYwrtBpd3", ""},    {"", "", nBinsVtxXY, -vtxXYmax/10, vtxXYmax/10, nBinsVtxXY, -vtxXYmax/10, vtxXYmax/10}},
    {{"trQlog20p", "trdEdx", ""},           {"", "", 500, -10, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trQlog20p", "trdEdx", "trGood"},     {"", "", 500, -10, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trLog20p", "trdEdx", "trGoodNegWeightNdEdx"},   {"", "", 200, 0, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trLog20p", "trdEdx", "trGoodPosWeightNdEdx"},   {"", "", 200, 0, 10, int(nBinsdEdx*2.5), 0, dEdxMax*2.5}},
    {{"trLog20p", "trdEdx", "deuteron"},    {"", "", 200, 0, 10, int(nBinsdEdx*2.5), 0, dEdxMax*2.5}},
    {{"trLog20p", "trdEdx", "proton"},      {"", "", 200, 0, 10, int(nBinsdEdx*2.5), 0, dEdxMax*2.5}},
    {{"trLog20p", "trdEdx", "pionpos"},     {"", "", 200, 0, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trLog20p", "trdEdx", "kaonpos"},     {"", "", 200, 0, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trLog20p", "trdEdx", "epos"},        {"", "", 200, 0, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trLog20p", "trdEdx", "pionneg"},     {"", "", 200, 0, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trLog20p", "trdEdx", "kaonneg"},     {"", "", 200, 0, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trLog20p", "trdEdx", "eneg"},        {"", "", 200, 0, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trLog20p", "trdEdx", "proton_dEdx_m2"},      {"", "", 200, 0, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trLog20p", "trdEdx", "pionpos_dEdx_m2"},     {"", "", 200, 0, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trLog20p", "trdEdx", "kaonpos_dEdx_m2"},     {"", "", 200, 0, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trLog20p", "trdEdx", "epos_dEdx_m2"},        {"", "", 200, 0, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trLog20p", "trdEdx", "pionneg_dEdx_m2"},     {"", "", 200, 0, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trLog20p", "trdEdx", "kaonneg_dEdx_m2"},     {"", "", 200, 0, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trLog20p", "trdEdx", "eneg_dEdx_m2"},        {"", "", 200, 0, 10, nBinsdEdx, 0, dEdxMax}},
    {{"trdEdx", "trTofM2", ""},             {"", "", nBinsdEdx, 0, dEdxMax, nBinsM2, M2Min, M2Max}},
    {{"trdEdx", "trTofM2", "deuteron"},     {"", "", nBinsdEdx, 0, dEdxMax, nBinsM2, M2Min, M2Max}},
    {{"trdEdx", "trTofM2", "proton"},       {"", "", nBinsdEdx, 0, dEdxMax, nBinsM2, M2Min, M2Max}},
    {{"trdEdx", "trTofM2", "pionpos"},      {"", "", nBinsdEdx, 0, dEdxMax, nBinsM2, M2Min, M2Max}},
    {{"trdEdx", "trTofM2", "kaonpos"},      {"", "", nBinsdEdx, 0, dEdxMax, nBinsM2, M2Min, M2Max}},
    {{"trdEdx", "trTofM2", "epos"},         {"", "", nBinsdEdx, 0, dEdxMax, nBinsM2, M2Min, M2Max}},
    {{"trdEdx", "trTofM2", "pionneg"},      {"", "", nBinsdEdx, 0, dEdxMax, nBinsM2, M2Min, M2Max}},
    {{"trdEdx", "trTofM2", "kaonneg"},      {"", "", nBinsdEdx, 0, dEdxMax, nBinsM2, M2Min, M2Max}},
    {{"trdEdx", "trTofM2", "eneg"},         {"", "", nBinsdEdx, 0, dEdxMax, nBinsM2, M2Min, M2Max}},
    {{"trdEdx", "trTofM2", "trGood"},       {"", "", nBinsdEdx, 0, dEdxMax, nBinsM2, M2Min, M2Max}},
    {{"trdEdx", "trTofM2", "trGoodNegWeightNdEdx"},    {"", "", nBinsdEdx, 0, dEdxMax, nBinsM2, M2Min, M2Max}},
    {{"trdEdx", "trTofM2", "trGoodPosWeightNdEdx"},    {"", "", nBinsdEdx, 0, dEdxMax, nBinsM2, M2Min, M2Max}},
    {{"trP", "trPt", "trGoodPos"},     {"", "", 900, 0, 150, 100, 0, 5}},
    {{"trP", "trPt", "trGoodNeg"},     {"", "", 900, 0, 150, 100, 0, 5}},
    {{"trY", "trPt", "pionneg_dEdx_m2"},    {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "pionpos_dEdx_m2"},    {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "proton_dEdx_m2"},     {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "pionneg_effPid"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "pionpos_effPid"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "proton_effPid"},             {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "pionneg_effTrPid"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "pionpos_effTrPid"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "proton_effTrPid"},             {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
  };

  vector <pair <vector<string>, TH2DModel>> h2simOnly = {
    {{"simY", "simPt", "sim_pionneg_prim"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"simY", "simPt", "sim_pionpos_prim"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"simY", "simPt", "sim_proton_prim"},             {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"simPhi", "simY", "sim_pionneg_prim"},            {"", "", nBinsPhi, -3.15, 3.15, nBinsY, yMin, yMax}},
    {{"simPhi", "simY", "sim_pionpos_prim"},            {"", "", nBinsPhi, -3.15, 3.15, nBinsY, yMin, yMax}},
    {{"simPhi", "simY", "sim_proton_prim"},             {"", "", nBinsPhi, -3.15, 3.15, nBinsY, yMin, yMax}},
    {{"simPhi", "simPt", "sim_pionneg_prim"},            {"", "", nBinsPhi, -3.15, 3.15, nBinsP, 0, ptMax}},
    {{"simPhi", "simPt", "sim_pionpos_prim"},            {"", "", nBinsPhi, -3.15, 3.15, nBinsP, 0, ptMax}},
    {{"simPhi", "simPt", "sim_proton_prim"},             {"", "", nBinsPhi, -3.15, 3.15, nBinsP, 0, ptMax}},
  };

  vector <pair <vector<string>, TProfile1DModel>> p1common = {
//    {{"vtxZ", "nSigmaVtxZ", ""},        {"", "", nBinsVtxXY, vtxZmin, vtxZmax}},
  };

  vector <pair <vector<string>, TProfile1DModel>> p1recOnly = {
//    {{"vtxZ", "nSigmaVtxZ", ""},        {"", "", nBinsVtxXY, vtxZmin, vtxZmax}},
  };

  vector <pair <vector<string>, TProfile1DModel>> p1simOnly = {
  };

  vector <pair <vector<string>, TProfile2DModel>> p2common = {
    {{"psdModX", "psdModY", "psdModE", ""}, {"", "", 100, -150, 150, 100, -150, 150}},
    {{"psdModX", "psdModY", "psdModE", "psdModSub0"}, {"", "", 100, -150, 150, 100, -150, 150}},
    {{"psdModX", "psdModY", "psdModE", "psdModSub1"}, {"", "", 100, -150, 150, 100, -150, 150}},
    {{"psdModX", "psdModY", "psdModE", "psdModSub2"}, {"", "", 100, -150, 150, 100, -150, 150}},
    {{"psdModX", "psdModY", "psdModE", "psdModSub3"}, {"", "", 100, -150, 150, 100, -150, 150}},
    {{"trY", "trPt", "pionneg_effTr", "pionneg"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "pionpos_effTr", "pionpos"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "proton_effTr", "proton"},             {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
  };

  vector <pair <vector<string>, TProfile2DModel>> p2recOnly = {
    {{"trY", "trPt", "pionneg_effPid", "pionneg"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "pionpos_effPid", "pionpos"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "proton_effPid", "proton"},             {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "pionneg_effTrPid", "pionneg"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "pionpos_effTrPid", "pionpos"},            {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
    {{"trY", "trPt", "proton_effTrPid", "proton"},             {"", "", nBinsY, yMin, yMax, nBinsP, 0, ptMax}},
//    {{"trP", "trdEdx", "purityProton", "trGoodPos"},      {"", "", 1000, 0.05, 150, 100, 0, 5}},
//    {{"trP", "trdEdx", "purityPion", "trGoodNeg"},        {"", "", 1000, 0.05, 150, 100, 0, 5}},
//    {{"trP", "trdEdx", "purityPion", "trGoodPos"},        {"", "", 1000, 0.05, 150, 100, 0, 5}},
  };

  vector <pair <vector<string>, TProfile2DModel>> p2simOnly = {
  };

  vector <pair <vector<string>, TH3DModel>> h3common = {
//    {{"trY", "trPt", "trPhi", "pionpos"},      {"", "", nBinsY/5, yMin, yMax, nBinsP/5, 0, ptMax, nBinsPhi/5, -3.15, 3.15}},
//    {{"trY", "trPt", "trPhi", "pionneg"},      {"", "", nBinsY/5, yMin, yMax, nBinsP/5, 0, ptMax, nBinsPhi/5, -3.15, 3.15}},
//    {{"trY", "trPt", "trPhi", "proton"},      {"", "", nBinsY/5, yMin, yMax, nBinsP/5, 0, ptMax, nBinsPhi/5, -3.15, 3.15}},
  };

  vector <pair <vector<string>, TH3DModel>> h3recOnly = {
//    {{"trP", "trPt", "trdEdx", "trGoodPosWeightNdEdx"},      {"", "", 400, 0, 20, 200, 0, 5, nBinsdEdx, 0, dEdxMax}},
//    {{"trP", "trPt", "trdEdx", "trGoodNegWeightNdEdx"},      {"", "", 400, 0, 20, 200, 0, 5, int(nBinsdEdx*2.5), 0, dEdxMax*2.5}},
//    {{"trLog20p", "trPt", "trdEdx", "trGoodPosWeightNdEdx"},      {"", "", 400, 0, 10, 200, 0, 5, nBinsdEdx, 0, dEdxMax}},
//    {{"trLog20p", "trPt", "trdEdx", "trGoodNegWeightNdEdx"},      {"", "", 400, 0, 10, 200, 0, 5, int(nBinsdEdx*2.5), 0, dEdxMax*2.5}},
  };

  vector <pair <vector<string>, TH3DModel>> h3simOnly = {
//    {{"simY", "simPt", "simPhi", "sim_pionneg_prim"},            {"", "", nBinsY/5, yMin, yMax, nBinsP/5, 0, ptMax, nBinsPhi/5, -3.15, 3.15}},
//    {{"simY", "simPt", "simPhi", "sim_pionpos_prim"},            {"", "", nBinsY/5, yMin, yMax, nBinsP/5, 0, ptMax, nBinsPhi/5, -3.15, 3.15}},
//    {{"simY", "simPt", "simPhi", "sim_proton_prim"},             {"", "", nBinsY/5, yMin, yMax, nBinsP/5, 0, ptMax, nBinsPhi/5, -3.15, 3.15}},
  };

  saveHists(dd, h1common,  outFile);
  saveHists(dd, p1common,  outFile);
  saveHists(dd, h2common,  outFile);
  saveHists(dd, p2common,  outFile);
  saveHists(dd, h3common,  outFile);

  if (!isSimulation)
  {
    saveHists(dd, h1recOnly,  outFile);
    saveHists(dd, p1recOnly,  outFile);
    saveHists(dd, h2recOnly,  outFile);
    saveHists(dd, p2recOnly,  outFile);
    saveHists(dd, h3recOnly,  outFile);
//    #include "histsKIT.h"
//    saveHists(dd, h1recOnlyKIT,  outFile);
//    saveHists(dd, h2recOnlyKIT,  outFile);
  }
  if (isSimulation)
  {
    saveHists(dd, h1simOnly,  outFile);
    saveHists(dd, p1simOnly,  outFile);
    saveHists(dd, h2simOnly,  outFile);
    saveHists(dd, p2simOnly,  outFile);
    saveHists(dd, h3simOnly,  outFile);
  }

// #include "histsDT.h"
//  saveHists(dd, h1commonDT,  *outFile.mkdir("reco_info"));
  outFile.Close();
}
