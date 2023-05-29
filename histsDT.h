vector <pair <vector<string>, TH1DModel>> h1commonDT = {
  {{"timesS1_1", ""},             {"hWFABeamHits", "", int(6e4), -3e4, 3e4}},
  {{"adcS1", ""},                 {"hS1", "", 100, 100, 200}},
  {{"adcS2", ""},                 {"hS2", "", 120, 130, 250}},
  {{"vtxX", ""},                  {"hVertexPosX", "", nBinsVtxXY, -vtxXYmax, vtxXYmax}},
  {{"vtxY", ""},                  {"hVertexPosY", "", nBinsVtxXY, -vtxXYmax, vtxXYmax}},
  {{"vtxZ", ""},                  {"hVertexPosZ", "", nBinsVtxZ, vtxZmin, vtxZmax}},
  {{"psdE", ""},                  {"hE", "", int(psdEmax/10), 0, psdEmax}},

//  {{"nTracks", ""},               {"hMreco", "", Mmax, 0, double(Mmax)}},
//  {{"trNdf", ""},                 {"hNdf", "", 500, 0, 500}},
//  {{"trCharge", ""},              {"hCharge", "", 5, -2, 2}},
//  {{"trChi2ndf", ""},             {"hChi2Ndf", "", 1000, 0, 100}},
//  {{"trNclustVTPC1", ""},         {"hNhits_VTPC1", "", 250, 0, 250}},
//  {{"trNclustVTPC2", ""},         {"hNhits_VTPC2", "", 250, 0, 250}},
//  {{"trNclustMTPC", ""},          {"hNhits_MTPC", "", 250, 0, 250}},
//  {{"trNclustTPC", ""},           {"hNhits_allTPC", "", 250, 0, 250}},
//  {{"trNclustPotVTPC1", ""},      {"hNhitsPot_VTPC1", "", 250, 0, 250}},
//  {{"trNclustPotVTPC2", ""},      {"hNhitsPot_VTPC2", "", 250, 0, 250}},
//  {{"trNclustPotMTPC", ""},       {"hNhitsPot_MTPC", "", 250, 0, 250}},
//  {{"trNclustPotTPC", ""},        {"hNhitsPot_allTPC", "", 250, 0, 250}},
//  {{"trPt", ""},                  {"hTrackPt", "", nBinsP, 0, ptMax}},
//  {{"trEta", ""},                 {"hTrackEta", "", nBinsEta, 0, etaMax}},
//  {{"trPhi", ""},                 {"hTrackPhi", "", nBinsPhi, -3.15, 3.15}},
//  {{"trDcaX", ""},                {"hTrackDCAX", "", nBinsPhi, -3.15, 3.15}},
//  {{"trDcaY", ""},                {"hTrackDCAY", "", nBinsPhi, -3.15, 3.15}},

  {{"Mgood", ""},                 {"hMreco", "", Mmax, 0, double(Mmax)}},
  {{"trNdf", "trGood"},           {"hNdf", "", 500, 0, 500}},
  {{"trCharge", "trGood"},        {"hCharge", "", 5, -2, 2}},
  {{"trChi2ndf", "trGood"},       {"hChi2Ndf", "", 1000, 0, 100}},
  {{"trNclustVTPC1", "trGood"},   {"hNhits_VTPC1", "", 250, 0, 250}},
  {{"trNclustVTPC2", "trGood"},   {"hNhits_VTPC2", "", 250, 0, 250}},
  {{"trNclustMTPC", "trGood"},    {"hNhits_MTPC", "", 250, 0, 250}},
  {{"trNclust", "trGood"},        {"hNhits_allTPC", "", 250, 0, 250}},
  {{"trNclustPotVTPC1", "trGood"},{"hNhitsPot_VTPC1", "", 250, 0, 250}},
  {{"trNclustPotVTPC2", "trGood"},{"hNhitsPot_VTPC2", "", 250, 0, 250}},
  {{"trNclustPotMTPC", "trGood"}, {"hNhitsPot_MTPC", "", 250, 0, 250}},
  {{"trNclustPotTPC", "trGood"},  {"hNhitsPot_allTPC", "", 250, 0, 250}},
  {{"trPt", "trGood"},            {"hTrackPt", "", nBinsP, 0, ptMax}},
  {{"trEta", "trGood"},           {"hTrackEta", "", nBinsEta, 0, etaMax}},
  {{"trPhi", "trGood"},           {"hTrackPhi", "", nBinsPhi, -3.15, 3.15}},
  {{"trDcaX", "trGood"},          {"hTrackDCAX", "", nBinsPhi, -3.15, 3.15}},
  {{"trDcaY", "trGood"},          {"hTrackDCAY", "", nBinsPhi, -3.15, 3.15}},
  };
