vector <TFile*> makeFileVector (string &fileListPath);
vector <string> buildObjectList (vector <string> &histNames, TDirectory *folder, int depth=1);
void filterObjectList (vector <string> &objectList, string &pattern);
vector <float> getRunIds(vector<TFile*> &files, string &pattern);
void plot(vector <TFile*> &files, vector<string> &histnames, vector<float> &xAxis, const char *xAxisTitle, bool isTimeGraph, TDirectoryFile *outDir);

void plotTrends(string fileListPath="qa.list", const char *outFilePath="trends.root", string pattern=".*")
{
  vector <TFile*> files=makeFileVector(fileListPath);
  vector <string> histNames;
  buildObjectList(histNames, files.at(0));
  vector <float> runIds;
  vector <float> runTimes;
  runIds = getRunIds(files, ".*(\\d+)\\.qa\\.root");
  //runTimes = getRunTimes(files);
  filterObjectList(histNames, pattern);
  auto outFile = new TFile(outFilePath, "recreate");
  plot(files, histNames, runIds, "runId", false, outFile);
  if (runTimes.size()>0)
    plot(files, histNames, runTimes, "time", true, outFile);
  outFile->Close();
}

void plot(vector <TFile*> &files, vector<string> &histNames, vector<float> &xAxis, const char *xAxisTitle, bool isTimeGraph, TDirectoryFile *outDir)
{
  printf ("\nPlotting %s trends:\n", xAxisTitle);
  outDir->mkdir(xAxisTitle);
  outDir->cd(xAxisTitle);
  for (auto histName:histNames)
  {
    cout << histName << endl;
    auto writeName=histName;
    while (writeName.find("%")<writeName.size())
      writeName.erase(writeName.find("%"),writeName.find("i",writeName.find("%"))-writeName.find("%")+1);
    while (writeName.find("/")<writeName.size())
      writeName.replace(writeName.find("/"),1,"_");
    auto templateHist=(TH1*)files.at(0)->Get(Form(histName.c_str(), int(xAxis.at(0)), int(xAxis.at(0))));
    auto histTitle=templateHist->GetTitle();
    string yAxisTitle=templateHist->GetXaxis()->GetTitle();
    uint commaPosition=yAxisTitle.rfind(",");
    string units;
    if (commaPosition < yAxisTitle.size())
    {
      units=yAxisTitle.substr(commaPosition, yAxisTitle.size()-commaPosition);
      yAxisTitle.erase(commaPosition, yAxisTitle.size()-commaPosition);
    }
    vector <string> moments={"mean", "RMS"};
    for (int i=0;i<moments.size();i++)
    {
      vector <float> trend;
      auto moment=moments.at(i);
      for (int f=0;f<files.size();f++)
      {
        auto file=files.at(f);
        auto h=(TH1*)file->Get(Form(histName.c_str(), int(xAxis.at(f)), int(xAxis.at(f))));
        if (!h) 
        {
          printf("\nNot found in file %s\n", file->GetName());
          continue;
        }
        if (moment=="mean") 
          trend.push_back(h->GetMean()); 
        else if (moment=="RMS") 
          trend.push_back(h->GetRMS()); 
        else 
        {
          printf ("Undefined moment: %s", moment.c_str());
          assert(0);
        }
      }
      auto g=new TGraph(files.size(), &xAxis[0], &trend[0]);
      g->SetTitle(Form("%s (%s);%s;%s %s%s", histTitle, moment.c_str(), xAxisTitle, yAxisTitle.c_str(), moment.c_str(), units.c_str()));
      g->SetName(Form("%s_%s", writeName.c_str(), moment.c_str()));
      if (isTimeGraph) 
      {
        g->GetXaxis()->SetTimeDisplay(1);
        g->GetXaxis()->SetTimeFormat("%d/%m-%H:%M:%S");
        g->GetXaxis()->SetNdivisions(507);
      }
      g->Write();
    }
  }
}

vector <TFile*> makeFileVector (string &fileListPath)
{
  cout << "\nFile list:\n";
  vector <TFile*> files;
  if(fileListPath.rfind(".root") < fileListPath.size()) // root file
    files.push_back(new TFile(fileListPath.c_str()));
  else // file list
  {
    ifstream fileList;
    fileList.open(fileListPath);
    if(!fileList)
      printf("Error opening file: %s\n", fileListPath.c_str());
    string name;
    while(getline(fileList, name))
    {
      cout << name << endl;
      files.push_back(new TFile(name.c_str()));
    }
    fileList.close();
  }
  return files;
}

vector <float> getRunIds(vector<TFile*> &files, string &pattern)
{
  vector <float> runIds;
  for (auto f:files)
  { 
    smatch sm;
    string name=f->GetName();
    regex_search(name, sm, regex(pattern));
    string id=sm[1];
    runIds.push_back(stoi(id));
  }
  return runIds;
}

void filterObjectList (vector <string> &objectNames, string &pattern)
{
  if (pattern.size()==0)
    return;
  cout << "\nFiltering object list\n";
  regex re(pattern);
  vector <string> filteredNames;
  for (auto &name:objectNames)
  {
    smatch sm;
    if(regex_search(name, sm, regex(pattern)))
    {
      cout << name << endl;
      filteredNames.push_back(name);
    }
  }
  swap(filteredNames, objectNames);
  return;
}

vector <string> buildObjectList (vector <string> &objectNames, TDirectory *folder, int depth)
{
  int maxDepth=10000;
  string folderPath = folder -> GetPath();
  folderPath.erase (0,folderPath.rfind(":/")+2);
  while (folderPath.find("/")==0)
    folderPath.erase (0,1);
  TList *keys = folder -> GetListOfKeys ();
  if (depth==1)
    cout << "\nBuilding object list:\n";
  for (auto key : *keys)
  {
    string objectName = key -> GetName();
    if (objectName.find("run")==0)
      continue;
    auto object = dynamic_cast <TKey*> (key) -> ReadObj ();
    string className = object -> ClassName();
    if (!(className.find("3") < className.size()) && 
      (className.find("TH1")==0 || className.find("TProfile")==0 || className.find("Graph")==0))
    {
      string objectPath=Form("%s/%s", folderPath.c_str(), objectName.c_str());
      if(objectPath.find("/")==0)
        objectPath.erase(0,1);
      cout << objectPath << endl;
      objectNames.push_back(objectPath);
    }
    else if (className == "TDirectoryFile" && depth < maxDepth)
      buildObjectList (objectNames, dynamic_cast <TDirectoryFile*> (object), depth + 1);
    delete object;
  }
  if (depth==1)
    cout << objectNames.size() << " objects added.\n";
  return objectNames; 
}
