// 
void addMultiToyDataset(const char * workspaceName, const char * toyFileName, const char * toyName, const char* outputFileName, const int ntoys=1) {
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");

  // the workspace with the required POIs
  TFile workspace(workspaceName);
  RooWorkspace * myWS = workspace.Get("w");

  // the file with the toy 
  TFile fileWithToy(toyFileName);
//   RooAbsData * toyDataset = (RooAbsData *) fileWithToy.Get(Form("toys/%s",toyName));
  RooAbsData * toyDataset;
  
  for (int it = 1; it<=ntoys; it++){
    TString toy_tag;
//     TString toy_tag = (ntoys > 1 ) ? Form("toys/toy_%d",it) : TString(toyName);
    if (ntoys > 1 ) toy_tag = Form("toys/toy_%d",it);
    else toy_tag = TString(toyName);
    cout << toy_tag <<endl;
    
    toyDataset = (RooAbsData *) fileWithToy.Get(toy_tag);
    if (!toyDataset) {
    cout << "Error: toyDataset " << toy_tag
         << " not found in file " << toyFileName << endl;
    return;
  }
  toyDataset->SetName(toy_tag);
  // import (the two for completness
  myWS->import(*toyDataset);
  }
  
  // write on a new workspace file
  myWS->writeToFile(outputFileName);
}


