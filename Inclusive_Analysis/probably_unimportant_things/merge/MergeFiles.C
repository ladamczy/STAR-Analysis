
//_____________________________________________________________________________
void MergeFiles(const string tmpnam, const string outfile) {

  //merge files in 'tmpnam' into 'outfile'

  //load to prevent error messages about missing dictionary class
  gSystem->Load("../build/libstar-upc.so");

  //load files to the merger
  TFileMerger *merg = new TFileMerger();
  cout << "Loading files for merging" << endl;
  ifstream in1(tmpnam.c_str());
  string line;
  Int_t nFiles = 0;
  while(getline(in1, line)) {
    merg->AddFile(line.c_str(), kFALSE);
    nFiles++;
  }
  //files loaded
  in1.close();
  if( nFiles == 0 ) {cout << "No files to merge, exiting." << endl; return;}
  cout << nFiles << " files loaded" << endl;

  //perform the merging
  merg->OutputFile( outfile.c_str() );
  Bool_t stat = merg->Merge();
  if (stat) {
    cout << "Successfully merged " << nFiles << " files" << endl;
  }
  else {
    cout << "Error in merging" << endl;
  }
  delete merg;

}//MergeFiles





















