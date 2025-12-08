std::map<int, int> date;
std::map<int, int> time;

// ----------------------------------------------------------------------------
bool loadTimestampFile(){
    ifstream fp;
    fp.open("./StRoot/run17_date_time.txt");
    if(!fp.is_open()){
        printf("Run timestamp file could not be loaded\n");
        return false;
    }
    int runNumb;
    int DateSt;
    int timeSt;
    int rowNumber = 0;
    while(fp.good()){
        std::string line;
        fp>>DateSt>>timeSt>>runNumb;
        date[runNumb] = DateSt;
        time[runNumb] = timeSt;
        rowNumber++;
    }
    printf("Imported %d rows\n", rowNumber);

    return true;
}
// ----------------------------------------------------------------------------
bool getRunStart(int runNumber, int* date_time_array){
    loadTimestampFile();
    //checking if run actually exists on a list
    auto found_date = date.find(runNumber);
    auto found_time = time.find(runNumber);
    if(found_date==date.end()||found_time==time.end()){
        printf("Run does not exist in the list, possibly it isn't from Run17?\n");
        return false;
    }
    date_time_array[0] = date[runNumber];
    date_time_array[1] = time[runNumber];
    printf("Date set:\t%d\nTime set:\t%d\n", date_time_array[0], date_time_array[1]);
}
// ----------------------------------------------------------------------------
std::string to_string(int input){
    TString inside;
    inside.Form("%d", input);
    std::string str(inside.Data());
    if(6>str.size())
        str.insert(0, 6-str.size(), '0');
    return str;
}
// ----------------------------------------------------------------------------
void bfc_run(std::string fileToProcess, int runNumber = 18091010){
    int date_time_array[2] = { 0,0 };
    if(runNumber==0)
        runNumber = 18091010;
    getRunStart(runNumber, date_time_array);
    std::string first_part = ".x bfc.C(1005,\"DbV20200225,y2017a,agml,StiCA,btof,mtd,mtdCalib,pp2pp,PicoVtxDefault,UseBTOFmatchOnly,VFStoreX,fmsDat,fmsPoint,fpsDat,OSpaceZ2,OGridLeak3D,fzin,usexgeom,l0,FieldOn,MakeEvent,Idst,BAna,Tree,Tree,logger,genvtx,TpcRS,TpxClu,tpcDb,bbcSim,btofsim,btofMixer,btofMatch,emcy2,eefs,geantout,-dstout,IdTruth,big,clearmem,CorrX,gen_T,geomT,sim_T,evout,TpcHitMover,tags,MiniMcMk,sdt";
    first_part = first_part+to_string(date_time_array[0])+"."+to_string(date_time_array[1]);
    first_part = first_part+",emcDY2,BtofDat,picoWrite\",\""+fileToProcess+"\")";
    printf("Processed line:\n%s\n", first_part.c_str());
    gROOT->ProcessLine(first_part.c_str());
}