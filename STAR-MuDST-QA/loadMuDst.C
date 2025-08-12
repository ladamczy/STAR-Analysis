cd ~/work/pp500_2017/ExclusiveAnalysis/

# Create the loadMuDst.C file
cat > loadMuDst.C << 'EOF'
void loadMuDst() {
    // For grid jobs, most libraries should already be loaded by root4star
    // But we can ensure critical ones are available
    
    cout << "Loading MuDst libraries for grid job..." << endl;
    
    // Load basic STAR libraries if not already loaded
    if (!gSystem->GetLibraries("St_base")) {
        gSystem->Load("St_base");
    }
    if (!gSystem->GetLibraries("StChain")) {
        gSystem->Load("StChain");
    }
    if (!gSystem->GetLibraries("StMuDSTMaker")) {
        gSystem->Load("StMuDSTMaker");
    }
    
    cout << "MuDst libraries loaded successfully." << endl;
}
EOF
