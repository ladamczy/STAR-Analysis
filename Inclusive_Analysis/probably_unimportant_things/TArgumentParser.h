
//_____________________________________________________________________________
class TArgumentParser: public TMap {

public:

  TArgumentParser() { SetOwner(); }

  void AddBool(const string &name, Bool_t *par) { Add(new TObjString(name.c_str()), new w_boolp(par)); }
  void AddInt(const string &name, Int_t *par) { Add(new TObjString(name.c_str()), new w_intp(par)); }
  void AddFloat(const string &name, Float_t *par) { Add(new TObjString(name.c_str()), new w_floatp(par)); }
  void AddFuncInt3(const string &name, void (*par)(UInt_t, Int_t, Int_t)) { Add(new TObjString(name.c_str()), new w_int3p(par)); }

  //_____________________________________________________________________________
  void Parse(const string &config) {

    //open config file
    ifstream inp(config.c_str());
    string line;

    //config file loop
    while( getline(inp, line) ) {

      if( line.empty() == true ) continue;

      //remove white space characters at the beginning
      line = line.substr(line.find_first_not_of(" \t"));
      //skip comments
      if( line.find("#") == 0 ) continue;

      //split the line to tokens
      TString strlin(line.c_str());
      TObjArray *arr = strlin.Tokenize(" \t");
      TIterator *tok_iter = arr->MakeIterator();
 
      //get the parameter name
      TString name = StringFromIter(tok_iter);

      //load value for the parameter
      LoadPar<w_boolp>(name, tok_iter, "bool");
      LoadPar<w_intp>(name, tok_iter, "int");
      LoadPar<w_floatp>(name, tok_iter, "float");
      LoadInt3(name, tok_iter);

      delete tok_iter;
      delete arr;

    }//config file loop

  }//Parse

private:

  //_____________________________________________________________________________
  template<class par_type> void LoadPar(const TString &name, TIterator *tok_iter, string type) {

    //load parameter of a given par_type

    par_type *valp = dynamic_cast<par_type*>( GetValue(name) );
    if( !valp ) return;

    //explicit type checking necessary in ROOT 5 interpreter
    if( valp->type != type ) return;

    //parameter found, get the value
    istringstream st(StringFromIter(tok_iter).Data());
    st >> *(valp->val);

  }//LoadPar

  //_____________________________________________________________________________
  void LoadInt3(const TString &name, TIterator *tok_iter) {

    //load pointer to void (UInt_t, Int_t, Int_t) function

    w_int3p *valp = dynamic_cast<w_int3p*>( GetValue(name) );
    if( !valp ) return;

    //explicit type checking necessary in ROOT 5 interpreter
    if( valp->type != "int3" ) return;

    //function pointer found, get arguments
    istringstream st(StringFromIter(tok_iter).Data());
    UInt_t arg0;
    st >> arg0;
    Int_t arg12[2] = {0, 0}; // default run range
    for(Int_t i=0; i<2; i++) {
      string inp = StringFromIter(tok_iter).Data();
      if( inp.empty() == true || inp.find("#") == 0 ) break;
      istringstream st(inp);
      st >> arg12[i];
    }

    //call the function
    (*(valp->val))(arg0, arg12[0], arg12[1]);

  }//LoadInt3

  //wrapper for bool pointer
  struct w_boolp : public TObject {
    Bool_t *val;
    string type;
    w_boolp(Bool_t *bp): val(bp), type("bool") {}
  };//w_boolp

  //wrapper for integer pointer
  struct w_intp : public TObject {
    Int_t *val;
    string type;
    w_intp(Int_t *ip): val(ip), type("int") {}
  };//w_intp

  //float pointer
  struct w_floatp : public TObject {
    Float_t *val;
    string type;
    w_floatp(Float_t *fp): val(fp), type("float") {}
  };//w_floatp

  //pointer to void (UInt_t, Int_t, Int_t) function
  struct w_int3p : public TObject {
    void (*val)(UInt_t, Int_t, Int_t);
    string type;
    w_int3p(void (*i3p)(UInt_t, Int_t, Int_t)): type("int3") {
      val = i3p;
    }
  };//w_intp

  //_____________________________________________________________________________
  TString StringFromIter(TIterator *tok_iter) {

    //retrieve string from token iterator

    TObject *obj = tok_iter->Next();
    if( !obj ) return TString("");

    return (dynamic_cast<TObjString*>(obj))->GetString();

  }//StringFromIter

};//ArgumentParser




















