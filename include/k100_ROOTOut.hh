/*==================k100_ROOTOut.hh==================================== 
   ƒ
      PROGRAMMER:  Shubham Pandey 4 March 2022  

      UPDATES:      
       

      PURPOSE: Class for writing stored data to ROOT file in the
               k100 simulation.
              
======================================================================*/

#ifndef k100_ROOTOut_h
#define k100_ROOTOut_h 1

// ------------------------------------------------

#include "globals.hh"
#include <fstream>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
//#include "vector"
using namespace std;
// ------------------------------------------------

class k100_ROOTOut
{
public :
  k100_ROOTOut();
  k100_ROOTOut(G4String, char**, G4int, G4int);
  ~k100_ROOTOut();

public :
  
  //void WriteHeaderToFile();
  void CloseDataFile();
  void clearVariables();
  void DumpToFile(G4double* , G4int, G4int);
  
private :

  G4String dataFileOutputName;
  G4int NumDet;
  G4int NumVars;
  char** DetNames;
  char** VariableNames;

  TFile* dataFileOutput;
  TTree* simtree;
  

  Long64_t EV;
  vector<Int_t> DT;
  vector<long> TS;
  vector<long> P;
  vector<Int_t> Type;
  vector<float> E1;
  vector<float> D3;
  vector<float> X1;
  vector<float> Y1;
  vector<float> Z1;
  vector<float> X3;
  vector<float> Y3;
  vector<float> Z3;
  vector<float> PX1;
  vector<float> PY1;
  vector<float> PZ1;
  vector<float> PX3;
  vector<float> PY3;
  vector<float> PZ3;
  vector<float> time1;
  vector<float> time3;
  vector<Int_t> nCap;
};

#endif
