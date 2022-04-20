/*=======================k100_DataStorage.hh============================ 
   
      PROGRAMMER:  Anthony Villano 10/18/12

      UPDATES:      
       

      PURPOSE: Class to specify the format for data storage and handle
               calls to the k100_ROOTOut class, which writes the
	       data.
              
======================================================================*/

#ifndef k100_DataStorage_H
#define k100_DataStorage_H 1

//#define ROOTOUT 1


#include "globals.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"

#include <list>

#define N_DATA 22 
#define N_LENGTH 200000

// #ifdef ROOTOUT
// class k100_ROOTOut;
// #endif

class k100_ROOTOut;
class k100_AsciiOut;
class k100_DataStorage
{
public:
  k100_DataStorage();
  k100_DataStorage(G4String, G4int, G4int, G4bool);
	
  ~k100_DataStorage();
	
private:
  G4int storeEvtPrimaries,lastGE;
  G4int n_data,n_length; 
  G4int n_adds,n_entries,n_files;
  
  G4int runID, randSeed;
  G4String outfilename;
  

  G4double* data;
  G4double* dataArray;

  char** DetNames;
  char** VariableNames;
  G4int NumDets;

  G4bool rootOutFlag;
  k100_ROOTOut* OutROOT;
  k100_AsciiOut* OutTEXT;
// #ifdef ROOTOUT
//   k100_ROOTOut* Out;
// #endif

public:		
  // member functions
  void setData(double*);
  G4double* getData();
  void addData();
  void printArray();
  void writeArray();
  G4bool overflowArray(G4int);

private:
  void resetArray();
  G4String NumtoStr(G4int, G4int);

public:

};

#endif
