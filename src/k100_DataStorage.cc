/*==================k100_DataStorage.cc============================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12

      UPDATES:      
       

      PURPOSE:  Code support for the data storage class and routines
                to handle the organization output data into various
		forms which default to ascii supported by k100_ROOTOut.
              
======================================================================*/

#include "k100_DataStorage.hh"

// #ifdef ROOTOUT
// #include "k100_ROOTOut.hh"
// #endif
#include "k100_AsciiOut.hh"
#include "k100_ROOTOut.hh"
#include "k100_DetectorConstruction.hh"

k100_DataStorage::k100_DataStorage()
{
  printf(" k100_DATA storage -- EMPTY constructor\n");
}
k100_DataStorage::k100_DataStorage(G4String filename, G4int run, G4int rseed, G4bool rootFileBool)
{  

  NumDets = 1;

  //  storeEvtPrimaries = 0; //// This should come from RunAction (which is set via RunMessenger)

  n_data    = N_DATA;
  n_length  = N_LENGTH;
  n_entries = 0;
  n_files   = 0;

  runID = run;
  randSeed = rseed;
  outfilename = filename;

  rootOutFlag = rootFileBool;
  data = new G4double [N_DATA];
  dataArray = new G4double[N_LENGTH*N_DATA];
 
  
  VariableNames = new char*[N_DATA];
  char myNames[N_DATA][20];
  sprintf(myNames[0], "EV"); sprintf(myNames[1], "DT"); sprintf(myNames[2], "TS");
  sprintf(myNames[3], "P"); sprintf(myNames[4], "Type"); sprintf(myNames[5], "E1");
  sprintf(myNames[6], "D3");
  sprintf(myNames[7], "PX3"); sprintf(myNames[8], "PY3"); sprintf(myNames[9], "PZ3");
  sprintf(myNames[10], "X3"); sprintf(myNames[11], "Y3"); sprintf(myNames[12], "Z3");
  sprintf(myNames[13], "time3");
  sprintf(myNames[14], "PX1"); sprintf(myNames[15], "PY1"); sprintf(myNames[16], "PZ1");
  sprintf(myNames[17], "X1"); sprintf(myNames[18], "Y1"); sprintf(myNames[19], "Z1");
  sprintf(myNames[20], "time1");
  sprintf(myNames[21], "nCap"); 
  for (int kk=0; kk<N_DATA; kk++)
    VariableNames[kk] = myNames[kk];


// #ifdef ROOTOUT

//   // Open the data file for writing 
//   outfilename = filename + G4String("_") + NumtoStr(runID,3) + G4String("_"); 
//   Out = new k100_ROOTOut(outfilename + NumtoStr(n_files,3) + ".root", VariableNames, n_data, NumDets);

// #endif

  outfilename = filename + G4String("_") + NumtoStr(runID,3) + G4String("_"); 
  if(!rootOutFlag) {
    OutTEXT = new k100_AsciiOut(outfilename + NumtoStr(n_files,3) + ".txt", VariableNames, n_data, NumDets);
  }
  else {
    OutROOT = new k100_ROOTOut(outfilename + NumtoStr(n_files,3) + ".root", VariableNames, n_data, NumDets);
  }

}
k100_DataStorage::~k100_DataStorage()
{

// #ifdef ROOTOUT
//   Out->DumpToFile(dataArray , n_entries, n_data);
//   delete Out;
// #endif

  if(!rootOutFlag) {
    OutTEXT->DumpToFile(dataArray , n_entries, n_data);
    delete OutTEXT;
  }
  else {
    OutROOT->DumpToFile(dataArray , n_entries, n_data);
    delete OutROOT;
  }

  //delete data;  
  delete dataArray;
  G4cout << "Deleted dataArray." << G4endl;
}
void k100_DataStorage::setData(G4double* input)
{
  data = input;

}
G4double* k100_DataStorage::getData()
{
  return data;
}
void k100_DataStorage::addData()
{
 
  if (n_entries>n_length)
    writeArray();

  for (G4int ii=0; ii<n_data; ii++) {
    dataArray[n_entries*n_data + ii] = data[ii];
  }
  n_entries++;  



}
void k100_DataStorage::resetArray()
{

  for (G4int jj=0; jj<n_entries; jj++) {
    for (G4int ii=0; ii<n_data; ii++) {
      dataArray[jj*n_data + ii] = -99999;
      }
  }
  n_entries = 0;  
}
void k100_DataStorage::writeArray()
{

  char myNames[N_DATA][20];
  sprintf(myNames[0], "EV"); sprintf(myNames[1], "DT"); sprintf(myNames[2], "TS");
  sprintf(myNames[3], "P"); sprintf(myNames[4], "Type"); sprintf(myNames[5], "E1");
  sprintf(myNames[6], "D3");
  sprintf(myNames[7], "PX3"); sprintf(myNames[8], "PY3"); sprintf(myNames[9], "PZ3");
  sprintf(myNames[10], "X3"); sprintf(myNames[11], "Y3"); sprintf(myNames[12], "Z3");
  sprintf(myNames[13], "time3");
  sprintf(myNames[14], "PX1"); sprintf(myNames[15], "PY1"); sprintf(myNames[16], "PZ1");
  sprintf(myNames[17], "X1"); sprintf(myNames[18], "Y1"); sprintf(myNames[19], "Z1");
  sprintf(myNames[20], "time1");

  for (int kk=0; kk<N_DATA; kk++)
    VariableNames[kk] = myNames[kk];



// #ifdef ROOTOUT
//   Out->DumpToFile(dataArray , n_entries, n_data);
//   delete Out;
//   n_files++;

//   Out = new k100_ROOTOut(outfilename + NumtoStr(n_files,3) + G4String(".root"), VariableNames, n_data, NumDets);
// #endif

  if(!rootOutFlag) {
    OutTEXT->DumpToFile(dataArray , n_entries, n_data);
    delete OutTEXT;
    n_files++;

    OutTEXT = new k100_AsciiOut(outfilename + NumtoStr(n_files,3) + G4String(".txt"), VariableNames, n_data, NumDets);
  }
  else {
    OutROOT->DumpToFile(dataArray , n_entries, n_data);
    delete OutROOT;
    n_files++;

    OutROOT = new k100_ROOTOut(outfilename + NumtoStr(n_files,3) + G4String(".root"), VariableNames, n_data, NumDets);
  }

  resetArray();
}
void k100_DataStorage::printArray()
{

  for (G4int jj=0; jj<n_entries; jj++) {
    G4cout << "******** DataStorage : -----" << G4endl;
    for (G4int ii=0; ii<n_data; ii++) {
      G4cout << "data["<< jj << "," << ii <<  "] = " << dataArray[jj*n_data + ii] << " - ";
      }
    G4cout << G4endl;
  }
  G4cout << G4endl;

}
G4String k100_DataStorage::NumtoStr(G4int number, G4int length)
{
  //G4String name;
  std::ostringstream oss;

  oss << std::setfill('0');
  oss << std::setw(length) << number;

  return oss.str();
}
G4bool k100_DataStorage::overflowArray(G4int hits)
{
  if (n_entries+hits > n_length) return TRUE;
  else return FALSE;
}

