/*==================k100_AsciiOut.hh==================================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Class for writing stored data to ascii file in the
               k100 simulation.
              
======================================================================*/

#ifndef k100_AsciiOut_h
#define k100_AsciiOut_h 1

// ------------------------------------------------

#include "globals.hh"
#include <fstream>


// ------------------------------------------------

class k100_AsciiOut
{
public :
  k100_AsciiOut();
  k100_AsciiOut(G4String, char**, G4int, G4int);
  ~k100_AsciiOut();

public :
  
  void WriteHeaderToFile();
  void CloseDataFile();
  void DumpToFile(G4double* , G4int, G4int);
  
private :

  G4String dataFileOutputName;
  G4int NumDet;
  G4int NumVars;
  char** DetNames;
  char** VariableNames;

  std::ofstream dataFileOutput;
};

#endif
