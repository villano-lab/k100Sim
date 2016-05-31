/*==================k100_AsciiOut.cc================================ 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Code support for the ascii output in the k100 simulation.
              
======================================================================*/

#include "G4RunManager.hh" 
#include "G4UnitsTable.hh"

#include "k100_AsciiOut.hh"


k100_AsciiOut::k100_AsciiOut()
{

}
k100_AsciiOut::k100_AsciiOut(G4String filename, char** vnames, G4int nvar, G4int dnum)
{
  dataFileOutputName = filename;
  ///////////////  HeaderString = header;
  NumVars = nvar;
  VariableNames = vnames;
  NumDet = dnum;
  //DetNames = dnames;
  
  dataFileOutput.open(dataFileOutputName);
  dataFileOutput.precision(10);
  WriteHeaderToFile();
}
k100_AsciiOut::~k100_AsciiOut()
{
  CloseDataFile();
}
void k100_AsciiOut::DumpToFile(G4double* dataArray , G4int nrows, G4int ncols)
{
  for (G4int jj=0; jj<nrows; jj++) {
    for (G4int ii=0; ii<ncols; ii++) {
      dataFileOutput << dataArray[jj*ncols + ii] << "\t";
      //G4cout << dataArray[jj*ncols + ii] << "          " << G4endl;
    }
    dataFileOutput << G4endl;
  }  
 
}
void k100_AsciiOut::WriteHeaderToFile()
{
  for (G4int ii=0; ii<NumDet; ii++) {
    //dataFileOutput << DetNames[ii] << "\t";
  }
  //dataFileOutput << G4endl;
  
  for (G4int ii=0; ii<NumVars; ii++) {
    dataFileOutput << VariableNames[ii] << "\t";
  }
  dataFileOutput << G4endl;  
}
void k100_AsciiOut::CloseDataFile()
{
  G4cout << "file name:  " << dataFileOutputName << G4endl;
  dataFileOutput.close();
}
 
// ------------------------------------------------

