/*==================k100_ROOTOut.cc================================ 
   
      PROGRAMMER:  Shubham Pandey 4 March 2022   

      UPDATES:      
       

      PURPOSE: Code support for the ROOT output in the k100 simulation.
              
======================================================================*/

#include "G4RunManager.hh" 
#include "G4UnitsTable.hh"

#include "k100_ROOTOut.hh"


k100_ROOTOut::k100_ROOTOut()
{

}
k100_ROOTOut::k100_ROOTOut(G4String filename, char** vnames, G4int nvar, G4int dnum)
{
  dataFileOutputName = filename;
  ///////////////  HeaderString = header;
  NumVars = nvar;
  VariableNames = vnames;
  NumDet = dnum;
  //DetNames = dnames;
  dataFileOutput = TFile::Open(dataFileOutputName,"RECREATE");
  simtree = new TTree("simtree","simulation output");

  simtree->Branch("EV",&EV,"EV/L");
  simtree->Branch("DT",&DT);
  simtree->Branch("TS",&TS);
  simtree->Branch("P",&P);
  simtree->Branch("Type",&Type);
  simtree->Branch("E1",&E1);
  simtree->Branch("D3",&D3);
  simtree->Branch("X1",&X1);
  simtree->Branch("Y1",&Y1);
  simtree->Branch("Z1",&Z1);
  simtree->Branch("X3",&X3);
  simtree->Branch("Y3",&Y3);
  simtree->Branch("Z3",&Z3);
  simtree->Branch("PX1",&PX1);
  simtree->Branch("PY1",&PY1);
  simtree->Branch("PZ1",&PZ1);
  simtree->Branch("PX3",&PX3);
  simtree->Branch("PY3",&PY3);
  simtree->Branch("PZ3",&PZ3);
  simtree->Branch("time1",&time1);
  simtree->Branch("time3",&time3);
  simtree->Branch("nCap",&nCap);

  // dataFileOutput.open(dataFileOutputName);
  // dataFileOutput.precision(10);
  // WriteHeaderToFile();
}
k100_ROOTOut::~k100_ROOTOut()
{
  CloseDataFile();
}
void k100_ROOTOut::DumpToFile(G4double* dataArray , G4int nrows, G4int ncols) {
  Long64_t evnow;
  for (G4int jj=0; jj<nrows; jj++) {
    
    if(jj == 0) evnow = dataArray[0];
    if(evnow != dataArray[jj*ncols]) {
      EV = evnow;
      simtree->Fill();
      clearVariables();
      evnow = dataArray[jj*ncols];
    } 

    DT.push_back(dataArray[jj*ncols + 1]);
    TS.push_back(dataArray[jj*ncols + 2]);
    P.push_back(dataArray[jj*ncols + 3]);
    Type.push_back(dataArray[jj*ncols + 4]);
    E1.push_back(dataArray[jj*ncols + 5]);
    D3.push_back(dataArray[jj*ncols + 6]);
    PX3.push_back(dataArray[jj*ncols + 7]);
    PY3.push_back(dataArray[jj*ncols + 8]);
    PZ3.push_back(dataArray[jj*ncols + 9]);
    X3.push_back(dataArray[jj*ncols + 10]);
    Y3.push_back(dataArray[jj*ncols + 11]);
    Z3.push_back(dataArray[jj*ncols + 12]);
    time1.push_back(dataArray[jj*ncols + 13]);
    PX1.push_back(dataArray[jj*ncols + 14]);
    PY1.push_back(dataArray[jj*ncols + 15]);
    PZ1.push_back(dataArray[jj*ncols + 16]);
    X1.push_back(dataArray[jj*ncols + 17]);
    Y1.push_back(dataArray[jj*ncols + 18]);
    Z1.push_back(dataArray[jj*ncols + 19]);
    time3.push_back(dataArray[jj*ncols + 20]);
    nCap.push_back(dataArray[jj*ncols + 21]);
  


    // for (G4int ii=0; ii<ncols; ii++) {
    //   dataFileOutput << dataArray[jj*ncols + ii] << "\t";
    //   //G4cout << dataArray[jj*ncols + ii] << "          " << G4endl;
    // }
    // dataFileOutput << G4endl;
  }  
  EV = evnow;
  simtree->Fill();
 
}
// void k100_ROOTOut::WriteHeaderToFile()
// {
//   for (G4int ii=0; ii<NumDet; ii++) {
//     //dataFileOutput << DetNames[ii] << "\t";
//   }
//   //dataFileOutput << G4endl;
  
//   for (G4int ii=0; ii<NumVars; ii++) {
//     dataFileOutput << VariableNames[ii] << "\t";
//   }
//   dataFileOutput << G4endl;  
// }

void k100_ROOTOut::clearVariables() {
  DT.clear();
  TS.clear();
  P.clear();
  Type.clear();
  E1.clear();
  D3.clear();
  X1.clear();
  Y1.clear();
  Z1.clear();
  X3.clear();
  Y3.clear();
  Z3.clear();
  PX1.clear();
  PY1.clear();
  PZ1.clear();
  PX3.clear();
  PY3.clear();
  PZ3.clear();
  time1.clear();
  time3.clear();
  nCap.clear();
}

void k100_ROOTOut::CloseDataFile()
{
  G4cout << "file name:  " << dataFileOutputName << G4endl;
  dataFileOutput->cd();
  simtree->Write();
  dataFileOutput->Write();
  dataFileOutput->Close();
}
 
// ------------------------------------------------

