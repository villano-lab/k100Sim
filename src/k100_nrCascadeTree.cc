/*==================k100_nrCascadeTree.cc================================ 
   
      PROGRAMMER:  Shubham Pandey 4 March 2022   

      UPDATES:      
       

      PURPOSE: Code support for the ROOT output in the k100 simulation.
              
======================================================================*/

#include "G4RunManager.hh" 
#include "G4UnitsTable.hh"
#include <iostream>
#include "k100_nrCascadeTree.hh"


k100_nrCascadeTree::k100_nrCascadeTree()
{

}
k100_nrCascadeTree::k100_nrCascadeTree(G4String filename)
{
  inFile = TFile::Open(filename);
  if(!inFile) {
    G4cout<<"Could not find file: "<<filename<<G4endl;
    exit(0);
  }
  tree = (TTree*)inFile->Get("combined");
  if(!tree) {
    G4cout<<"Could not open \'combined\' tree in file: "<<inFile<<G4endl;
    exit(0);
  }

  tree->SetMakeClass(1);

  tree->SetBranchAddress("Event", &Event, &b_Event);
  tree->SetBranchAddress("nCap_flag", &nCap_flag, &b_nCap_flag);
  //tree->SetBranchAddress("n_gamma", &n_gamma, &b_n_gamma);
  tree->SetBranchAddress("nE_gamma", &nE_gamma, &b_nE_gamma);
  tree->SetBranchAddress("E_gamma", E_gamma, &b_E_gamma);
  tree->SetBranchAddress("x_cap", &x_cap, &b_x_cap);
  tree->SetBranchAddress("y_cap", &y_cap, &b_y_cap);
  tree->SetBranchAddress("z_cap", &z_cap, &b_z_cap);

  // if(Init(tree)) {
  //   G4cout<<"Tree initialized successfully."<<G4endl;
  // }
  
}
k100_nrCascadeTree::~k100_nrCascadeTree()
{
  
}

 
// ------------------------------------------------

