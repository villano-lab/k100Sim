/*==================k100_nrCascadeTree.hh==================================== 
   ƒ
      PROGRAMMER:  Shubham Pandey 7 March 2022  

      UPDATES:      
       

      PURPOSE: Class for reading NR cascade ROOT file to be used as new particle source in
               k100 simulation.
              
======================================================================*/

#ifndef k100_NRCASCADETREE_h
#define k100_NRCASCADETREE_h 1

// ------------------------------------------------

#include "globals.hh"
#include <fstream>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
//#include "vector"
using namespace std;
// ------------------------------------------------

class k100_nrCascadeTree
{
public :
  k100_nrCascadeTree();
  k100_nrCascadeTree(G4String);
  ~k100_nrCascadeTree();
  //bool Init(TTree *inTree);


  
  //void WriteHeaderToFile();
  
  
  TFile* inFile;
  TTree* tree;

  Long64_t        Event;
  Short_t         nCap_flag;
  //Short_t         n_gamma;
  Int_t           nE_gamma;
  Double_t        E_gamma[4];   //[nE_gamma]
  Double_t        x_cap;
  Double_t        y_cap;
  Double_t        z_cap;

  // List of branches
  TBranch        *b_Event;   //!
  TBranch        *b_nCap_flag;   //!
  //TBranch        *b_n_gamma;   //!
  TBranch        *b_nE_gamma;   //!
  TBranch        *b_E_gamma;   //!
  TBranch        *b_x_cap;   //!
  TBranch        *b_y_cap;   //!
  TBranch        *b_z_cap;   //!

};

// bool k100_nrCascadeTree::Init(TTree *inTree)
// {
//    // The Init() function is called when the selector needs to initialize
//    // a new tree or chain. Typically here the branch addresses and branch
//    // pointers of the tree will be set.
//    // It is normally not necessary to make changes to the generated
//    // code, but the routine can be extended by the user if needed.
//    // Init() will be called many times when running on PROOF
//    // (once per file to be processed).

//    // Set branch addresses and branch pointers
   
//    inTree->SetMakeClass(1);

//    inTree->SetBranchAddress("Event", &Event, &b_Event);
//    inTree->SetBranchAddress("nCap_flag", &nCap_flag, &b_nCap_flag);
//    inTree->SetBranchAddress("n_gamma", &n_gamma, &b_n_gamma);
//    inTree->SetBranchAddress("nE_gamma", &nE_gamma, &b_nE_gamma);
//    inTree->SetBranchAddress("E_gamma", E_gamma, &b_E_gamma);
//    inTree->SetBranchAddress("x_cap", &x_cap, &b_x_cap);
//    inTree->SetBranchAddress("y_cap", &y_cap, &b_y_cap);
//    inTree->SetBranchAddress("z_cap", &z_cap, &b_z_cap);
//    return true;
// }
#endif
