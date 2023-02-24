/*==================k100Sim.cc========================================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE:  Create a simulation package for the k100 fridge including
                all external shielding.
              
      INPUT:      

      OUTPUT:   
              
======================================================================*/
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <sys/time.h>
#include <time.h>

//command line options
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include "Randomize.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "G4PhysListFactory.hh"

#include "QGSP_BERT_HP.hh"
#include "k100_DetectorConstruction.hh"
#include "k100_PrimaryGeneratorAction.hh"
#include "k100_RunAction.hh"
#include "k100_EventAction.hh"
#include "k100_SteppingAction.hh"
#include "Shielding_ComptonsUpdate.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//The name of this program. 
const char* program_name;

//Prints usage information for this program to STREAM (typically
//stdout or stderr), and exit the program with EXIT_CODE.  Does not
//return. 

void print_usage (FILE* stream, int exit_code)
{
  fprintf (stream, "Usage:  %s options [ inputfile(s) ]\n", program_name);
  fprintf (stream,
	   //"\n"
           "  -c, --customgen                    use custom particle generator for neutron captures\n"
           "  -i, --infile       <filename>      name the input nrCasacadeSim file \n"
           "  -o, --outfile       <filename>     name the output file \n"
           "  -p, --only-ncap                    restrict event output to those including ncap \n"
           "  -r, --rootoutput                   generate output in root tree format \n"
           "  -q, --quiet         <level>        quiet printing \n"
           "                                     optional argument integer > 0 \n"
           "                                     0,default: no non-event output\n"
           "                                     1,: no standard event output\n"
           "  -s, --silent                       silent, no standard out \n"
           "  -v, --verbose       <level>        Print verbose messages at level <level>\n"
           "  -V, --version                      print version and exit\n");

  exit (exit_code);
}

#define no_argument 0
#define required_argument 1 
#define optional_argument 2

//using namespace std;

int main(int argc, char** argv) {

  //***********Begin get input*********************//
 
  //set parameters for input system
  G4String inputfile;
  std::string outputfile;
  uint verbosity=0;
  uint quietness=0;
  bool quiet=false;
  bool dataquiet=false;
  bool customgen=false;
  bool onlyncap=false;
  bool rootOutput=false;
   
  const struct option longopts[] =
  {
    {"customgen",     no_argument,  0, 'c'},
    {"infile",     required_argument,  0, 'i'},
    {"rootOutput",     no_argument,  0, 'r'},
    {"outfile",     required_argument,  0, 'o'},
    {"only-ncap",     no_argument,  0, 'p'},
    {"quiet",     optional_argument,  0, 'q'},
    {"silent",    no_argument,        0, 's'},
    {"verbose",   optional_argument,  0, 'v'},
    {"version",   no_argument,        0, 'V'},
    {0,0,0,0},
  };

  int index;
  int iarg=0;

  //turn off getopt error message
  opterr=1; 
  bool infile = false;
  while(iarg != -1)
  {
    iarg = getopt_long(argc, argv, "+co:ri:pq::sv::V", longopts, &index);

    switch (iarg)
    {

      case 'c':
        customgen = true;
        break;

      case 'i':
        inputfile = optarg;
        infile = true;
        break;
      
      case 'r':
        rootOutput = true;
        break;


      case 'o':
        outputfile = optarg;
        break;

      case 'p':
        onlyncap=true;
        break;

      case 'q':
        if(optarg)
          quietness = atoi(optarg);

        if(quietness == 0)
          quiet = true;
        else if(quietness == 1)
          dataquiet = true;
        else if(quietness == 2){
          dataquiet = true;
          quiet = true;
        }
        else{
       	  std::cerr << "ERROR eventAlignCheck: invalid quietness value" << std::endl;
          print_usage(stderr,1);
        }
        break;

      case 's':
        quiet = true;
        dataquiet = true;
        break;

      case 'v':
        if(optarg)
          verbosity = atoi(optarg);
        else
          verbosity = 1;
          break;

      case 'V':
        printf("Version: %s\n", __GIT_VERSION);
        return 0;
        break;

      //case for long only options
      case 0:
        break;

      case '?':
        print_usage(stderr,1);
        break;
    }
  } 

  //get the rest of the filenames in a vector of strings
  if(optind==argc){
    std::cerr << "eventAlignCheck: ERROR! no file supplied" << std::endl;
    exit(1);
  }
  if(customgen && !infile) {
    std::cout<<"Please provide input file if customgen flag is ON. Exiting..."<<std::endl;
    exit(0);
  }
  std::vector<std::string> filenames;
  for(int i = optind; i < argc; i++){
    if(verbosity>=1){
      //print the filenames on a line
      printf("%s\n",argv[i]);
    }
    filenames.push_back(argv[i]);
  }

  // std::cout<<"rootFileoutput = "<<rootOutput<<std::endl;
  // exit(0);
  //***********End get input*********************//

  for(int i=0;i<filenames.size();i++){
    std::cout << filenames[i] << std::endl;
  }

#ifdef G4VIS_USE
  // Visualisation, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  //Run Manager
  G4RunManager * runManager = new G4RunManager;

  //make a physics list factory
  G4PhysListFactory *g4Factory = new G4PhysListFactory();
  char* physics_list = new char[50];
  //sprintf(physics_list,"QGSP_BERT_HP");
  sprintf(physics_list,"Shielding");
  // User Initializaton classes (mandatory)
  runManager->SetUserInitialization(new k100_DetectorConstruction());
  //runManager->SetUserInitialization(new Shielding_ComptonsUpdate);
  //try standard shielding: "Shielding_EMZ"
  //runManager->SetUserInitialization(g4Factory->GetReferencePhysList("Shielding_EMZ"));
  //revert back to standard to try to fix infinite loop bug: N-MISC-17-003 pg 10
  //runManager->SetUserInitialization(g4Factory->GetReferencePhysList("Shielding"));
  
  runManager->SetUserInitialization(g4Factory->GetReferencePhysList(physics_list));
  G4cout<<"Using physics list = "<<physics_list<<G4endl;

  // UserAction Classes============
  // event generator
  
  k100_PrimaryGeneratorAction* myPrimaryEventGenerator=new k100_PrimaryGeneratorAction(customgen,inputfile); //sourceGun is the custom generator  
  runManager->SetUserAction(myPrimaryEventGenerator);

  //run action
  k100_RunAction* pRunAction = new k100_RunAction(rootOutput);
  runManager->SetUserAction(pRunAction);
  k100_EventAction* pEventAction = new k100_EventAction(pRunAction,onlyncap);
  runManager->SetUserAction(pEventAction);
  //runManager->SetUserAction(new k100_EventAction(pRunAction,onlyncap));
  runManager->SetUserAction(new k100_SteppingAction(pEventAction));

  // Initialize G4 kernel
  runManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager * UI = G4UImanager::GetUIpointer();

  // open G4UIterminal interactive terminal 
  G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
  session = new G4UIterminal(new G4UItcsh);
  G4cout << "\n---> I am in a smart tcsh terminal ;-) ." << G4endl;
#else
  session = new G4UIterminal();
  G4cout << "\n---> I am in a dumb terminal ;-( ." << G4endl;
#endif

  //start session
  for(int i=0;i<filenames.size();i++){
    G4String command = "/control/execute ";
    G4String macroFileName = filenames[i];
    G4cout << " Macro " << macroFileName << " is being run " << G4endl;
    UI->ApplyCommand(command+macroFileName);
  }
  //session->SessionStart();

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}
