#ifndef MeeClass_h
#define MeeClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class MeeClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   

   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   Double_t	   mll_EE;
   Double_t        mll_ET;
   Double_t        mll_TT;
   Double_t        B_mass;
   Double_t        PUweight;
   Double_t        deltaPhi;
   Double_t	   deltaEta;
   Double_t        deltaR;

	

   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_mll_EE;    //!
   TBranch        *b_mll_ET;    //!
   TBranch        *b_mll_TT;    //!
   TBranch        *b_B_mass;    //!
   TBranch        *b_PUweight;    //!
   TBranch        *b_deltaPhi; //!
   TBranch        *b_deltaEta;
   TBranch        *b_deltaR;
// Fixed size dimensions of array or collections stored in the TTree if any.
   virtual void     Init();
   virtual void     InitTree(TTree *tree);
};

#endif
