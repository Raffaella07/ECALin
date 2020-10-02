//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 15 17:19:59 2020 by ROOT version 6.12/07
// from TTree Events/Events
// found on file: crab_BuToKJpsi_ToEE__ext_0000.root
//////////////////////////////////////////////////////////

#ifndef DiLeptonMassClass_h
#define DiLeptonMassClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class DiLeptonMassClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   Float_t         genWeight;
   UInt_t          nPSWeight;
   Float_t         PSWeight[1];   //[nPSWeight]
   UChar_t         HLT_Mu7_IP4;
   UChar_t         HLT_Mu8_IP6;
   UChar_t         HLT_Mu8_IP5;
   UChar_t         HLT_Mu8_IP900;
   UChar_t         HLT_Mu8p5_IP3p5;
   UChar_t         HLT_Mu9_IP6;
   UChar_t         HLT_Mu9_IP5;
   UChar_t         HLT_Mu9_IP4;
   UChar_t         HLT_Mu10p5_IP3p5;
   UChar_t         HLT_Mu12_IP6;
   Int_t           PV_npvs;
   UInt_t          nDiLepton;
   Double_t        DiLepton_mll_fullfit[900];   //[nDiLepton]
   Double_t        DiLepton_l1Idx[900];   //[nDiLepton]
   Double_t        DiLepton_l2Idx[900];   //[nDiLepton]
   Double_t        DiLepton_l1isPF[900];   //[nDiLepton]
   Double_t        DiLepton_l1_pt[900];   //[nDiLepton]
   Double_t        DiLepton_l1_eta[900];   //[nDiLepton]
   Double_t        DiLepton_l1_phi[900];   //[nDiLepton]
   Double_t        DiLepton_l2isPF[900];   //[nDiLepton]
   Double_t        DiLepton_l2_pt[900];   //[nDiLepton]
   Double_t        DiLepton_l2_eta[900];   //[nDiLepton]
   Double_t        DiLepton_l2_phi[900];   //[nDiLepton]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_nPSWeight;   //!
   TBranch        *b_PSWeight;   //!
   TBranch        *b_HLT_Mu7_IP4;   //!
   TBranch        *b_HLT_Mu8_IP6;   //!
   TBranch        *b_HLT_Mu8_IP5;   //!
   TBranch        *b_HLT_Mu8_IP3;   //!
   TBranch        *b_HLT_Mu8p5_IP3p5;   //!
   TBranch        *b_HLT_Mu9_IP6;   //!
   TBranch        *b_HLT_Mu9_IP5;   //!
   TBranch        *b_HLT_Mu9_IP4;   //!
   TBranch        *b_HLT_Mu10p5_IP3p5;   //!
   TBranch        *b_HLT_Mu12_IP6;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_nDiLepton;   //!
   TBranch        *b_DiLepton_mll_fullfit;   //!
   TBranch        *b_DiLepton_l1Idx;   //!
   TBranch        *b_DiLepton_l2Idx;   //!
   TBranch        *b_DiLepton_l1isPF;   //!
   TBranch        *b_DiLepton_l1_pt;   //!
   TBranch        *b_DiLepton_l1_eta;   //!
   TBranch        *b_DiLepton_l1_phi;   //!
   TBranch        *b_DiLepton_l2isPF;   //!
   TBranch        *b_DiLepton_l2_pt;   //!
   TBranch        *b_DiLepton_l2_eta;   //!
   TBranch        *b_DiLepton_l2_phi;   //!

   virtual void     Init(TTree *tree);
};

#endif


