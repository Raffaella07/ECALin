//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 15 11:16:29 2020 by ROOT version 6.12/07
// from TTree Events/Events
// found on file: /eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Oct14_withB/ParkingBPH1/crab_data_Run2018D_part1/201014_151829/0000/BParkNANO_data_2020Oct14_withB_23.root
//////////////////////////////////////////////////////////

#ifndef DileptonGeneralClass_h
#define DileptonGeneralClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class DileptonGeneralClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          nelectronPair;
   Float_t         electronPair_eta[69];   //[nelectronPair]
   Float_t         electronPair_l1_SC_Energy[69];   //[nelectronPair]
   Float_t         electronPair_l1_SC_RegEnergy[69];   //[nelectronPair]
   Float_t         electronPair_l1_SC_rawEnergy[69];   //[nelectronPair]
   Float_t         electronPair_l1_SCeta[69];   //[nelectronPair]
   Float_t         electronPair_l1_SCphi[69];   //[nelectronPair]
   Float_t         electronPair_l2_SC_Energy[69];   //[nelectronPair]
   Float_t         electronPair_l2_SC_RegEnergy[69];   //[nelectronPair]
   Float_t         electronPair_l2_SC_rawEnergy[69];   //[nelectronPair]
   Float_t         electronPair_l2_SCeta[69];   //[nelectronPair]
   Float_t         electronPair_l2_SCphi[69];   //[nelectronPair]
   Float_t         electronPair_mass[69];   //[nelectronPair]
   Float_t         electronPair_mll_llfit[69];   //[nelectronPair]
   Float_t         electronPair_mll_raw[69];   //[nelectronPair]
   Float_t         electronPair_phi[69];   //[nelectronPair]
   Float_t         electronPair_pt[69];   //[nelectronPair]
   Int_t           electronPair_charge[69];   //[nelectronPair]
   Int_t           electronPair_l1Idx[69];   //[nelectronPair]
   Int_t           electronPair_l2Idx[69];   //[nelectronPair]
   Int_t           electronPair_pdgId[69];   //[nelectronPair]
   UInt_t          nElectron;
   Float_t         Electron_dxy[6];   //[nElectron]
   Float_t         Electron_dxyErr[6];   //[nElectron]
   Float_t         Electron_dz[6];   //[nElectron]
   Float_t         Electron_dzErr[6];   //[nElectron]
   Float_t         Electron_energySC[6];   //[nElectron]
   Float_t         Electron_eta[6];   //[nElectron]
   Float_t         Electron_fBrem[6];   //[nElectron]
   Float_t         Electron_ip3d[6];   //[nElectron]
   Float_t         Electron_mass[6];   //[nElectron]
   Float_t         Electron_mvaId[6];   //[nElectron]
   Float_t         Electron_pfRelIso[6];   //[nElectron]
   Float_t         Electron_pfmvaId[6];   //[nElectron]
   Float_t         Electron_phi[6];   //[nElectron]
   Float_t         Electron_pt[6];   //[nElectron]
   Float_t         Electron_ptBiased[6];   //[nElectron]
   Float_t         Electron_sip3d[6];   //[nElectron]
   Float_t         Electron_trkRelIso[6];   //[nElectron]
   Float_t         Electron_unBiased[6];   //[nElectron]
   Float_t         Electron_vx[6];   //[nElectron]
   Float_t         Electron_vy[6];   //[nElectron]
   Float_t         Electron_vz[6];   //[nElectron]
   Int_t           Electron_charge[6];   //[nElectron]
   Int_t           Electron_pdgId[6];   //[nElectron]
   Bool_t          Electron_convVeto[6];   //[nElectron]
   Bool_t          Electron_isLowPt[6];   //[nElectron]
   Bool_t          Electron_isPF[6];   //[nElectron]
   Bool_t          Electron_isPFoverlap[6];   //[nElectron]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nelectronPair;   //!
   TBranch        *b_electronPair_eta;   //!
   TBranch        *b_electronPair_l1_SC_Energy;   //!
   TBranch        *b_electronPair_l1_SC_RegEnergy;   //!
   TBranch        *b_electronPair_l1_SC_rawEnergy;   //!
   TBranch        *b_electronPair_l1_SCeta;   //!
   TBranch        *b_electronPair_l1_SCphi;   //!
   TBranch        *b_electronPair_l2_SC_Energy;   //!
   TBranch        *b_electronPair_l2_SC_RegEnergy;   //!
   TBranch        *b_electronPair_l2_SC_rawEnergy;   //!
   TBranch        *b_electronPair_l2_SCeta;   //!
   TBranch        *b_electronPair_l2_SCphi;   //!
   TBranch        *b_electronPair_mass;   //!
   TBranch        *b_electronPair_mll_llfit;   //!
   TBranch        *b_electronPair_mll_raw;   //!
   TBranch        *b_electronPair_phi;   //!
   TBranch        *b_electronPair_pt;   //!
   TBranch        *b_electronPair_charge;   //!
   TBranch        *b_electronPair_l1Idx;   //!
   TBranch        *b_electronPair_l2Idx;   //!
   TBranch        *b_electronPair_pdgId;   //!
   TBranch        *b_nElectron;   //!
   TBranch        *b_Electron_dxy;   //!
   TBranch        *b_Electron_dxyErr;   //!
   TBranch        *b_Electron_dz;   //!
   TBranch        *b_Electron_dzErr;   //!
   TBranch        *b_Electron_energySC;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_fBrem;   //!
   TBranch        *b_Electron_ip3d;   //!
   TBranch        *b_Electron_mass;   //!
   TBranch        *b_Electron_mvaId;   //!
   TBranch        *b_Electron_pfRelIso;   //!
   TBranch        *b_Electron_pfmvaId;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_ptBiased;   //!
   TBranch        *b_Electron_sip3d;   //!
   TBranch        *b_Electron_trkRelIso;   //!
   TBranch        *b_Electron_unBiased;   //!
   TBranch        *b_Electron_vx;   //!
   TBranch        *b_Electron_vy;   //!
   TBranch        *b_Electron_vz;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_pdgId;   //!
   TBranch        *b_Electron_convVeto;   //!
   TBranch        *b_Electron_isLowPt;   //!
   TBranch        *b_Electron_isPF;   //!
   TBranch        *b_Electron_isPFoverlap;   //!

   virtual void     Init(TTree *tree);
};

#endif
