//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct  9 11:19:49 2020 by ROOT version 6.22/02
// from TTree Events/Events
// found on file: /eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Oct08_withB/EGamma/crab_data_EGamma/201008_171241/0000/BParkNANO_data_2020Oct08_withB_544.root
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
   UInt_t          nElectron;
   Float_t         Electron_dxy[7];   //[nElectron]
   Float_t         Electron_dxyErr[7];   //[nElectron]
   Float_t         Electron_dz[7];   //[nElectron]
   Float_t         Electron_dzErr[7];   //[nElectron]
   Float_t         Electron_energySC[7];   //[nElectron]
   Float_t         Electron_eta[7];   //[nElectron]
   Float_t         Electron_fBrem[7];   //[nElectron]
   Float_t         Electron_ip3d[7];   //[nElectron]
   Float_t         Electron_mass[7];   //[nElectron]
   Float_t         Electron_mvaId[7];   //[nElectron]
   Float_t         Electron_pfRelIso[7];   //[nElectron]
   Float_t         Electron_pfmvaId[7];   //[nElectron]
   Float_t         Electron_phi[7];   //[nElectron]
   Float_t         Electron_pt[7];   //[nElectron]
   Float_t         Electron_ptBiased[7];   //[nElectron]
   Float_t         Electron_sip3d[7];   //[nElectron]
   Float_t         Electron_trkRelIso[7];   //[nElectron]
   Float_t         Electron_unBiased[7];   //[nElectron]
   Float_t         Electron_vx[7];   //[nElectron]
   Float_t         Electron_vy[7];   //[nElectron]
   Float_t         Electron_vz[7];   //[nElectron]
   Int_t           Electron_charge[7];   //[nElectron]
   Int_t           Electron_pdgId[7];   //[nElectron]
   Bool_t          Electron_convVeto[7];   //[nElectron]
   Bool_t          Electron_isLowPt[7];   //[nElectron]
   Bool_t          Electron_isPF[7];   //[nElectron]
   Bool_t          Electron_isPFoverlap[7];   //[nElectron]
   UInt_t          nelectronPair;
   Float_t         electronPair_l1_SC_Energy[6];   //[nelectronPair]
   Float_t         electronPair_l1_SC_RegEnergy[6];   //[nelectronPair]
   Float_t         electronPair_l1_SC_rawEnergy[6];   //[nelectronPair]
   Float_t         electronPair_l1_SCeta[6];   //[nelectronPair]
   Float_t         electronPair_l1_SCphi[6];   //[nelectronPair]
   Float_t         electronPair_l2_SC_Energy[6];   //[nelectronPair]
   Float_t         electronPair_l2_SC_RegEnergy[6];   //[nelectronPair]
   Float_t         electronPair_l2_SC_rawEnergy[6];   //[nelectronPair]
   Float_t         electronPair_l2_SCeta[6];   //[nelectronPair]
   Float_t         electronPair_l2_SCphi[6];   //[nelectronPair]
   Float_t         electronPair_mll_llfit[6];   //[nelectronPair]
   Float_t         electronPair_mll_raw[6];   //[nelectronPair]
   Int_t           electronPair_l1Idx[6];   //[nelectronPair]
   Int_t           electronPair_l2Idx[6];   //[nelectronPair]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
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
   TBranch        *b_nelectronPair;   //!
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
   TBranch        *b_electronPair_mll_llfit;   //!
   TBranch        *b_electronPair_mll_raw;   //!
   TBranch        *b_electronPair_l1Idx;   //!
   TBranch        *b_electronPair_l2Idx;   //!

   virtual void     Init(TTree *tree);
};

#endif

