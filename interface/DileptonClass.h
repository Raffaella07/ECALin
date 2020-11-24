//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 18 11:51:38 2020 by ROOT version 6.12/07
// from TTree Events/Events
// found on file: /eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Nov17_withB/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_AmcatNLO_EGamma/201117_095235/0000/BParkNANO_data_2020Nov17_withB_30.root
//////////////////////////////////////////////////////////

#ifndef DileptonClass_h
#define DileptonClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class DileptonClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          nElectron;
   Float_t         Electron_dxy[5];   //[nElectron]
   Float_t         Electron_dxyErr[5];   //[nElectron]
   Float_t         Electron_dz[5];   //[nElectron]
   Float_t         Electron_dzErr[5];   //[nElectron]
   Float_t         Electron_energySC[5];   //[nElectron]
   Float_t         Electron_eta[5];   //[nElectron]
   Float_t         Electron_fBrem[5];   //[nElectron]
   Float_t         Electron_ip3d[5];   //[nElectron]
   Float_t         Electron_mass[5];   //[nElectron]
   Float_t         Electron_mvaId[5];   //[nElectron]
   Float_t         Electron_pfRelIso[5];   //[nElectron]
   Float_t         Electron_pfmvaId[5];   //[nElectron]
   Float_t         Electron_phi[5];   //[nElectron]
   Float_t         Electron_pt[5];   //[nElectron]
   Float_t         Electron_ptBiased[5];   //[nElectron]
   Float_t         Electron_sip3d[5];   //[nElectron]
   Float_t         Electron_trkRelIso[5];   //[nElectron]
   Float_t         Electron_unBiased[5];   //[nElectron]
   Float_t         Electron_vx[5];   //[nElectron]
   Float_t         Electron_vy[5];   //[nElectron]
   Float_t         Electron_vz[5];   //[nElectron]
   Int_t           Electron_charge[5];   //[nElectron]
   Int_t           Electron_pdgId[5];   //[nElectron]
   Bool_t          Electron_convVeto[5];   //[nElectron]
   Bool_t          Electron_isLowPt[5];   //[nElectron]
   Bool_t          Electron_isPF[5];   //[nElectron]
   Bool_t          Electron_isPFoverlap[5];   //[nElectron]
   UInt_t          nelectronPair;
   Float_t         electronPair_l1_SC_Energy[6];   //[nelectronPair]
   Float_t         electronPair_l1_SC_RegEnergy[6];   //[nelectronPair]
   Float_t         electronPair_l1_SC_rawEnergy[6];   //[nelectronPair]
   Float_t         electronPair_l1_SCeta[6];   //[nelectronPair]
   Float_t         electronPair_l1_SCphi[6];   //[nelectronPair]
   Float_t         electronPair_l1_SmearE[6];   //[nelectronPair]
   Float_t         electronPair_l2_SC_Energy[6];   //[nelectronPair]
   Float_t         electronPair_l2_SC_RegEnergy[6];   //[nelectronPair]
   Float_t         electronPair_l2_SC_rawEnergy[6];   //[nelectronPair]
   Float_t         electronPair_l2_SCeta[6];   //[nelectronPair]
   Float_t         electronPair_l2_SCphi[6];   //[nelectronPair]
   Float_t         electronPair_l2_SmearE[6];   //[nelectronPair]
   Float_t         electronPair_mll_llfit[6];   //[nelectronPair]
   Float_t         electronPair_mll_raw[6];   //[nelectronPair]
   Int_t           electronPair_l1Idx[6];   //[nelectronPair]
   Int_t           electronPair_l2Idx[6];   //[nelectronPair]
   UInt_t          nOtherPV;
   Float_t         OtherPV_z[3];   //[nOtherPV]
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   UInt_t          nSV;
   Float_t         SV_dlen[14];   //[nSV]
   Float_t         SV_dlenSig[14];   //[nSV]
   Float_t         SV_pAngle[14];   //[nSV]
   Float_t         SV_chi2[14];   //[nSV]
   Float_t         SV_eta[14];   //[nSV]
   Float_t         SV_mass[14];   //[nSV]
   Float_t         SV_ndof[14];   //[nSV]
   Float_t         SV_phi[14];   //[nSV]
   Float_t         SV_pt[14];   //[nSV]
   Float_t         SV_x[14];   //[nSV]
   Float_t         SV_y[14];   //[nSV]
   Float_t         SV_z[14];   //[nSV]

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
   TBranch        *b_electronPair_l1_SmearE;   //!
   TBranch        *b_electronPair_l2_SC_Energy;   //!
   TBranch        *b_electronPair_l2_SC_RegEnergy;   //!
   TBranch        *b_electronPair_l2_SC_rawEnergy;   //!
   TBranch        *b_electronPair_l2_SCeta;   //!
   TBranch        *b_electronPair_l2_SCphi;   //!
   TBranch        *b_electronPair_l2_SmearE;   //!
   TBranch        *b_electronPair_mll_llfit;   //!
   TBranch        *b_electronPair_mll_raw;   //!
   TBranch        *b_electronPair_l1Idx;   //!
   TBranch        *b_electronPair_l2Idx;   //!
   TBranch        *b_nOtherPV;   //!
   TBranch        *b_OtherPV_z;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_nSV;   //!
   TBranch        *b_SV_dlen;   //!
   TBranch        *b_SV_dlenSig;   //!
   TBranch        *b_SV_pAngle;   //!
   TBranch        *b_SV_chi2;   //!
   TBranch        *b_SV_eta;   //!
   TBranch        *b_SV_mass;   //!
   TBranch        *b_SV_ndof;   //!
   TBranch        *b_SV_phi;   //!
   TBranch        *b_SV_pt;   //!
   TBranch        *b_SV_x;   //!
   TBranch        *b_SV_y;   //!
   TBranch        *b_SV_z;   //!

   virtual void     Init(TTree *tree);
};

#endif

