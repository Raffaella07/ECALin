//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov  6 10:45:28 2020 by ROOT version 6.12/07
// from TTree Events/Events
// found on file: /eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Nov04_withB/ParkingBPH1/crab_data_Run2018C_part1/201104_184338/0000/BParkNANO_data_2020Nov04_withB_45.root
//////////////////////////////////////////////////////////

#ifndef DiLeptonClass_h
#define DiLeptonClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class DiLeptonClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          nelectronPair;
   Float_t         electronPair_eta[49];   //[nelectronPair]
   Float_t         electronPair_l1_SC_Energy[49];   //[nelectronPair]
   Float_t         electronPair_l1_SC_RegEnergy[49];   //[nelectronPair]
   Float_t         electronPair_l1_SC_rawEnergy[49];   //[nelectronPair]
   Float_t         electronPair_l1_SCeta[49];   //[nelectronPair]
   Float_t         electronPair_l1_SCphi[49];   //[nelectronPair]
   Float_t         electronPair_l2_SC_Energy[49];   //[nelectronPair]
   Float_t         electronPair_l2_SC_RegEnergy[49];   //[nelectronPair]
   Float_t         electronPair_l2_SC_rawEnergy[49];   //[nelectronPair]
   Float_t         electronPair_l2_SCeta[49];   //[nelectronPair]
   Float_t         electronPair_l2_SCphi[49];   //[nelectronPair]
   Float_t         electronPair_mass[49];   //[nelectronPair]
   Float_t         electronPair_mll_llfit[49];   //[nelectronPair]
   Float_t         electronPair_mll_raw[49];   //[nelectronPair]
   Float_t         electronPair_phi[49];   //[nelectronPair]
   Float_t         electronPair_pt[49];   //[nelectronPair]
   Int_t           electronPair_charge[49];   //[nelectronPair]
   Int_t           electronPair_l1Idx[49];   //[nelectronPair]
   Int_t           electronPair_l2Idx[49];   //[nelectronPair]
   Int_t           electronPair_pdgId[49];   //[nelectronPair]
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
   Float_t         SV_dlen[10];   //[nSV]
   Float_t         SV_dlenSig[10];   //[nSV]
   Float_t         SV_pAngle[10];   //[nSV]
   Float_t         SV_chi2[10];   //[nSV]
   Float_t         SV_eta[10];   //[nSV]
   Float_t         SV_mass[10];   //[nSV]
   Float_t         SV_ndof[10];   //[nSV]
   Float_t         SV_phi[10];   //[nSV]
   Float_t         SV_pt[10];   //[nSV]
   Float_t         SV_x[10];   //[nSV]
   Float_t         SV_y[10];   //[nSV]
   Float_t         SV_z[10];   //[nSV]

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



