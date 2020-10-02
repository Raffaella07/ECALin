//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct  1 18:52:22 2020 by ROOT version 6.12/07
// from TTree Events/Events
// found on file: /eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Oct01_withB/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_Toee_ext/201001_141013/0000/BParkNANO_mc_2020Oct01_withB_233.root
//////////////////////////////////////////////////////////

#ifndef DiLeptonFromBClass_h
#define DiLeptonFromBClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class DiLeptonFromBClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          nBToKEE;
   Float_t         BToKEE_b_iso03[68];   //[nBToKEE]
   Float_t         BToKEE_b_iso04[68];   //[nBToKEE]
   Float_t         BToKEE_cos2D[68];   //[nBToKEE]
   Float_t         BToKEE_eta[68];   //[nBToKEE]
   Float_t         BToKEE_fit_cos2D[68];   //[nBToKEE]
   Float_t         BToKEE_fit_eta[68];   //[nBToKEE]
   Float_t         BToKEE_fit_k_eta[68];   //[nBToKEE]
   Float_t         BToKEE_fit_k_phi[68];   //[nBToKEE]
   Float_t         BToKEE_fit_k_pt[68];   //[nBToKEE]
   Float_t         BToKEE_fit_l1_eta[68];   //[nBToKEE]
   Float_t         BToKEE_fit_l1_phi[68];   //[nBToKEE]
   Float_t         BToKEE_fit_l1_pt[68];   //[nBToKEE]
   Float_t         BToKEE_fit_l2_eta[68];   //[nBToKEE]
   Float_t         BToKEE_fit_l2_phi[68];   //[nBToKEE]
   Float_t         BToKEE_fit_l2_pt[68];   //[nBToKEE]
   Float_t         BToKEE_fit_mass[68];   //[nBToKEE]
   Float_t         BToKEE_fit_massErr[68];   //[nBToKEE]
   Float_t         BToKEE_fit_phi[68];   //[nBToKEE]
   Float_t         BToKEE_fit_pt[68];   //[nBToKEE]
   Float_t         BToKEE_k_iso03[68];   //[nBToKEE]
   Float_t         BToKEE_k_iso04[68];   //[nBToKEE]
   Float_t         BToKEE_l1_SC_EOverP[68];   //[nBToKEE]
   Float_t         BToKEE_l1_SC_Energy[68];   //[nBToKEE]
   Float_t         BToKEE_l1_SC_RegEnergy[68];   //[nBToKEE]
   Float_t         BToKEE_l1_SC_rawEnergy[68];   //[nBToKEE]
   Float_t         BToKEE_l1_SCeta[68];   //[nBToKEE]
   Float_t         BToKEE_l1_SCphi[68];   //[nBToKEE]
   Float_t         BToKEE_l1_iso03[68];   //[nBToKEE]
   Float_t         BToKEE_l1_iso04[68];   //[nBToKEE]
   Float_t         BToKEE_l2_SC_EOverP[68];   //[nBToKEE]
   Float_t         BToKEE_l2_SC_Energy[68];   //[nBToKEE]
   Float_t         BToKEE_l2_SC_RegEnergy[68];   //[nBToKEE]
   Float_t         BToKEE_l2_SC_rawEnergy[68];   //[nBToKEE]
   Float_t         BToKEE_l2_SCeta[68];   //[nBToKEE]
   Float_t         BToKEE_l2_SCphi[68];   //[nBToKEE]
   Float_t         BToKEE_l2_iso03[68];   //[nBToKEE]
   Float_t         BToKEE_l2_iso04[68];   //[nBToKEE]
   Float_t         BToKEE_l_xy[68];   //[nBToKEE]
   Float_t         BToKEE_l_xy_unc[68];   //[nBToKEE]
   Float_t         BToKEE_mass[68];   //[nBToKEE]
   Float_t         BToKEE_maxDR[68];   //[nBToKEE]
   Float_t         BToKEE_minDR[68];   //[nBToKEE]
   Float_t         BToKEE_mllErr_llfit[68];   //[nBToKEE]
   Float_t         BToKEE_mll_fullfit[68];   //[nBToKEE]
   Float_t         BToKEE_mll_llfit[68];   //[nBToKEE]
   Float_t         BToKEE_mll_raw[68];   //[nBToKEE]
   Float_t         BToKEE_phi[68];   //[nBToKEE]
   Float_t         BToKEE_pt[68];   //[nBToKEE]
   Float_t         BToKEE_svprob[68];   //[nBToKEE]
   Float_t         BToKEE_vtx_ex[68];   //[nBToKEE]
   Float_t         BToKEE_vtx_ey[68];   //[nBToKEE]
   Float_t         BToKEE_vtx_ez[68];   //[nBToKEE]
   Float_t         BToKEE_vtx_x[68];   //[nBToKEE]
   Float_t         BToKEE_vtx_y[68];   //[nBToKEE]
   Float_t         BToKEE_vtx_z[68];   //[nBToKEE]
   Int_t           BToKEE_charge[68];   //[nBToKEE]
   Int_t           BToKEE_kIdx[68];   //[nBToKEE]
   Int_t           BToKEE_l1Idx[68];   //[nBToKEE]
   Int_t           BToKEE_l2Idx[68];   //[nBToKEE]
   Int_t           BToKEE_n_k_used[68];   //[nBToKEE]
   Int_t           BToKEE_n_l1_used[68];   //[nBToKEE]
   Int_t           BToKEE_n_l2_used[68];   //[nBToKEE]
   Int_t           BToKEE_pdgId[68];   //[nBToKEE]
   UInt_t          nElectron;
   Float_t         Electron_dxy[4];   //[nElectron]
   Float_t         Electron_dxyErr[4];   //[nElectron]
   Float_t         Electron_dz[4];   //[nElectron]
   Float_t         Electron_dzErr[4];   //[nElectron]
   Float_t         Electron_energySC[4];   //[nElectron]
   Float_t         Electron_eta[4];   //[nElectron]
   Float_t         Electron_fBrem[4];   //[nElectron]
   Float_t         Electron_ip3d[4];   //[nElectron]
   Float_t         Electron_mass[4];   //[nElectron]
   Float_t         Electron_mvaId[4];   //[nElectron]
   Float_t         Electron_pfRelIso[4];   //[nElectron]
   Float_t         Electron_pfmvaId[4];   //[nElectron]
   Float_t         Electron_phi[4];   //[nElectron]
   Float_t         Electron_pt[4];   //[nElectron]
   Float_t         Electron_ptBiased[4];   //[nElectron]
   Float_t         Electron_sip3d[4];   //[nElectron]
   Float_t         Electron_trkRelIso[4];   //[nElectron]
   Float_t         Electron_unBiased[4];   //[nElectron]
   Float_t         Electron_vx[4];   //[nElectron]
   Float_t         Electron_vy[4];   //[nElectron]
   Float_t         Electron_vz[4];   //[nElectron]
   Int_t           Electron_charge[4];   //[nElectron]
   Int_t           Electron_pdgId[4];   //[nElectron]
   Bool_t          Electron_convVeto[4];   //[nElectron]
   Bool_t          Electron_isLowPt[4];   //[nElectron]
   Bool_t          Electron_isPF[4];   //[nElectron]
   Bool_t          Electron_isPFoverlap[4];   //[nElectron]
   UInt_t          nGenPart;
   Float_t         GenPart_eta[113];   //[nGenPart]
   Float_t         GenPart_mass[113];   //[nGenPart]
   Float_t         GenPart_phi[113];   //[nGenPart]
   Float_t         GenPart_pt[113];   //[nGenPart]
   Float_t         GenPart_vx[113];   //[nGenPart]
   Float_t         GenPart_vy[113];   //[nGenPart]
   Float_t         GenPart_vz[113];   //[nGenPart]
   Int_t           GenPart_genPartIdxMother[113];   //[nGenPart]
   Int_t           GenPart_pdgId[113];   //[nGenPart]
   Int_t           GenPart_status[113];   //[nGenPart]
   Int_t           GenPart_statusFlags[113];   //[nGenPart]
   Float_t         Generator_binvar;
   Float_t         Generator_scalePDF;
   Float_t         Generator_weight;
   Float_t         Generator_x1;
   Float_t         Generator_x2;
   Float_t         Generator_xpdf1;
   Float_t         Generator_xpdf2;
   Int_t           Generator_id1;
   Int_t           Generator_id2;
   Float_t         genWeight;
   UInt_t          nPSWeight;
   Float_t         PSWeight[1];   //[nPSWeight]
   UInt_t          nMuon;
   Float_t         Muon_dxy[7];   //[nMuon]
   Float_t         Muon_dxyErr[7];   //[nMuon]
   Float_t         Muon_dz[7];   //[nMuon]
   Float_t         Muon_dzErr[7];   //[nMuon]
   Float_t         Muon_eta[7];   //[nMuon]
   Float_t         Muon_ip3d[7];   //[nMuon]
   Float_t         Muon_mass[7];   //[nMuon]
   Float_t         Muon_pfRelIso03_all[7];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[7];   //[nMuon]
   Float_t         Muon_phi[7];   //[nMuon]
   Float_t         Muon_pt[7];   //[nMuon]
   Float_t         Muon_ptErr[7];   //[nMuon]
   Float_t         Muon_sip3d[7];   //[nMuon]
   Float_t         Muon_vx[7];   //[nMuon]
   Float_t         Muon_vy[7];   //[nMuon]
   Float_t         Muon_vz[7];   //[nMuon]
   Int_t           Muon_charge[7];   //[nMuon]
   Int_t           Muon_isTriggering[7];   //[nMuon]
   Int_t           Muon_pdgId[7];   //[nMuon]
   Bool_t          Muon_isGlobal[7];   //[nMuon]
   Bool_t          Muon_isPFcand[7];   //[nMuon]
   Bool_t          Muon_isTracker[7];   //[nMuon]
   Bool_t          Muon_mediumId[7];   //[nMuon]
   UChar_t         Muon_pfIsoId[7];   //[nMuon]
   Bool_t          Muon_softId[7];   //[nMuon]
   Bool_t          Muon_tightId[7];   //[nMuon]
   UChar_t         Muon_tkIsoId[7];   //[nMuon]
   Bool_t          Muon_triggerIdLoose[7];   //[nMuon]
   UInt_t          nTriggerMuon;
   Float_t         TriggerMuon_eta[2];   //[nTriggerMuon]
   Float_t         TriggerMuon_mass[2];   //[nTriggerMuon]
   Float_t         TriggerMuon_phi[2];   //[nTriggerMuon]
   Float_t         TriggerMuon_pt[2];   //[nTriggerMuon]
   Float_t         TriggerMuon_vx[2];   //[nTriggerMuon]
   Float_t         TriggerMuon_vy[2];   //[nTriggerMuon]
   Float_t         TriggerMuon_vz[2];   //[nTriggerMuon]
   Int_t           TriggerMuon_charge[2];   //[nTriggerMuon]
   Int_t           TriggerMuon_pdgId[2];   //[nTriggerMuon]
   Int_t           TriggerMuon_trgMuonIndex[2];   //[nTriggerMuon]
   Float_t         Pileup_nTrueInt;
   Float_t         Pileup_pudensity;
   Float_t         Pileup_gpudensity;
   Int_t           Pileup_nPU;
   Int_t           Pileup_sumEOOT;
   Int_t           Pileup_sumLOOT;
   Float_t         fixedGridRhoFastjetAll;
   Float_t         fixedGridRhoFastjetCentral;
   Float_t         fixedGridRhoFastjetCentralCalo;
   Float_t         fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         fixedGridRhoFastjetCentralNeutral;
   UInt_t          nProbeTracks;
   Float_t         ProbeTracks_DCASig[308];   //[nProbeTracks]
   Float_t         ProbeTracks_dxy[308];   //[nProbeTracks]
   Float_t         ProbeTracks_dxyS[308];   //[nProbeTracks]
   Float_t         ProbeTracks_dz[308];   //[nProbeTracks]
   Float_t         ProbeTracks_dzS[308];   //[nProbeTracks]
   Float_t         ProbeTracks_eta[308];   //[nProbeTracks]
   Float_t         ProbeTracks_mass[308];   //[nProbeTracks]
   Float_t         ProbeTracks_phi[308];   //[nProbeTracks]
   Float_t         ProbeTracks_pt[308];   //[nProbeTracks]
   Float_t         ProbeTracks_vx[308];   //[nProbeTracks]
   Float_t         ProbeTracks_vy[308];   //[nProbeTracks]
   Float_t         ProbeTracks_vz[308];   //[nProbeTracks]
   Int_t           ProbeTracks_charge[308];   //[nProbeTracks]
   Int_t           ProbeTracks_isLostTrk[308];   //[nProbeTracks]
   Int_t           ProbeTracks_isPacked[308];   //[nProbeTracks]
   Int_t           ProbeTracks_nValidHits[308];   //[nProbeTracks]
   Int_t           ProbeTracks_pdgId[308];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToEle[308];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToLooseMuon[308];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToMediumMuon[308];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToMuon[308];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToSoftMuon[308];   //[nProbeTracks]
   UInt_t          nTrigObj;
   Float_t         TrigObj_pt[3];   //[nTrigObj]
   Float_t         TrigObj_eta[3];   //[nTrigObj]
   Float_t         TrigObj_phi[3];   //[nTrigObj]
   Float_t         TrigObj_l1pt[3];   //[nTrigObj]
   Float_t         TrigObj_l1pt_2[3];   //[nTrigObj]
   Float_t         TrigObj_l2pt[3];   //[nTrigObj]
   Int_t           TrigObj_id[3];   //[nTrigObj]
   Int_t           TrigObj_l1iso[3];   //[nTrigObj]
   Int_t           TrigObj_l1charge[3];   //[nTrigObj]
   Int_t           TrigObj_filterBits[3];   //[nTrigObj]
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
   Float_t         SV_dlen[7];   //[nSV]
   Float_t         SV_dlenSig[7];   //[nSV]
   Float_t         SV_pAngle[7];   //[nSV]
   Int_t           Electron_genPartIdx[4];   //[nElectron]
   Int_t           Electron_genPartFlav[4];   //[nElectron]
   Int_t           Muon_genPartIdx[7];   //[nMuon]
   Int_t           Muon_genPartFlav[7];   //[nMuon]
   Float_t         SV_chi2[7];   //[nSV]
   Float_t         SV_eta[7];   //[nSV]
   Float_t         SV_mass[7];   //[nSV]
   Float_t         SV_ndof[7];   //[nSV]
   Float_t         SV_phi[7];   //[nSV]
   Float_t         SV_pt[7];   //[nSV]
   Float_t         SV_x[7];   //[nSV]
   Float_t         SV_y[7];   //[nSV]
   Float_t         SV_z[7];   //[nSV]
   Int_t           ProbeTracks_genPartIdx[308];   //[nProbeTracks]
   Int_t           ProbeTracks_genPartFlav[308];   //[nProbeTracks]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nBToKEE;   //!
   TBranch        *b_BToKEE_b_iso03;   //!
   TBranch        *b_BToKEE_b_iso04;   //!
   TBranch        *b_BToKEE_cos2D;   //!
   TBranch        *b_BToKEE_eta;   //!
   TBranch        *b_BToKEE_fit_cos2D;   //!
   TBranch        *b_BToKEE_fit_eta;   //!
   TBranch        *b_BToKEE_fit_k_eta;   //!
   TBranch        *b_BToKEE_fit_k_phi;   //!
   TBranch        *b_BToKEE_fit_k_pt;   //!
   TBranch        *b_BToKEE_fit_l1_eta;   //!
   TBranch        *b_BToKEE_fit_l1_phi;   //!
   TBranch        *b_BToKEE_fit_l1_pt;   //!
   TBranch        *b_BToKEE_fit_l2_eta;   //!
   TBranch        *b_BToKEE_fit_l2_phi;   //!
   TBranch        *b_BToKEE_fit_l2_pt;   //!
   TBranch        *b_BToKEE_fit_mass;   //!
   TBranch        *b_BToKEE_fit_massErr;   //!
   TBranch        *b_BToKEE_fit_phi;   //!
   TBranch        *b_BToKEE_fit_pt;   //!
   TBranch        *b_BToKEE_k_iso03;   //!
   TBranch        *b_BToKEE_k_iso04;   //!
   TBranch        *b_BToKEE_l1_SC_EOverP;   //!
   TBranch        *b_BToKEE_l1_SC_Energy;   //!
   TBranch        *b_BToKEE_l1_SC_RegEnergy;   //!
   TBranch        *b_BToKEE_l1_SC_rawEnergy;   //!
   TBranch        *b_BToKEE_l1_SCeta;   //!
   TBranch        *b_BToKEE_l1_SCphi;   //!
   TBranch        *b_BToKEE_l1_iso03;   //!
   TBranch        *b_BToKEE_l1_iso04;   //!
   TBranch        *b_BToKEE_l2_SC_EOverP;   //!
   TBranch        *b_BToKEE_l2_SC_Energy;   //!
   TBranch        *b_BToKEE_l2_SC_RegEnergy;   //!
   TBranch        *b_BToKEE_l2_SC_rawEnergy;   //!
   TBranch        *b_BToKEE_l2_SCeta;   //!
   TBranch        *b_BToKEE_l2_SCphi;   //!
   TBranch        *b_BToKEE_l2_iso03;   //!
   TBranch        *b_BToKEE_l2_iso04;   //!
   TBranch        *b_BToKEE_l_xy;   //!
   TBranch        *b_BToKEE_l_xy_unc;   //!
   TBranch        *b_BToKEE_mass;   //!
   TBranch        *b_BToKEE_maxDR;   //!
   TBranch        *b_BToKEE_minDR;   //!
   TBranch        *b_BToKEE_mllErr_llfit;   //!
   TBranch        *b_BToKEE_mll_fullfit;   //!
   TBranch        *b_BToKEE_mll_llfit;   //!
   TBranch        *b_BToKEE_mll_raw;   //!
   TBranch        *b_BToKEE_phi;   //!
   TBranch        *b_BToKEE_pt;   //!
   TBranch        *b_BToKEE_svprob;   //!
   TBranch        *b_BToKEE_vtx_ex;   //!
   TBranch        *b_BToKEE_vtx_ey;   //!
   TBranch        *b_BToKEE_vtx_ez;   //!
   TBranch        *b_BToKEE_vtx_x;   //!
   TBranch        *b_BToKEE_vtx_y;   //!
   TBranch        *b_BToKEE_vtx_z;   //!
   TBranch        *b_BToKEE_charge;   //!
   TBranch        *b_BToKEE_kIdx;   //!
   TBranch        *b_BToKEE_l1Idx;   //!
   TBranch        *b_BToKEE_l2Idx;   //!
   TBranch        *b_BToKEE_n_k_used;   //!
   TBranch        *b_BToKEE_n_l1_used;   //!
   TBranch        *b_BToKEE_n_l2_used;   //!
   TBranch        *b_BToKEE_pdgId;   //!
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
   TBranch        *b_nGenPart;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_mass;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_vx;   //!
   TBranch        *b_GenPart_vy;   //!
   TBranch        *b_GenPart_vz;   //!
   TBranch        *b_GenPart_genPartIdxMother;   //!
   TBranch        *b_GenPart_pdgId;   //!
   TBranch        *b_GenPart_status;   //!
   TBranch        *b_GenPart_statusFlags;   //!
   TBranch        *b_Generator_binvar;   //!
   TBranch        *b_Generator_scalePDF;   //!
   TBranch        *b_Generator_weight;   //!
   TBranch        *b_Generator_x1;   //!
   TBranch        *b_Generator_x2;   //!
   TBranch        *b_Generator_xpdf1;   //!
   TBranch        *b_Generator_xpdf2;   //!
   TBranch        *b_Generator_id1;   //!
   TBranch        *b_Generator_id2;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_nPSWeight;   //!
   TBranch        *b_PSWeight;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_dxyErr;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_dzErr;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_ip3d;   //!
   TBranch        *b_Muon_mass;   //!
   TBranch        *b_Muon_pfRelIso03_all;   //!
   TBranch        *b_Muon_pfRelIso04_all;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_ptErr;   //!
   TBranch        *b_Muon_sip3d;   //!
   TBranch        *b_Muon_vx;   //!
   TBranch        *b_Muon_vy;   //!
   TBranch        *b_Muon_vz;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_isTriggering;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_isPFcand;   //!
   TBranch        *b_Muon_isTracker;   //!
   TBranch        *b_Muon_mediumId;   //!
   TBranch        *b_Muon_pfIsoId;   //!
   TBranch        *b_Muon_softId;   //!
   TBranch        *b_Muon_tightId;   //!
   TBranch        *b_Muon_tkIsoId;   //!
   TBranch        *b_Muon_triggerIdLoose;   //!
   TBranch        *b_nTriggerMuon;   //!
   TBranch        *b_TriggerMuon_eta;   //!
   TBranch        *b_TriggerMuon_mass;   //!
   TBranch        *b_TriggerMuon_phi;   //!
   TBranch        *b_TriggerMuon_pt;   //!
   TBranch        *b_TriggerMuon_vx;   //!
   TBranch        *b_TriggerMuon_vy;   //!
   TBranch        *b_TriggerMuon_vz;   //!
   TBranch        *b_TriggerMuon_charge;   //!
   TBranch        *b_TriggerMuon_pdgId;   //!
   TBranch        *b_TriggerMuon_trgMuonIndex;   //!
   TBranch        *b_Pileup_nTrueInt;   //!
   TBranch        *b_Pileup_pudensity;   //!
   TBranch        *b_Pileup_gpudensity;   //!
   TBranch        *b_Pileup_nPU;   //!
   TBranch        *b_Pileup_sumEOOT;   //!
   TBranch        *b_Pileup_sumLOOT;   //!
   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_fixedGridRhoFastjetCentral;   //!
   TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_nProbeTracks;   //!
   TBranch        *b_ProbeTracks_DCASig;   //!
   TBranch        *b_ProbeTracks_dxy;   //!
   TBranch        *b_ProbeTracks_dxyS;   //!
   TBranch        *b_ProbeTracks_dz;   //!
   TBranch        *b_ProbeTracks_dzS;   //!
   TBranch        *b_ProbeTracks_eta;   //!
   TBranch        *b_ProbeTracks_mass;   //!
   TBranch        *b_ProbeTracks_phi;   //!
   TBranch        *b_ProbeTracks_pt;   //!
   TBranch        *b_ProbeTracks_vx;   //!
   TBranch        *b_ProbeTracks_vy;   //!
   TBranch        *b_ProbeTracks_vz;   //!
   TBranch        *b_ProbeTracks_charge;   //!
   TBranch        *b_ProbeTracks_isLostTrk;   //!
   TBranch        *b_ProbeTracks_isPacked;   //!
   TBranch        *b_ProbeTracks_nValidHits;   //!
   TBranch        *b_ProbeTracks_pdgId;   //!
   TBranch        *b_ProbeTracks_isMatchedToEle;   //!
   TBranch        *b_ProbeTracks_isMatchedToLooseMuon;   //!
   TBranch        *b_ProbeTracks_isMatchedToMediumMuon;   //!
   TBranch        *b_ProbeTracks_isMatchedToMuon;   //!
   TBranch        *b_ProbeTracks_isMatchedToSoftMuon;   //!
   TBranch        *b_nTrigObj;   //!
   TBranch        *b_TrigObj_pt;   //!
   TBranch        *b_TrigObj_eta;   //!
   TBranch        *b_TrigObj_phi;   //!
   TBranch        *b_TrigObj_l1pt;   //!
   TBranch        *b_TrigObj_l1pt_2;   //!
   TBranch        *b_TrigObj_l2pt;   //!
   TBranch        *b_TrigObj_id;   //!
   TBranch        *b_TrigObj_l1iso;   //!
   TBranch        *b_TrigObj_l1charge;   //!
   TBranch        *b_TrigObj_filterBits;   //!
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
   TBranch        *b_Electron_genPartIdx;   //!
   TBranch        *b_Electron_genPartFlav;   //!
   TBranch        *b_Muon_genPartIdx;   //!
   TBranch        *b_Muon_genPartFlav;   //!
   TBranch        *b_SV_chi2;   //!
   TBranch        *b_SV_eta;   //!
   TBranch        *b_SV_mass;   //!
   TBranch        *b_SV_ndof;   //!
   TBranch        *b_SV_phi;   //!
   TBranch        *b_SV_pt;   //!
   TBranch        *b_SV_x;   //!
   TBranch        *b_SV_y;   //!
   TBranch        *b_SV_z;   //!
   TBranch        *b_ProbeTracks_genPartIdx;   //!
   TBranch        *b_ProbeTracks_genPartFlav;   //!

   virtual void     Init(TTree *tree);
};

#endif


