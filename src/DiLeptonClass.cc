#include "../interface/DiLeptonClass.h"

void DiLeptonClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nelectronPair", &nelectronPair, &b_nelectronPair);
   fChain->SetBranchAddress("electronPair_eta", electronPair_eta, &b_electronPair_eta);
   fChain->SetBranchAddress("electronPair_l1_SC_Energy", electronPair_l1_SC_Energy, &b_electronPair_l1_SC_Energy);
   fChain->SetBranchAddress("electronPair_l1_SC_RegEnergy", electronPair_l1_SC_RegEnergy, &b_electronPair_l1_SC_RegEnergy);
   fChain->SetBranchAddress("electronPair_l1_SC_rawEnergy", electronPair_l1_SC_rawEnergy, &b_electronPair_l1_SC_rawEnergy);
   fChain->SetBranchAddress("electronPair_l1_SCeta", electronPair_l1_SCeta, &b_electronPair_l1_SCeta);
   fChain->SetBranchAddress("electronPair_l1_SCphi", electronPair_l1_SCphi, &b_electronPair_l1_SCphi);
   fChain->SetBranchAddress("electronPair_l2_SC_Energy", electronPair_l2_SC_Energy, &b_electronPair_l2_SC_Energy);
   fChain->SetBranchAddress("electronPair_l2_SC_RegEnergy", electronPair_l2_SC_RegEnergy, &b_electronPair_l2_SC_RegEnergy);
   fChain->SetBranchAddress("electronPair_l2_SC_rawEnergy", electronPair_l2_SC_rawEnergy, &b_electronPair_l2_SC_rawEnergy);
   fChain->SetBranchAddress("electronPair_l2_SCeta", electronPair_l2_SCeta, &b_electronPair_l2_SCeta);
   fChain->SetBranchAddress("electronPair_l2_SCphi", electronPair_l2_SCphi, &b_electronPair_l2_SCphi);
   fChain->SetBranchAddress("electronPair_mass", electronPair_mass, &b_electronPair_mass);
   fChain->SetBranchAddress("electronPair_mll_llfit", electronPair_mll_llfit, &b_electronPair_mll_llfit);
   fChain->SetBranchAddress("electronPair_mll_raw", electronPair_mll_raw, &b_electronPair_mll_raw);
   fChain->SetBranchAddress("electronPair_phi", electronPair_phi, &b_electronPair_phi);
   fChain->SetBranchAddress("electronPair_pt", electronPair_pt, &b_electronPair_pt);
   fChain->SetBranchAddress("electronPair_charge", electronPair_charge, &b_electronPair_charge);
   fChain->SetBranchAddress("electronPair_l1Idx", electronPair_l1Idx, &b_electronPair_l1Idx);
   fChain->SetBranchAddress("electronPair_l2Idx", electronPair_l2Idx, &b_electronPair_l2Idx);
   fChain->SetBranchAddress("electronPair_pdgId", electronPair_pdgId, &b_electronPair_pdgId);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dxyErr", Electron_dxyErr, &b_Electron_dxyErr);
   fChain->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron_dzErr", Electron_dzErr, &b_Electron_dzErr);
   fChain->SetBranchAddress("Electron_energySC", Electron_energySC, &b_Electron_energySC);
   fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_fBrem", Electron_fBrem, &b_Electron_fBrem);
   fChain->SetBranchAddress("Electron_ip3d", Electron_ip3d, &b_Electron_ip3d);
   fChain->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
   fChain->SetBranchAddress("Electron_mvaId", Electron_mvaId, &b_Electron_mvaId);
   fChain->SetBranchAddress("Electron_pfRelIso", Electron_pfRelIso, &b_Electron_pfRelIso);
   fChain->SetBranchAddress("Electron_pfmvaId", Electron_pfmvaId, &b_Electron_pfmvaId);
   fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_ptBiased", Electron_ptBiased, &b_Electron_ptBiased);
   fChain->SetBranchAddress("Electron_sip3d", Electron_sip3d, &b_Electron_sip3d);
   fChain->SetBranchAddress("Electron_trkRelIso", Electron_trkRelIso, &b_Electron_trkRelIso);
   fChain->SetBranchAddress("Electron_unBiased", Electron_unBiased, &b_Electron_unBiased);
   fChain->SetBranchAddress("Electron_vx", Electron_vx, &b_Electron_vx);
   fChain->SetBranchAddress("Electron_vy", Electron_vy, &b_Electron_vy);
   fChain->SetBranchAddress("Electron_vz", Electron_vz, &b_Electron_vz);
   fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_pdgId", Electron_pdgId, &b_Electron_pdgId);
   fChain->SetBranchAddress("Electron_convVeto", Electron_convVeto, &b_Electron_convVeto);
   fChain->SetBranchAddress("Electron_isLowPt", Electron_isLowPt, &b_Electron_isLowPt);
   fChain->SetBranchAddress("Electron_isPF", Electron_isPF, &b_Electron_isPF);
   fChain->SetBranchAddress("Electron_isPFoverlap", Electron_isPFoverlap, &b_Electron_isPFoverlap);
   fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   fChain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   fChain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   fChain->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
   fChain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   fChain->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
   fChain->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
   fChain->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
   fChain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   fChain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   fChain->SetBranchAddress("SV_z", SV_z, &b_SV_z);
}
