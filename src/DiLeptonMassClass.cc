#define DiLeptonMassClass_cxx
#include "DiLeptonMassClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void DiLeptonMassClass::Init(TTree *tree)
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
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("nPSWeight", &nPSWeight, &b_nPSWeight);
   fChain->SetBranchAddress("PSWeight", PSWeight, &b_PSWeight);
   fChain->SetBranchAddress("HLT_Mu7_IP4", &HLT_Mu7_IP4, &b_HLT_Mu7_IP4);
   fChain->SetBranchAddress("HLT_Mu8_IP6", &HLT_Mu8_IP6, &b_HLT_Mu8_IP6);
   fChain->SetBranchAddress("HLT_Mu8_IP5", &HLT_Mu8_IP5, &b_HLT_Mu8_IP5);
   fChain->SetBranchAddress("HLT_Mu8p5_IP3p5", &HLT_Mu8p5_IP3p5, &b_HLT_Mu8p5_IP3p5);
   fChain->SetBranchAddress("HLT_Mu9_IP6", &HLT_Mu9_IP6, &b_HLT_Mu9_IP6);
   fChain->SetBranchAddress("HLT_Mu9_IP5", &HLT_Mu9_IP5, &b_HLT_Mu9_IP5);
   fChain->SetBranchAddress("HLT_Mu9_IP4", &HLT_Mu9_IP4, &b_HLT_Mu9_IP4);
   fChain->SetBranchAddress("HLT_Mu10p5_IP3p5", &HLT_Mu10p5_IP3p5, &b_HLT_Mu10p5_IP3p5);
   fChain->SetBranchAddress("HLT_Mu12_IP6", &HLT_Mu12_IP6, &b_HLT_Mu12_IP6);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("nDiLepton", &nDiLepton, &b_nDiLepton);
   fChain->SetBranchAddress("DiLepton_mll_fullfit", DiLepton_mll_fullfit, &b_DiLepton_mll_fullfit);
   fChain->SetBranchAddress("DiLepton_l1Idx", DiLepton_l1Idx, &b_DiLepton_l1Idx);
   fChain->SetBranchAddress("DiLepton_l2Idx", DiLepton_l2Idx, &b_DiLepton_l2Idx);
   fChain->SetBranchAddress("DiLepton_l1isPF", DiLepton_l1isPF, &b_DiLepton_l1isPF);
   fChain->SetBranchAddress("DiLepton_l1_pt", DiLepton_l1_pt, &b_DiLepton_l1_pt);
   fChain->SetBranchAddress("DiLepton_l1_eta", DiLepton_l1_eta, &b_DiLepton_l1_eta);
   fChain->SetBranchAddress("DiLepton_l1_phi", DiLepton_l1_phi, &b_DiLepton_l1_phi);
   fChain->SetBranchAddress("DiLepton_l2isPF", DiLepton_l2isPF, &b_DiLepton_l2isPF);
   fChain->SetBranchAddress("DiLepton_l2_pt", DiLepton_l2_pt, &b_DiLepton_l2_pt);
   fChain->SetBranchAddress("DiLepton_l2_eta", DiLepton_l2_eta, &b_DiLepton_l2_eta);
   fChain->SetBranchAddress("DiLepton_l2_phi", DiLepton_l2_phi, &b_DiLepton_l2_phi);
}

