#include "../interface/MeeClass.h"
void MeeClass::InitTree(TTree *tree)
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

   fChain->SetBranchAddress("run", &run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock);
   fChain->SetBranchAddress("event", &event);
   fChain->SetBranchAddress("mll_EE", &mll_EE);
   fChain->SetBranchAddress("mll_ET", &mll_ET);
   fChain->SetBranchAddress("mll_TT", &mll_TT);
   fChain->SetBranchAddress("B_mass", &B_mass);
   fChain->SetBranchAddress("PUweight", &PUweight);
   fChain->SetBranchAddress("deltaPhi", &deltaPhi);
   fChain->SetBranchAddress("deltaEta", &deltaEta);
   fChain->SetBranchAddress("deltaR", &deltaR);

}
void MeeClass::Init()
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   TTree *tree = new TTree("tree","tree for unbinned fit");
   //if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->Branch("run", &run);
   fChain->Branch("luminosityBlock", &luminosityBlock);
   fChain->Branch("event", &event);
   fChain->Branch("mll_EE", &mll_EE);
   fChain->Branch("mll_ET", &mll_ET);
   fChain->Branch("mll_TT", &mll_TT);
   fChain->Branch("B_mass", &B_mass);
   fChain->Branch("PUweight", &PUweight);
   fChain->Branch("deltaPhi", &deltaPhi);
   fChain->Branch("deltaEta", &deltaEta);
   fChain->Branch("deltaR", &deltaR);


}
