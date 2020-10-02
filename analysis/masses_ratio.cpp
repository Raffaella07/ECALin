
#include "../interface/BParkTools.h"
#include "TTree.h"
//#include <ROOT/RDataFrame.hxx>
//#include <ROOT/RVec.hxx>
//#include "TStopwatch.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <math.h>
#include <string>
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TAxis.h"
#include "TTree.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TString.h"
#include "TLorentzVector.h"



void FillHistos(DiLeptonMassClass* sgl,TH2D* hist);


int main(){

	std::ifstream bkg_input("pattern.txt");
	std::string data_path = "/eos/home-r/ratramon/ECAL_lin/Run2018D_part1.root";
	std::string mc_path = "/eos/home-r/ratramon/ECAL_lin/crab_BuToKJpsi_ToEE__ext_0000.root";
/*	std::string bkg_path = "/eos/home-r/ratramon/BToKEE/";
	if(!bkg_input) {
		std::cout << "Cannot open input file." << std::endl ;
		return 1;
	}*/

	TChain *mc = new TChain("Events");
	TChain *data = new TChain("Events");

	mc->Add(mc_path.c_str());
	data->Add(data_path.c_str());

	//std::cout << "Loaded " << mc->GetEntries() << " entries for signal" << std::endl;
	//std::cout << "Loaded " << data->GetEntries() << " entries for signal" << std::endl;
	//std::cout << "Loaded " << background->GetEntries() << " entries for background" << std::endl;

	DiLeptonMassClass sgl_mc;
	DiLeptonMassClass sgl_data;
	sgl_mc.Init(mc);
	sgl_data.Init(data);
	
	int nPtbin = 10;	
	int i;
	double min_pt = 0;
	double max_pt = 30;
	

	TH2D* Dimass_pt_mc = new TH2D("diLepton mass vs lepton pt ","diLepton mass vs lepton pt ",50,2.5,3.5,10,min_pt,max_pt);
	TH2D* Dimass_pt_data = new TH2D("diLepton mass vs lepton pt ","diLepton mass vs lepton pt ",50,2.5,3.5,10,min_pt,max_pt);
	


	
	
	FillHistos(&sgl_mc,Dimass_pt_mc);	
	FillHistos(&sgl_data,Dimass_pt_data);	
	
	SavePlot2D("m_{e^{+}e^{-}}(GeV)","electron_pt (GeV)",Dimass_pt_mc,"Dimass_pt_mc",0,0);
	SavePlot2D("m_{e^{+}e^{-}}(GeV)","electron_pt (GeV)",Dimass_pt_data,"Dimass_pt_data",0,0);
	
	DoubleSlicer("/plots/",4,min_pt, max_pt,"m_{e^{+}e^{-}}",Dimass_pt_mc,Dimass_pt_data,"LeptonPairMassSlice",1);	
	delete Dimass_pt_mc;	
	delete Dimass_pt_data;	
	return 0;
}




void FillHistos(DiLeptonMassClass* sgl,TH2D* hist){

	int i;


	std::cout << "Loaded " << sgl->fChain->GetEntries() << " entries for signal" << std::endl;
	
	for(i=0; i<sgl->fChain->GetEntries();i++){
		
		if (i%10000==0 ) std::cout << "at entry" << i << std::endl;
		sgl->fChain->GetEntry(i);	

		for (int j=0; j<sgl->nDiLepton; j++){
			
			 
	       		 if (fabs(sgl->DiLepton_l1_eta[j])<1.5) hist->Fill(sgl->DiLepton_mll_fullfit[j],sgl->DiLepton_l1_pt[j]);
	       		 if (fabs(sgl->DiLepton_l2_eta[j])<1.5) hist->Fill(sgl->DiLepton_mll_fullfit[j],sgl->DiLepton_l2_pt[j]);
	
		}
	}






}
