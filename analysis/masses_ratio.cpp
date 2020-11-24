
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
#include "TCut.h"
#include "TLorentzVector.h"



void FillHistos(DileptonGeneralClass* sgl,TH2D* hist, int i, bool isMC, bool TwoTracki, bool OuterBarrel);
void FillHistosFromB(DiLeptonFromBClass* sgl,TH2D* hist, int i, bool isMC, bool TwoTrack,bool OuterBarrel);


int main(int argc, char** argv ){

	bool TwoTrack = atoi(argv[1]);
	std::ifstream bkg_input("pattern.txt");
	std::string data_path = "../data/Dielectron_small_*.root";
	std::string mc_path = "/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Oct01_withB/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_Toee_ext/201001_141013/0000/*.root";
/*	std::string bkg_path = "/eos/home-r/ratramon/BToKEE/";
	if(!bkg_input) {
		std::cout << "Cannot open input file." << std::endl ;
		return 1;
	}*/

	TChain *mc = new TChain("Events");
	TChain *data = new TChain("Events");

	mc->Add(mc_path.c_str());
	data->Add((data_path).c_str());
/*	for (int i=10; i<546;i++){
	
	
	}*/

	//std::cout << "Loaded " << mc->GetEntries() << " entries for signal" << std::endl;
	//std::cout << "Loaded " << data->GetEntries() << " entries for signal" << std::endl;
	//std::cout << "Loaded " << background->GetEntries() << " entries for background" << std::endl;

	DiLeptonFromBClass sgl_mc;
	DileptonGeneralClass sgl_data;
	sgl_mc.Init(mc);
	sgl_data.Init(data);
	
	int nPtbin = 6;
	int nMassbin = 20;
	double min_mass = 2.2;
	double max_mass = 4;
		
	int i;
	double min_pt = 2;
	double max_pt = 20;

	int EtaBin = 2;
	std::string etaRegion[EtaBin];
	etaRegion[0]= "BarrelCentral";
	etaRegion[1]= "BarrelEdge";

	
	

	TH2D* Dimass_pt_mc[EtaBin];
	TH2D* Dimass_pt_data[EtaBin];
	//TH2D* Dimass_pt_LBlock = new TH2D("diLepton mass vs lumiblock ","diLepton mass vs luiblock ",nMassbin,min_mass,max_mass,nLBlockbin,min_Lblock,max_Lblock);
	
	for (i=0;i<EtaBin; i++){

	Dimass_pt_mc[i] = new TH2D(("diLepton mass vs lepton pt "+etaRegion[i]).c_str(),("diLepton mass vs lepton pt "+etaRegion[i]).c_str(),nMassbin,min_mass,max_mass,nPtbin,min_pt,max_pt);
	Dimass_pt_data[i] = new TH2D(("diLepton mass vs lepton pt "+etaRegion[i]).c_str(),("diLepton mass vs lepton pt "+etaRegion[i]).c_str(),nMassbin,min_mass,max_mass,nPtbin,min_pt,max_pt);
	
	}	


	std::cout << "------------------------------------" << std::endl;
	std::cout << "Loaded " << sgl_data.fChain->GetEntries() << " entries for data" << std::endl;
	std::cout << "Loaded " << sgl_mc.fChain->GetEntries() << " entries for MC" << std::endl;
	std::cout << "------------------------------------" << std::endl;
	
	for(i=0; i<sgl_mc.fChain->GetEntries();i++){
	
	FillHistosFromB(&sgl_mc,Dimass_pt_mc[0],i,1,TwoTrack,0);	
	FillHistosFromB(&sgl_mc,Dimass_pt_mc[1],i,1,TwoTrack,1);	

	}


	for(i=0; i<sgl_data.fChain->GetEntries();i++){

	FillHistos(&sgl_data,Dimass_pt_data[0],i,0,TwoTrack,0);	
	FillHistos(&sgl_data,Dimass_pt_data[1],i,0,TwoTrack,1);	

	}
	
	SavePlot2D("m_{e^{+}e^{-}}(GeV)","electron_pt (GeV)",Dimass_pt_mc[0],"Dimass_pt_mc",0,0);
	SavePlot2D("m_{e^{+}e^{-}}(GeV)","electron_pt (GeV)",Dimass_pt_data[0],"Dimass_pt_data",0,0);
	std::string combination;
	if (TwoTrack) combination = "Twotrack_";
	else combination = "ECAL-track_";
	DoubleSlicer("/plots/",nPtbin,min_pt, max_pt,"m_{e^{+}e^{-}}",Dimass_pt_mc[0],Dimass_pt_data[0],(combination+"MassSlice_BarrelCentral").c_str(),1,etaRegion[0]);	
	DoubleSlicer("/plots/",nPtbin,min_pt, max_pt,"m_{e^{+}e^{-}}",Dimass_pt_mc[1],Dimass_pt_data[1],(combination+"MassSlice_BarrelEdge").c_str(),1,etaRegion[1]);      
	delete Dimass_pt_mc[0];	
	delete Dimass_pt_mc[1];	
	delete Dimass_pt_data[0];	
	delete Dimass_pt_data[1];	
	return 0;
}




void FillHistos(DileptonGeneralClass* sgl,TH2D* hist, int i,bool isMC, bool Twotrack, bool EtaOuter){

		
		sgl->fChain->GetEntry(i);	
		bool mass; 
		bool barrel;
		bool eta;
		bool ele_ID;

		double weight;
		
		
		if (isMC) weight = 1;
		else weight =1.;

		if (i%10000==0 ) std::cout << "at entry" << i << std::endl;
			
		if(sgl->nelectronPair == 0) return;	
		int j =0;
		barrel = std::max(fabs(sgl->Electron_eta[sgl->electronPair_l1Idx[j]]),fabs(sgl->Electron_eta[sgl->electronPair_l2Idx[j]]))>1.5;
		mass = (sgl->electronPair_mass[j]<5.2 || sgl->electronPair_mass[j]>5.4);	
	
		if (!EtaOuter)	eta =  std::max(fabs(sgl->Electron_eta[sgl->electronPair_l1Idx[j]]),fabs(sgl->Electron_eta[sgl->electronPair_l2Idx[j]]))>1;
		else 	eta = !(std::max(fabs(sgl->Electron_eta[sgl->electronPair_l1Idx[j]]),fabs(sgl->Electron_eta[sgl->electronPair_l2Idx[j]]))>1);
		ele_ID = sgl->Electron_pfmvaId[sgl->electronPair_l1Idx[j]]<0;
		//	 if((sgl->electronPair_mass[j]<5.2 || sgl->electronPair_mass[j]>5.4) ||  fabs(sgl->Electron_eta[sgl->electronPair_l1Idx[j]])> 1 || fabs(sgl->Electron_eta[sgl->electronPair_l2Idx[j]])> 1 || sgl->Electron_pfmvaId[sgl->electronPair_l1Idx[j]]<0  ) return;
			if (mass||barrel||eta||ele_ID) return;
			double pt1_AngEcal,pt1_AngTrack;
			double pt2_AngEcal,pt2_AngTrack;
	
			pt1_AngEcal = sqrt(pow(sgl->electronPair_l1_SC_RegEnergy[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl->electronPair_l1_SCeta[j])));
			pt1_AngTrack = sqrt(pow(sgl->electronPair_l1_SC_RegEnergy[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl->Electron_eta[sgl->electronPair_l1Idx[j]])));
			pt2_AngEcal = sqrt(pow(sgl->electronPair_l2_SC_RegEnergy[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl->electronPair_l2_SCeta[j])));
			pt2_AngTrack = sqrt(pow(sgl->electronPair_l2_SC_RegEnergy[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl->Electron_eta[sgl->electronPair_l2Idx[j]])));
			
			TLorentzVector l1_angTrack, l2_angTrack,l1_angEcal,l2_angEcal,l1,l2,e1,e2;
//			std::cout << "in vector initialization" << std::endl;
			l1_angTrack.SetPtEtaPhiM(pt1_AngTrack,sgl->Electron_eta[sgl->electronPair_l1Idx[j]],sgl->Electron_phi[sgl->electronPair_l1Idx[j]],0.000511);
			l2_angTrack.SetPtEtaPhiM(pt2_AngTrack,sgl->Electron_eta[sgl->electronPair_l2Idx[j]],sgl->Electron_phi[sgl->electronPair_l2Idx[j]],0.000511);
			l1_angEcal.SetPtEtaPhiM(pt1_AngEcal,sgl->electronPair_l1_SCeta[j],sgl->electronPair_l1_SCphi[j],0.000511);
			l2_angEcal.SetPtEtaPhiM(pt2_AngEcal,sgl->electronPair_l2_SCeta[j],sgl->electronPair_l2_SCphi[j],0.000511);
			l1.SetPtEtaPhiM(sgl->Electron_pt[sgl->electronPair_l1Idx[j]],sgl->Electron_eta[sgl->electronPair_l1Idx[j]],sgl->Electron_phi[sgl->electronPair_l1Idx[j]],0.000511);
			l2.SetPtEtaPhiM(sgl->Electron_pt[sgl->electronPair_l2Idx[j]],sgl->Electron_eta[sgl->electronPair_l2Idx[j]],sgl->Electron_phi[sgl->electronPair_l2Idx[j]],0.000511);

			if(Twotrack){
				e1 = l1;
				e2 = l2;
			  hist->Fill((e1+e2).M(),sgl->Electron_pt[sgl->electronPair_l1Idx[j]],weight);
			  hist->Fill((e1+e2).M(),sgl->Electron_pt[sgl->electronPair_l2Idx[j]],weight);
	       		}else{	


		          e1 = l1_angTrack;
			  e2 = l2;
			  hist->Fill((e1+e2).M(),sgl->Electron_pt[sgl->electronPair_l1Idx[j]],weight);
		          e1 = l1;
			  e2 = l2_angTrack;
			  hist->Fill((e1+e2).M(),sgl->Electron_pt[sgl->electronPair_l2Idx[j]],weight);

		}

		
	

}
void FillHistosFromB(DiLeptonFromBClass* sgl,TH2D* hist, int i,bool isMC, bool Twotrack, bool EtaOuter){

		
		sgl->fChain->GetEntry(i);	

		double weight;
		bool mass; 
		bool barrel;
		bool eta;
		bool ele_ID;
		
		
		if (isMC) weight = 1;
		else weight =1.;

		if (i%10000==0 ) std::cout << "at entry" << i << std::endl;
			
		if(sgl->nBToKEE == 0) return;	
		int j =0;
		barrel = std::max(fabs(sgl->Electron_eta[sgl->BToKEE_l1Idx[j]]),fabs(sgl->Electron_eta[sgl->BToKEE_l2Idx[j]]))>1.5;
		mass = (sgl->BToKEE_mass[j]<5.2 || sgl->BToKEE_mass[j]>5.4);	
		if (!EtaOuter) eta =  std::max(fabs(sgl->Electron_eta[sgl->BToKEE_l1Idx[j]]),fabs(sgl->Electron_eta[sgl->BToKEE_l2Idx[j]]))>1;
		else eta = !( std::max(fabs(sgl->Electron_eta[sgl->BToKEE_l1Idx[j]]),fabs(sgl->Electron_eta[sgl->BToKEE_l2Idx[j]]))>1);
		ele_ID = sgl->Electron_pfmvaId[sgl->BToKEE_l1Idx[j]]<0;
			// if((sgl->BToKEE_mass[j]<5.2 || sgl->BToKEE_mass[j]>5.4) ||  fabs(sgl->Electron_eta[sgl->BToKEE_l1Idx[j]])> 1 || fabs(sgl->Electron_eta[sgl->BToKEE_l2Idx[j]])> 1 || sgl->Electron_pfmvaId[sgl->BToKEE_l1Idx[j]]<0 ) return;
			if (mass||barrel||eta||ele_ID) return;
			double pt1_AngEcal,pt1_AngTrack;
			double pt2_AngEcal,pt2_AngTrack;
	
			pt1_AngEcal = sqrt(pow(sgl->BToKEE_l1_SC_RegEnergy[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl->BToKEE_l1_SCeta[j])));
			pt1_AngTrack = sqrt(pow(sgl->BToKEE_l1_SC_RegEnergy[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl->Electron_eta[sgl->BToKEE_l1Idx[j]])));
			pt2_AngEcal = sqrt(pow(sgl->BToKEE_l2_SC_RegEnergy[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl->BToKEE_l2_SCeta[j])));
			pt2_AngTrack = sqrt(pow(sgl->BToKEE_l2_SC_RegEnergy[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl->Electron_eta[sgl->BToKEE_l2Idx[j]])));
			
			TLorentzVector l1_angTrack, l2_angTrack,l1_angEcal,l2_angEcal,l1,l2,e1,e2;
//			std::cout << "in vector initialization" << std::endl;
			l1_angTrack.SetPtEtaPhiM(pt1_AngTrack,sgl->Electron_eta[sgl->BToKEE_l1Idx[j]],sgl->Electron_phi[sgl->BToKEE_l1Idx[j]],0.000511);
			l2_angTrack.SetPtEtaPhiM(pt2_AngTrack,sgl->Electron_eta[sgl->BToKEE_l2Idx[j]],sgl->Electron_phi[sgl->BToKEE_l2Idx[j]],0.000511);
			l1_angEcal.SetPtEtaPhiM(pt1_AngEcal,sgl->BToKEE_l1_SCeta[j],sgl->BToKEE_l1_SCphi[j],0.000511);
			l2_angEcal.SetPtEtaPhiM(pt2_AngEcal,sgl->BToKEE_l2_SCeta[j],sgl->BToKEE_l2_SCphi[j],0.000511);
			l1.SetPtEtaPhiM(sgl->Electron_pt[sgl->BToKEE_l1Idx[j]],sgl->Electron_eta[sgl->BToKEE_l1Idx[j]],sgl->Electron_phi[sgl->BToKEE_l1Idx[j]],0.000511);
			l2.SetPtEtaPhiM(sgl->Electron_pt[sgl->BToKEE_l2Idx[j]],sgl->Electron_eta[sgl->BToKEE_l2Idx[j]],sgl->Electron_phi[sgl->BToKEE_l2Idx[j]],0.000511);

			if(Twotrack){
				e1 = l1;
				e2 = l2;
			  hist->Fill((e1+e2).M(),sgl->Electron_pt[sgl->BToKEE_l1Idx[j]],weight);
			  hist->Fill((e1+e2).M(),sgl->Electron_pt[sgl->BToKEE_l2Idx[j]],weight);
	       		}else{	


		          e1 = l1_angTrack;
			  e2 = l2;
			  hist->Fill((e1+e2).M(),sgl->Electron_pt[sgl->BToKEE_l1Idx[j]],weight);
		          e1 = l1;
			  e2 = l2_angTrack;
			  hist->Fill((e1+e2).M(),sgl->Electron_pt[sgl->BToKEE_l2Idx[j]],weight);

		}

		
	

}
