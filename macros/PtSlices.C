#include "../src/BParkTools.cc"
#include "../src/DiLeptonFromBClass.cc"
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
#include "TColor.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TString.h"
#include "TLorentzVector.h"



void FillHistos(DiLeptonFromBClass *sgl,int i, TH2D* hist, double, double,double,int);


int BinP(){

	//std::ifstream bkg_input("pattern.txt");
	//std::string Notrack_path = "/eos/home-r/ratramon/ECAL_lin/Run2018D_part1.root";
	std::string Onetrack_path = "/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Oct01_withB/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_Toee_ext/201001_141013/0000/*.root";

/*	std::string bkg_path = "/eos/home-r/ratramon/BToKEE/";
	if(!bkg_input) {
		std::cout << "Cannot open input file." << std::endl ;
		return 1;
	}*/

//	TChain *Notrack = new TChain("Events");
	TChain *Onetrack = new TChain("Events");

//	Notrack->Add(Notrack_path.c_str());
	Onetrack->Add(Onetrack_path.c_str());
	
	int n[10];
	setStyle();
	gStyle->SetPalette(91,n);
	//std::cout << "Loaded " << mc->GetEntries() << " entries for signal" << std::endl;
	//std::cout << "Loaded " << data->GetEntries() << " entries for signal" << std::endl;
	//std::cout << "Loaded " << background->GetEntries() << " entries for background" << std::endl;

//	DiLeptonFromBClass sgl_notrack;
	DiLeptonFromBClass sgl_Onetrack;
	DiLeptonFromBClass* sgl_Onetrack_ptr;
//	sgl_notrack.Init(Notrack);
	sgl_Onetrack.Init(Onetrack);
	
	sgl_Onetrack_ptr = &sgl_Onetrack;
	
	int nPbin = 5;	
	int i;
	double min_p = 0;
	double max_p = 100;
	double min_E =0;
	double max_E = 120;
	double BMassMin = 4;
	double BMassMax = 6;
	double pMin = 60;
	std::vector<TH1D*> slices_onetrack;
	std::vector<TH1D*> slices_twotrack;
	std::vector<TH1D*> slices_notrack_angEcal;
	std::vector<TH1D*> slices_notrack_angtrack;
	std::vector<TH1D*> slices_onetrack_angtrack;

	TH2D* DimassP_notrack_AngEcal = new TH2D("diLepton mass vs lepton p-notrack ecal ","diLepton mass vs lepton pt-notrack ecal ",100,0,5,nPbin,min_p,max_p);
	TH2D* DimassP_notrack_Angtrack = new TH2D("diLepton mass vs lepton p-notrack ","diLepton mass vs lepton pt-notrack ",100,0,5,nPbin,min_p,max_p);
	TH2D* DimassP_onetrack = new TH2D("diLepton mass vs lepton pt-1track ","diLepton mass vs lepton pt-1track ",100,0,5,nPbin,min_p,max_p);
	TH2D* DimassP_onetrack_Angtrack = new TH2D("diLepton mass vs lepton pt-1tracki-ecal trk ang","diLepton mass vs lepton pt-1track-ecal trk ang ",100,0,5,nPbin,min_p,max_p);
	TH2D* DimassP_twotrack = new TH2D("diLepton mass vs lepton pt-2tracks ","diLepton mass vs lepton pt-2tracks ",100,0,5,nPbin,min_p,max_p);
	

	std::cout << "Loaded " << sgl_Onetrack.fChain->GetEntries() << " entries for signal" << std::endl;

	for(i=0; i<sgl_Onetrack.fChain->GetEntries();i++){
	
	
	FillHistos(&sgl_Onetrack,i,DimassP_notrack_AngEcal,pMin,BMassMin,BMassMax,3);	
	FillHistos(&sgl_Onetrack,i,DimassP_onetrack,pMin,BMassMin,BMassMax,0);	
	FillHistos(&sgl_Onetrack,i,DimassP_twotrack,pMin,BMassMin,BMassMax,1);	
	FillHistos(&sgl_Onetrack,i,DimassP_notrack_Angtrack,pMin,BMassMin,BMassMax,2);	
	FillHistos(&sgl_Onetrack,i,DimassP_onetrack_Angtrack,pMin,BMassMin,BMassMax,4);

	}	
	
//	SavePlot2D("m_{e^{+}e^{-}}(GeV)","electron_p (GeV)",DimassP_notrack,"DimassP_notrack",0,0);
	//SavePlot2D("m_{e^{+}e^{-}}(GeV)","electron_pt (GeV)",Dimass_pt_data,"Dimass_pt_data",0,0);
	
	Slicer(&slices_notrack_angEcal,"/plots/",nPbin,min_p, max_p,"m_{e^{+}e^{-}}",DimassP_notrack_AngEcal,"mll_PSlices_notrackEcalAng",0);	
	Slicer(&slices_notrack_angtrack,"/plots/",nPbin,min_p, max_p,"m_{e^{+}e^{-}}",DimassP_notrack_Angtrack,"mll_PSlices_notracktrackAng",0);	
	Slicer(&slices_onetrack,"/plots/",nPbin,min_p, max_p,"m_{e^{+}e^{-}}",DimassP_onetrack,"mll_PSlices_1track",0);	
	Slicer(&slices_twotrack,"/plots/",nPbin,min_p, max_p,"m_{e^{+}e^{-}}(GeV)",DimassP_twotrack,"mll_PSlices_2track",0);
	Slicer(&slices_onetrack_angtrack,"/plots/",nPbin,min_p, max_p,"m_{e^{+}e^{-}}",DimassP_onetrack_Angtrack,"mll_PSlices_1tracktrackAng",0);	

	std::cout << slices_twotrack.at(0)->GetEntries()<< std::endl;
	for (i=0;i<nPbin;i++){
	double ymax =  std::max(slices_twotrack.at(i)->GetMaximum()*1.2,slices_onetrack_angtrack.at(i)->GetMaximum()*1.2);
	TH2D* plotter = new TH2D("plot","plot",10,0,5,10,0,ymax);
	plotter->GetXaxis()->SetTitle("m_{e^{+}e^{-}}(GeV)");
	plotter->GetYaxis()->SetTitle("entries");
	TCanvas* canva = new TCanvas(("Prange_"+std::to_string(i)).c_str(),("Prange_"+std::to_string(i)).c_str(),800,600);
	slices_twotrack.at(i)->SetLineColor(kTeal+9);
	slices_onetrack.at(i)->SetLineColor(kBlue-3);
	slices_onetrack_angtrack.at(i)->SetLineColor(kOrange-3);
	slices_notrack_angEcal.at(i)->SetLineColor(kRed-3);
	slices_notrack_angtrack.at(i)->SetLineColor(kMagenta-5);
	plotter->Draw("");
	slices_twotrack.at(i)->Draw("same");
	slices_onetrack.at(i)->Draw("same");
	slices_onetrack_angtrack.at(i)->Draw("same");
	slices_notrack_angEcal.at(i)->Draw("same");
	slices_notrack_angtrack.at(i)->Draw("same");

	TLatex* l = new TLatex();
	l->SetTextAlign(13);
	l->SetTextSize(0.03);
	l->SetTextColor(slices_twotrack.at(i)->GetLineColor());
	l->DrawLatex(0.5,ymax*0.8,"track-track");
	l->SetTextColor(slices_onetrack.at(i)->GetLineColor());
	l->DrawLatex(0.5,ymax*0.75,"track-ECAL w/ ECAL #eta, #phi");
	l->SetTextColor(slices_onetrack_angtrack.at(i)->GetLineColor());
	l->DrawLatex(0.5,ymax*0.7,"ECAL-track w/ track #eta, #phi");
	l->SetTextColor(slices_notrack_angEcal.at(i)->GetLineColor());
	l->DrawLatex(0.5,ymax*0.65,"ECAL-ECAL w/ ECAL #eta, #phi");
	l->SetTextColor(slices_notrack_angtrack.at(i)->GetLineColor());
	l->DrawLatex(0.5,ymax*0.6,"ECAL-ECAL w/ track #eta, #phi");
		
	canva->SaveAs(("mEE_prange"+std::to_string(i)+".png").c_str());	
	}	
//	Slicer("/plots/",4,min_p, max_p,"m_{e^{+}e^{-}}",DimassP_1track,"mll_PSlices_1track",0);	
//	Slicer("/plots/",4,min_p, max_p,"m_{e^{+}e^{-}}",DimassP_2track,"mll_PSlices_2tracks",0);	
	delete DimassP_notrack_AngEcal;	
	delete DimassP_notrack_Angtrack;	
	delete DimassP_onetrack;	
	delete DimassP_twotrack;	
	return 0;
}




void FillHistos(DiLeptonFromBClass *sgl,int i, TH2D* hist,double pMin,double BmassMin,double BmassMax, int TwoTrack){

//	int i;
//	TH1D* h_diff = new TH1D("energy diff","energy diff",50,-20,20);
	
		sgl->fChain->GetEntry(i);
		if (i%10000==0)std::cout << "at entry" << i << std::endl;
	//	std::cout << "at entry " << i <<"mass " <<  sgl->nBToKEE << " " << std::endl;
		if(sgl->nBToKEE == 0) return;	
		int j=0;
	//	for (int j=0; j<1; j++){	

			 if((sgl->BToKEE_mass[j]< 4 || sgl->BToKEE_mass[j]>6) ||  fabs(sgl->Electron_eta[sgl->BToKEE_l1Idx[j]])> 1.5 || fabs(sgl->Electron_eta[sgl->BToKEE_l2Idx[j]])> 1.5) return;

			double p_lead, p_sublead;
//			std::cout << "in b selection" << std::endl;
			 p_lead = fabs(sgl->Electron_pt[sgl->BToKEE_l1Idx[j]]/sin(2*atan(exp(-sgl->Electron_eta[sgl->BToKEE_l1Idx[j]]))));
			 p_sublead = fabs(sgl->Electron_pt[sgl->BToKEE_l2Idx[j]]/sin(2*atan(exp(-sgl->Electron_eta[sgl->BToKEE_l2Idx[j]]))));


			if(p_lead>p_sublead){ 
			
			TLorentzVector l1_angTrack, l2_angTrack,l1_angEcal,l2_angEcal, l1,l2;
			
			double pt1_AngEcal,pt1_AngTrack;
			double pt2_AngEcal,pt2_AngTrack;
	
			pt1_AngEcal = sqrt(pow(sgl->BToKEE_l1_SC_RegEnergy[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl->BToKEE_l1_SCeta[j])));
			pt1_AngTrack = sqrt(pow(sgl->BToKEE_l1_SC_RegEnergy[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl->Electron_eta[sgl->BToKEE_l1Idx[j]])));
			pt2_AngEcal = sqrt(pow(sgl->BToKEE_l2_SC_RegEnergy[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl->BToKEE_l2_SCeta[j])));
			pt2_AngTrack = sqrt(pow(sgl->BToKEE_l2_SC_RegEnergy[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl->Electron_eta[sgl->BToKEE_l2Idx[j]])));
			
//			std::cout << "in vector initialization" << std::endl;
			l1_angTrack.SetPtEtaPhiM(pt1_AngTrack,sgl->Electron_eta[sgl->BToKEE_l1Idx[j]],sgl->Electron_phi[sgl->BToKEE_l1Idx[j]],0.000511);
			l2_angTrack.SetPtEtaPhiM(pt2_AngTrack,sgl->Electron_eta[sgl->BToKEE_l2Idx[j]],sgl->Electron_phi[sgl->BToKEE_l2Idx[j]],0.000511);
			l1_angEcal.SetPtEtaPhiM(pt1_AngEcal,sgl->BToKEE_l1_SCeta[j],sgl->BToKEE_l1_SCphi[j],0.000511);
			l2_angEcal.SetPtEtaPhiM(pt2_AngEcal,sgl->BToKEE_l2_SCeta[j],sgl->BToKEE_l2_SCphi[j],0.000511);
			l1.SetPtEtaPhiM(sgl->Electron_pt[sgl->BToKEE_l1Idx[j]],sgl->Electron_eta[sgl->BToKEE_l1Idx[j]],sgl->Electron_phi[sgl->BToKEE_l1Idx[j]],0.000511);
			l2.SetPtEtaPhiM(sgl->Electron_pt[sgl->BToKEE_l2Idx[j]],sgl->Electron_eta[sgl->BToKEE_l2Idx[j]],sgl->Electron_phi[sgl->BToKEE_l2Idx[j]],0.000511);
			

//		std::cout << " " <<  std::endl;
//	std::cout << "check p E lead " << p_lead << " " << sgl->Electron_energySC[sgl->BToKEE_l1Idx[j]] << " " << sgl->Electron_eta[sgl->BToKEE_l1Idx[j]] << " "<< sgl->Electron_phi[sgl->BToKEE_l1Idx[j]] <<  std::endl;
//	std::cout << "check p E sublead" << p_sublead << " " << sgl->Electron_energySC[sgl->BToKEE_l2Idx[j]] <<  " " << sgl->Electron_eta[sgl->BToKEE_l2Idx[j]] << " "<< sgl->Electron_phi[sgl->BToKEE_l2Idx[j]] << std::endl;
//	std::cout << "check masses " << sgl->BToKEE_mll_raw[j] << " " << sgl->BToKEE_mll_SCraw[j] << " " << sgl->BToKEE_mll_SC1track[j]<<  std::endl;
//		std::cout << " " <<  std::endl;
	


				if (TwoTrack==1)	 hist->Fill(sgl->BToKEE_mll_llfit[j],sqrt(pow(sgl->Electron_energySC[sgl->BToKEE_l2Idx[j]],2)-0.000511*0.000511));
				else if (TwoTrack==2)    hist->Fill((l1_angTrack+l2_angTrack).M(),sqrt(pow(sgl->Electron_energySC[sgl->BToKEE_l2Idx[j]],2)-0.000511*0.000511));
				else if (TwoTrack==3)    hist->Fill((l1_angEcal+l2_angEcal).M(),sqrt(pow(sgl->Electron_energySC[sgl->BToKEE_l2Idx[j]],2)-0.000511*0.000511));
                                else if (TwoTrack==0)	 hist->Fill((l1_angEcal+l2).M(),sqrt(pow(sgl->Electron_energySC[sgl->BToKEE_l2Idx[j]],2)-0.000511*0.000511));
                                else if (TwoTrack==4)	 hist->Fill((l1_angTrack+l2).M(),sqrt(pow(sgl->Electron_energySC[sgl->BToKEE_l2Idx[j]],2)-0.000511*0.000511));
		//	break;			

		//	}
		}

//	SavePlot("E_{SC}-E_{trk}",h_diff,"energyDiff",false,NULL,false);
	}







