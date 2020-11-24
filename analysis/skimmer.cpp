//#include "../interface/BParkTools.h"
#include "../interface/MeeClass.h"
#include "../interface/DileptonClass.h"
#include "../interface/DiLeptonClass.h"
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


TH1D* ComputePUweight(bool , bool , std::string);
double PUweight(TH1D*, int);

int main(int argc, char** argv ){


	int i,j;
	bool isMC = true;
	std::string letter = argv[1];
	int  part = atoi(argv[2]);
	bool isZ = atoi(argv[4]);
	std::string IDir = argv[5];
	int subdir = atoi(argv[3]);
	std::string datapath,datapath1,datapath2;
	TChain *Onetrack = new TChain("Events");
	if(isMC){
		if (isZ){
			datapath ="/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Nov17_withB/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_AmcatNLO_EGamma/201117_095235/0000/BParkNANO_data_2020Nov17_withB_*.root";	
			//datapath ="/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Nov06_withB/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_AmcatNLO_EGamma/201106_142211/0000/BParkNANO_data_2020Nov06_withB_*.root";	
		//	datapath1 = "/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Nov02_withB/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_AmcatNLO_EGamma/201102_155913/0001/*.root";	
		//	datapath2 = "/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Nov02_withB/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_AmcatNLO_EGamma/201102_155913/0002/*.root";	
//	Onetrack->Add(datapath1.c_str());
//	Onetrack->Add(datapath2.c_str());
		}else{

			datapath= "/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Oct22_withB/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_Toee_ext/201022_132632/0000/*.root";

	}
	std::cout << datapath << std::endl;
	Onetrack->Add(datapath.c_str());
	}else{
		if(isZ){
			datapath = "/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Nov06_withB/EGamma/crab_data_EGamma/201106_095929/000"+std::to_string(subdir)+"/BParkNANO_data_2020Nov06_withB_*.root";
		}else{
			datapath= "/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Nov04_withB/ParkingBPH"+std::to_string(part)+"/crab_data_Run2018"+letter+"_part"+std::to_string(part)+"/"+IDir+"/000"+std::to_string(subdir)+"/BParkNANO_data_2020Nov04_withB_*.root";
			}
	Onetrack->Add(datapath.c_str());
	std::cout << datapath << std::endl;
	}
	DileptonClass sgl;
	sgl.Init(Onetrack);
	std::string outpath = "/afs/cern.ch/work/r/ratramon/ECAL_linearity/ECALinAnalysis/ECALlinearity/data/";	
	
	//std::cout << "HEREEE "  << std::endl;
	int nPtbin, EtaBin;
	double min_pt, max_pt;
	std::string resonance;

	if (isZ){

	nPtbin = 10;
	min_pt = 10;
	max_pt = 110;
	resonance = "Z_smear";	
	


	}else{

	nPtbin = 18;
	min_pt = 2;
	max_pt = 20;
	resonance = "JPsi";	

	}

	EtaBin = 2;
	bool eta[EtaBin];
		
	std::string etaRegion[EtaBin];
	etaRegion[0]= "BarrelCentral";
	etaRegion[1]= "BarrelEdge";

	double ptMin[nPtbin];
	double ptMax[nPtbin];
	double min_edge[nPtbin]= {10,40,50,75};
	double max_edge[nPtbin]= {40,50,75,150};

	bool mass; 
	bool barrel;
	bool ele_ID;
	for(i=0;i<nPtbin;i++){
	/*	if (isZ){
	 
			ptMin[i] = min_edge[i];
			ptMax[i] = max_edge[i];
		}else{*/
			ptMin[i] = min_pt+ (max_pt-min_pt)/nPtbin*(i);
			ptMax[i] = min_pt+ (max_pt-min_pt)/nPtbin*(i+1);

	//	}
	}
	TFile* OutFile[2][nPtbin];
	TTree* newTree[2][nPtbin];
	MeeClass newClass[2][nPtbin];

	for(i=0;i<2;i++){
	    for(j=0;j<nPtbin;j++){
		
		if(!isMC)OutFile[i][j] = new TFile((outpath+resonance+std::to_string((int)ptMin[j])+"_"+std::to_string((int)ptMax[j])+"_"+etaRegion[i]+"/mllTree_PtMin"+std::to_string((int)ptMin[j])+"PtMax"+std::to_string((int)ptMax[j])+etaRegion[i]+"_"+std::to_string(part)+"_"+letter+"_"+std::to_string(subdir)+".root").c_str(),"RECREATE");
	//	if(!isMC && isZ)OutFile[i][j] = new TFile((outpath+resonance+std::to_string((int)ptMin[j])+"_"+std::to_string((int)ptMax[j])+"/mllTree_PtMin"+std::to_string((int)ptMin[j])+"PtMax"+std::to_string((int)ptMax[j])+"_"+std::to_string(part)+"_"+letter+"_"+std::to_string(subdir)+".root").c_str(),"RECREATE");
		else if (isMC) OutFile[i][j] = new TFile((outpath+"MC"+resonance+"/mllTreeMC_PtMin"+std::to_string((int)ptMin[j])+"PtMax"+std::to_string((int)ptMax[j])+etaRegion[i]+".root").c_str(),"RECREATE");
		if(isMC)newTree[i][j] = new TTree("treeMC",(std::to_string(i+j)).c_str());
		else newTree[i][j] = new TTree("tree",(std::to_string(i+j)).c_str());
		newClass[i][j].InitTree(newTree[i][j]);
	std::cout << "outfile " << OutFile[i][j]->GetName() << std::endl;
		}
	}
	std::cout << "loading " << sgl.fChain->GetEntries() << std::endl;
	TH1D* h_PUweight;
	if(isMC)h_PUweight = ComputePUweight(1,isZ,IDir);
	for(i=0; i<sgl.fChain->GetEntries();i++){
			if (i%1000000==0) std::cout << "on entry " << i << std::endl;
			sgl.fChain->GetEntry(i);
			if (sgl.nelectronPair==0)continue;
			int j =0;
						
			barrel = std::min(fabs(sgl.Electron_eta[sgl.electronPair_l1Idx[j]]),fabs(sgl.Electron_eta[sgl.electronPair_l2Idx[j]]))>1.5;
	//		if (isZ){
			mass = false;
			ele_ID = false;
	//		}
	//		else{ mass = (sgl.electronPair_mass[j]<4 || sgl.electronPair_mass[j]>6);	
	//		ele_ID = sgl.Electron_pfmvaId[sgl.electronPair_l1Idx[j]]<0;

	//		}
			eta[0]=  std::max(fabs(sgl.Electron_eta[sgl.electronPair_l1Idx[j]]),fabs(sgl.Electron_eta[sgl.electronPair_l2Idx[j]]))<1;
			eta[1] = !eta[0];
			if (isZ && barrel ) continue;
			if (!isZ && (mass||barrel||ele_ID)) continue;
			double pt1_AngTrack;
			double pt2_AngTrack;
	
//			pt1_AngTrack = sqrt(pow(sgl.electronPair_l1_SC_RegEnergy[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl.Electron_eta[sgl.electronPair_l1Idx[j]])));
//			pt2_AngTrack = sqrt(pow(sgl.electronPair_l2_SC_RegEnergy[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl.Electron_eta[sgl.electronPair_l2Idx[j]])));
			pt1_AngTrack = sqrt(pow(sgl.electronPair_l1_SmearE[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl.Electron_eta[sgl.electronPair_l1Idx[j]])));
			pt2_AngTrack = sqrt(pow(sgl.electronPair_l2_SmearE[j]+0.000511,2)-pow(0.000511,2))*sin(2*atan(exp(-sgl.Electron_eta[sgl.electronPair_l2Idx[j]])));
			
			TLorentzVector l1_angTrack, l2_angTrack,l1_angEcal,l2_angEcal,l1,l2,e1,e2;
			l1_angTrack.SetPtEtaPhiM(pt1_AngTrack,sgl.Electron_eta[sgl.electronPair_l1Idx[j]],sgl.Electron_phi[sgl.electronPair_l1Idx[j]],0.000511);
			l2_angTrack.SetPtEtaPhiM(pt2_AngTrack,sgl.Electron_eta[sgl.electronPair_l2Idx[j]],sgl.Electron_phi[sgl.electronPair_l2Idx[j]],0.000511);
			l1.SetPtEtaPhiM(sgl.Electron_pt[sgl.electronPair_l1Idx[j]],sgl.Electron_eta[sgl.electronPair_l1Idx[j]],sgl.Electron_phi[sgl.electronPair_l1Idx[j]],0.000511);
			l2.SetPtEtaPhiM(sgl.Electron_pt[sgl.electronPair_l2Idx[j]],sgl.Electron_eta[sgl.electronPair_l2Idx[j]],sgl.Electron_phi[sgl.electronPair_l2Idx[j]],0.000511);
			
			for (int k=0;k<nPtbin;k++){					
				for(int h=0;h<2;h++){
					if(eta[h]){
					if(isZ){

					if ( (sgl.Electron_pt[sgl.electronPair_l1Idx[j]]> ptMin[k] && sgl.Electron_pt[sgl.electronPair_l1Idx[j]]< ptMax[k]) && (sgl.Electron_pt[sgl.electronPair_l2Idx[j]]> ptMin[k] && sgl.Electron_pt[sgl.electronPair_l2Idx[j]]< ptMax[k])){
					newClass[h][k].run = sgl.run;
					newClass[h][k].luminosityBlock = sgl.luminosityBlock;
					newClass[h][k].event = sgl.event;
				//	newClass[h][k].B_mass = sgl.electronPair_mass[j];
					newClass[h][k].B_mass = -1;
					newClass[h][k].mll_TT = (l1+l2).M();
					newClass[h][k].mll_ET = (l1_angTrack+l2).M();
					newClass[h][k].mll_EE = (l1_angTrack+l2_angTrack).M();
					if(isMC)newClass[h][k].PUweight = PUweight(h_PUweight,sgl.PV_npvs);
					else newClass[h][k].PUweight = -1;
					//newClass[h][k].fChain->Fill();
					newTree[h][k]->Fill();
					}
					}else{

					if ( (sgl.Electron_pt[sgl.electronPair_l1Idx[j]]> ptMin[k] && sgl.Electron_pt[sgl.electronPair_l1Idx[j]]< ptMax[k])){
	//				std::cout << "check 1" <<std::endl;
					newClass[h][k].run = sgl.run;
					newClass[h][k].luminosityBlock = sgl.luminosityBlock;
					newClass[h][k].event = sgl.event;
				//	newClass[h][k].B_mass = sgl.electronPair_mass[j];
			//		newClass[h][k].B_mass = sgl.electronPair_mass[j];
					newClass[h][k].mll_TT = (l1+l2).M();
					newClass[h][k].mll_ET = (l1_angTrack+l2).M();
					newClass[h][k].mll_EE = (l1_angTrack+l2_angTrack).M();
					if(isMC)newClass[h][k].PUweight = PUweight(h_PUweight,sgl.PV_npvs);
					else newClass[h][k].PUweight = -1;
					//newClass[h][k].fChain->Fill();
					newTree[h][k]->Fill();
					}
					if ( sgl.Electron_pt[sgl.electronPair_l2Idx[j]]> ptMin[k] && sgl.Electron_pt[sgl.electronPair_l2Idx[j]]< ptMax[k]){
	//				std::cout << "check 2" << std::endl;
					newClass[h][k].run = sgl.run;
					newClass[h][k].luminosityBlock = sgl.luminosityBlock;
					newClass[h][k].event = sgl.event;
			//		newClass[h][k].B_mass = sgl.electronPair_mass[j];
					newClass[h][k].mll_TT = (l1+l2).M();
					newClass[h][k].mll_ET = (l1+l2_angTrack).M();
					newClass[h][k].mll_EE = (l1_angTrack+l2_angTrack).M();
					if (isMC)newClass[h][k].PUweight = PUweight(h_PUweight,sgl.PV_npvs);
					else newClass[h][k].PUweight = -1;
					newTree[h][k]->Fill();
					}
					}
			}	
			}
			}
				
				}

 			

	
	for(i=0;i<2;i++){
	    for(j=0;j<nPtbin;j++){
	   // newTree[i][j]->SetDirectory(OutFile[i][j]);  
	    OutFile[i][j]->cd();
	    newTree[i][j]->Write();  
            OutFile[i][j]->Close();
           // delete newTree[i][j];
          //  delete OutFile[i][j];
           // OutFile[i][j]->Close();
		}
	}
	
	return 0;
 }




TH1D* ComputePUweight(bool prompt, bool isZ, std::string IDir){

	std::string mcpath, datapath;	
	TChain* data= new TChain("Events");
	TChain* mc= new TChain("Events");
		if (isZ){
			mcpath ="/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Nov17_withB/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_AmcatNLO_EGamma/201117_095235/0000/BParkNANO_data_2020Nov17_withB_*.root";	
			datapath = "/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Nov06_withB/EGamma/crab_data_EGamma/201106_095929/0000/BParkNANO_data_2020Nov06_withB_*.root";
		//	datapath1 = "/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Nov02_withB/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_AmcatNLO_EGamma/201102_155913/0001/*.root";	
		//	datapath2 = "/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Nov02_withB/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_AmcatNLO_EGamma/201102_155913/0002/*.root";	
//	Onetrack->Add(datapath1.c_str());
//	Onetrack->Add(datapath2.c_str());
		}else{

			mcpath= "/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Oct22_withB/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_BuToKJpsi_Toee_ext/201022_132632/0000/*.root";
			datapath= "/eos/cms/store/group/phys_bphys/ratramon/ECAL_linNano/2020Oct22_withB/ParkingBPH1/crab_data_Run2018D_part1/"+IDir+"/0000/BParkNANO_data_2020Oct22_withB_*.root";

	}
	data->Add(datapath.c_str());
	mc->Add(mcpath.c_str());
	DileptonClass sgl_mc;
	DiLeptonClass sgl_data;
	sgl_mc.Init(mc);
	sgl_data.Init(data);

	//std::cout << "mc entries " <<  mc->GetEntries()  << std::endl;	
//	std::cout << "data entries " << data->GetEntries()  << std::endl;	
	TH1D* h_data =new TH1D ("hdata","hdata",50,0,100);
	TH1D* h_mc = new TH1D("hmc","hmc",50,0,100);
	TH1D* h_mcWeighted = new TH1D("hmcWeigh","hmcWeigh",50,0,100);
	for(int i=0; i<sgl_mc.fChain->GetEntries();i++){
	sgl_mc.fChain->GetEntry(i);
	if(i%1000000==0)std::cout << "on mc entry  " <<  i  << std::endl;	
	h_mc->Fill(sgl_mc.PV_npvs);
	}
	for(int i=0; i<sgl_data.fChain->GetEntries();i++){
	sgl_data.fChain->GetEntry(i);
	if(i%1000000==0)std::cout << "on data entry  " <<  i  << std::endl;	
	
//	std::cout << "data n vertices" << i << " " <<  sgl_data.PV_npvs   << std::endl;	
	h_data->Fill(sgl_data.PV_npvs);
	}
	
//	std::cout << "mc entries " <<  h_mc->GetBinContent(3)  << std::endl;	
//	std::cout << "data entries " <<  h_data->GetBinContent(3)  << std::endl;	
//
	h_data->SaveAs("h_data.root");
	h_mc->Scale(1/h_mc->GetEntries());	
	h_data->Scale(1/h_data->GetEntries());	
	//h_mc->SaveAs("h_data.root");
//	std::cout << "mc entries " <<  h_mc->GetBinContent(15)  << std::endl;	
//	std::cout << "data entries " <<  h_data->GetBinContent(15)  << std::endl;	
	
	h_mc->Divide(h_data);
	
	h_mc->SaveAs("h_mc.root");
	for(int i=0; i<sgl_mc.fChain->GetEntries();i++){
	sgl_mc.fChain->GetEntry(i);
	if(i%1000000==0)std::cout << "on mc entry  " <<  i  << std::endl;
	//std::cout << sgl_mc.PV_npvs << " " << 1/h_mc->GetBinContent(h_mc->GetXaxis()->FindBin(sgl_mc.PV_npvs)) << " " << h_mc->GetXaxis()->FindBin(sgl_mc.PV_npvs) << std::endl;
	if(h_mc->GetBinContent(h_mc->GetXaxis()->FindBin(sgl_mc.PV_npvs))!=0)h_mcWeighted->Fill(sgl_mc.PV_npvs,1/(h_mc->GetBinContent(h_mc->GetXaxis()->FindBin(sgl_mc.PV_npvs))));
	}
	h_mcWeighted->SaveAs("h_weighted.root");
	//std::cout << "mc entries " <<  h_mc->GetBinContent(15)  << std::endl;	
	//std::cout << "mc entries " <<  h_mc->GetBinContent(20)  << std::endl;	
	//std::cout << "mc entries " <<  h_mc->GetBinContent(30)  << std::endl;	
	//std::cout << "mc entries " <<  h_mc->GetBinContent(45)  << std::endl;	
	return h_mc;

	}	


double PUweight(TH1D* h_weight,int n_vertices){
	double weight=-1;
	for (int i=0; i< 14;i++){

	int min_edge = i;
	int max_edge = i+1;

//	std::cout << h_weight->GetBinContent(i+1) << " "<< i  <<  " " << n_vertices  << std::endl;	
	if (n_vertices == max_edge){
	weight = h_weight->GetBinContent(i+1);
	//std::cout << weight << " " << n_vertices  << std::endl;	
	}

	}
	
	if (weight != 0) return 1/weight;
	else return -1;
}
