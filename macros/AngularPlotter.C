#include <iostream>
#include <utility>
#include <vector>
#include <math.h>
#include <string>
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TAxis.h"
#include "RooRealVar.h"
#include "RooBernstein.h"
#include "RooAbsPdf.h"
#include "RooBinning.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooDataSet.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include <time.h>
#include "TBox.h"
//#include "TASImage.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "../src/MeeClass.cc"
#include "../src/BParkTools.cc"




int AngularPlotter(std::string filename){


	setStyle();
	TH1D* DPhi[2],*DEta[2],*Dr[2];	
	
	int nbin;
	double phi_min, phi_max,eta_min,eta_max,dr_min,dr_max;

	nbin = 100;
	phi_min = -4;
	phi_max = 4;
	eta_min = -2;
	eta_max = 2;
	dr_min = 0;
	dr_max = 4;
	

	TChain * data = new TChain("treeMC");
	data->Add((filename).c_str());

	MeeClass sgl;
	sgl.InitTree(data);

	for (int i =0; i<2 ; i++){

	DPhi[i]=new TH1D(("phi"+std::to_string(i)).c_str(),("phi"+std::to_string(i)).c_str(),nbin,phi_min, phi_max);
	DEta[i]=new TH1D(("eta"+std::to_string(i)).c_str(),("eta"+std::to_string(i)).c_str(),nbin,eta_min, eta_max);
	Dr[i]=new TH1D(("dr"+std::to_string(i)).c_str(),("dr"+std::to_string(i)).c_str(),nbin,dr_min, dr_max);

	}

	
	for (int i =0; i<sgl.fChain->GetEntries() ; i++){

	sgl.fChain->GetEntry(i);

	//std::cout << "on entry " << sgl.mll_EE << " phi " << sgl.deltaPhi << " eta " << sgl.deltaEta << " dr " << sgl.deltaR << std::endl;


	if(sgl.mll_EE<100){

	DPhi[0]->Fill(sgl.deltaPhi);
	DEta[0]->Fill(sgl.deltaEta);
	Dr[0]->Fill(sgl.deltaR);



	}else{


	DPhi[1]->Fill(sgl.deltaPhi);
	DEta[1]->Fill(sgl.deltaEta);
	Dr[1]->Fill(sgl.deltaR);


	}



	}

	superpos("#Delta#phi",DPhi[0],DPhi[1],"DeltaPhi_PeakShoulder",0,0);
	superpos("#Delta#eta",DEta[0],DEta[1],"DeltaEta_PeakShoulder",0,0);
	superpos("#DeltaR",Dr[0],Dr[1],"DeltaR_PeakShoulder",0,0);

	return 1;
	



}
