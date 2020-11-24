#include <iostream>
#include <utility>
#include <vector>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
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

#define MASS_JPSI 3.0969
#define BR_NORES  5.50E-07
#define BR_RES  0.0000618

TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold
float lumiTextFont   = 42;  // default is helvetic-bold

bool writeExtraText = true;
TString extraText   = "Preliminary";
TString MCText   = "Simulation";
float extraTextFont = 52;  // default is helvetica-italics

// temll_ETt sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.6;
float lumiTextOffset   = 0.2;
float cmsTextSize      = 0.75;
float cmsTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and emll_ETtra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13TeV = "20.1 fb^{-1}";
TString lumi_8TeV  = "19.7 fb^{-1}";
TString lumi_7TeV  = "5.1 fb^{-1}";
TString lumi_sqrtS = "";

bool drawLogo      = false;

double SignalWeight(){
	double lumi,n_events;
	lumi = 5.187706291*1e15;
	n_events= 45.16919828*1e15;
	std::string MCpath;

	double crossec_mc;
	crossec_mc = 543100000*1e-12; 

	double weight = lumi/n_events;

	return weight;


}

void CMS_lumi( TPad* pad, int iPeriod, int iPosX, double lumi )
{            
	bool outOfFrame    = false;
	if( iPosX/10==0 ) 
	{
		outOfFrame = true;
	}
	int alignY_=3;
	int alignX_=2;
	if( iPosX/10==0 ) alignX_=1;
	if( iPosX==0    ) alignX_=1;
	if( iPosX==0    ) alignY_=1;
	if( iPosX/10==1 ) alignX_=1;
	if( iPosX/10==2 ) alignX_=2;
	if( iPosX/10==3 ) alignX_=3;
	//if( iPosX == 0  ) relPosX = 0.12;
	int align_ = 10*alignX_ + alignY_;

	float H = pad->GetWh();
	float W = pad->GetWw();
	float l = pad->GetLeftMargin();
	float t = pad->GetTopMargin();
	float r = pad->GetRightMargin();
	float b = pad->GetBottomMargin();
	//  float e = 0.025;

	pad->cd();

	TString lumiText;
	if( iPeriod==1 )
	{
		lumiText += lumi_7TeV;
		lumiText += " (7 TeV)";
	}
	else if ( iPeriod==2 )
	{
		lumiText += lumi_8TeV;
		lumiText += " (8 TeV)";
	}
	else if( iPeriod==3 ) 
	{
		lumiText = lumi_8TeV; 
		lumiText += " (8 TeV)";
		lumiText += " + ";
		lumiText += lumi_7TeV;
		lumiText += " (7 TeV)";
	}
	else if ( iPeriod==4 )
	{
		lumiText += lumi_13TeV;
		lumiText += " (13 TeV)";
	}
	else if ( iPeriod==7 )
	{ 
		if( outOfFrame ) lumiText += "#scale[0.85]{";
		lumiText += lumi_13TeV; 
		lumiText += " (13 TeV)";
		lumiText += " + ";
		lumiText += lumi_8TeV; 
		lumiText += " (8 TeV)";
		lumiText += " + ";
		lumiText += lumi_7TeV;
		lumiText += " (7 TeV)";
		if( outOfFrame) lumiText += "}";
	}
	else if ( iPeriod==12 )
	{
		lumiText += "8 TeV";
	}
	else if ( iPeriod==0 )
	{
		lumiText += lumi_sqrtS;
	}
	else if ( iPeriod==5 )
	{
		lumiText += lumi;
		lumiText +=  " fb^{-1}";
	}

	std::cout << lumiText << std::endl;

	TLatex latex;
	latex.SetNDC();
	latex.SetTextAngle(0);
	latex.SetTextColor(kBlack);    

	float extraTextSize = extraOverCmsTextSize*cmsTextSize;

	latex.SetTextFont(42);
	latex.SetTextAlign(31); 
	latex.SetTextSize(lumiTextSize*t);    
//	latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

	outOfFrame = false;
	if( outOfFrame )
	{
		latex.SetTextFont(cmsTextFont);
		latex.SetTextAlign(11); 
		latex.SetTextSize(cmsTextSize*t);    
		latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);
	}

	pad->cd();

	float posX_=0;
	if( iPosX%10<=1 )
	{
		posX_ =   l + relPosX*(1-l-r);
	}
	else if( iPosX%10==2 )
	{
		posX_ =  l + 0.5*(1-l-r);
	}
	else if( iPosX%10==3 )
	{
		posX_ =  1-r - relPosX*(1-l-r);
	}
	float posY_ = 1-t - relPosY*(1-t-b);
	if( !outOfFrame )
	{
		posX_ = r+6*r ;
		posY_ = 1-t - 0.02*(1-t-b);
		if( drawLogo )
		{
		//	posX_ =   r - 0.1*(1-l-r)*W/H;
		//	posY_ = 1-t - 0.045*(1-t-b);
			float xl_0 = posX_;
			float yl_0 = posY_ - 0.15;
			float xl_1 = posX_ + 0.15*H/W;
			float yl_1 = posY_;
			//	TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
			TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
			pad_logo->Draw();
			pad_logo->cd();
			//		CMS_logo->Draw("X");
			pad_logo->Modified();
			pad->cd();
		}
		else
		{
			latex.SetTextFont(cmsTextFont);
			latex.SetTextSize(cmsTextSize*t);
			latex.SetTextAlign(align_);
			latex.DrawLatex(posX_, posY_, cmsText);
			latex.DrawLatex(posX_, posY_-posY_*0.49, cmsText);
			if( writeExtraText ) 
			{
				latex.SetTextFont(extraTextFont);
				latex.SetTextAlign(align_);
				latex.SetTextSize(extraTextSize*t);
			//	latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, cmsText);
				latex.DrawLatex(posX_, posY_-0.03*posY_, MCText);
				latex.DrawLatex(posX_, posY_-posY_*0.52, extraText);
				latex.SetTextFont(lumiTextFont);
				latex.DrawLatex(posX_*0.95, posY_-posY_*0.55,lumiText+" Bparking");
	outOfFrame = false;

			}
		}
	}
	else if( writeExtraText )
	{
		if( iPosX==0) 
		{
			posX_ =   l +  relPosX*(1-l-r);
			posY_ =   1-t+lumiTextOffset*t;
		}
		latex.SetTextFont(extraTextFont);
		latex.SetTextSize(extraTextSize*t);
		latex.SetTextAlign(align_);
		latex.DrawLatex(posX_, posY_, extraText);      
	}
	return;
}



void fit(std::string filename,std::string lable, int  sigPDF,int bkgPDF, bool isZ, bool full,double * resu){

	std::cout << "in fit________________________ " << std::endl;
	RooWorkspace wspace("w");
	gStyle->SetOptFit(0000);
	gROOT->SetBatch(true);
	gROOT->SetStyle("Plain");
	gStyle->SetGridStyle(3);
	gStyle->SetOptStat(000000);
	gStyle->SetOptTitle(0);

	bool blind=true ;
	double eff_res, eff_nores,eff_misreco;
	double xmin, xmax;


	
	//mll_ETmin = hist->GetXaxis()->GetXmin()/*+2*hist->GetXaxis()->GetBinWidth(1)*/;
//	mll_ETmax = hist->GetXaxis()->GetXmax();//-2*hist->GetXaxis()->GetBinWidth(1);
//	std::cout << "bin width " << hist->GetXamll_ETis()->GetBinWidth(1) << std::endl;
//	TH1D* mB_res = new TH1D("sig_norm","sig_norm",(int)((mll_ETmax-xmin)/hist->GetXaxis()->GetBinWidth(1)), xmin, xmax); 
//	TH1D* mB_nores = new TH1D("sig_norm","sig_norm",(int)((mll_ETmax-xmin)/hist->GetXaxis()->GetBinWidth(1)), xmin, xmax); 

//	std::cout << "limits" << mll_ETmin << "   " << xmax << std::endl;

	//wspace.factory(("mll_ET[5.0,"+std::to_string(xmin)+","+std::to_string(xmax)+"]").c_str());
	TFile* DataIn, *DataIn_prompt,*DataIn_Noprompt;
	TTree* data,*data_prompt,*data_noprompt;
	TFile* mcIn;
	if(isZ)
	{
	DataIn = TFile::Open(("../data/mllTree_ZFullStat"+filename).c_str());
	DataIn->GetObject("tree",data);
	}else{
	if (full ){
		DataIn_prompt = TFile::Open(("../data/mllTree_"+filename).c_str());
		DataIn_Noprompt = TFile::Open(("../data/mllTree_JPsiNoPrompt"+filename).c_str());
		DataIn_prompt->GetObject("tree",data_prompt);
		DataIn_Noprompt->GetObject("tree",data_noprompt);
		}
	else {DataIn = TFile::Open(("../data/mllTree_"+lable+filename).c_str());
	DataIn->GetObject("tree",data);
	}
	}
	
	if(isZ) mcIn = TFile::Open(("../data/MCZ_fullStat/mllTreeMC_"+filename).c_str());
	else    mcIn = TFile::Open(("../data/MCJPsi_NoPrompt/mllTreeMC_"+filename).c_str());
	TTree* mc;
	mcIn->GetObject("treeMC",mc);

	if (full) lable = "JPsiTot";
	if (lable == "") lable ="JPsiPrompt";
	std::string PNGPATH = "/eos/home-r/ratramon/www/plots";
        gSystem->Exec(("mkdir "+PNGPATH+"/"+lable).c_str());		
        gSystem->Exec(("cp  "+PNGPATH+"/index.php " + PNGPATH+"/"+lable+"/.").c_str());		
	wspace.factory("nbkg[1000,100,100000000]");
	double lumi = 10.31;
	std::string x ;	
	if (isZ){
	 x = "mll_EE";
	wspace.factory((x+"[80,70,110]").c_str());
	}else{
	x = "mll_ET";
	wspace.factory((x+"[3,2.2,4]").c_str());
	}
	std::cout << "debug 1 " << std::endl;
//	std::cout << templ->GetEntries() << std::endl;
	wspace.factory("nsig[10000,0,10000000000]");
	wspace.factory("norm[100,20,1000000]");
	wspace.factory("NormC1[100,-0.000001,1000000]");
	wspace.factory("B_mass[5.3,5.2,5.4]");
	wspace.factory("NormC2[100,-0.000001,1000000]");
	wspace.factory("nmisreco[1000,0,1000000]");
	wspace.factory("nores[10,0,1000000]");
	wspace.factory("nVoigt[10,0.00,1000000]");
	wspace.factory("PUweight[0,-1.01,1000000]");
	wspace.var("mll_EE");
//	wspace.var("mll_ET1");
	wspace.var("mll_ET");
	wspace.var("B_mass");
	wspace.var("PUweight");
	wspace.var("mean");
	wspace.var("sigma");
	wspace.var("nsig");
	wspace.var("nbkg");
	wspace.var("Nnorm");
	wspace.var("nVoigt");
	wspace.var("NormC1");
	wspace.var("NOrmC2");
	wspace.var("nmisreco");
	wspace.var("alpha1");
	wspace.var("alpha2");
	wspace.var("n1");
	wspace.var("n2");
	wspace.var("nores");
	//RooDoubleCBShape DoubleCB("DoubleCBShape","DoubleCBShape",*wspace.var("mll_ET"),*wspace.var("mean"),*wspace.var("sigma"),*wspace.var("alpha1"),*wspace.var("n1"),*wspace.var("alpha2"),*wspace.var("n2"));
	//	wspace.import(DoubleCB);		
	//wspace.var("npred");
	//
//	RooRealVar mll_ET("mll_ET","mll_ET",2.2,4);
//	RooRealVar B_mass("B_mass","B_mass",5.2,5.4);
	//wspace.defineSet("TreeSet","mll_EE");
	std::string x_ext;
	if (!isZ)x_ext = x+",B_mass,PUweight";
	else x_ext = x+",PUweight";
	wspace.defineSet("TreeSet",(x_ext).c_str());
	std::string cut;
	if(isZ){
	cut = "";
	}else{
	
	cut = "(B_mass<5.4 && B_mass>5.2) ";
	}
	RooDataSet mcSet("mcSet","mcSet",*wspace.set("TreeSet"),RooFit::Import(*mc),RooFit::Cut((cut+" && PUweight!=-1 && PUweight<100 ").c_str()),RooFit::WeightVar("PUweight"));
	RooDataSet* dataSet;
	if(full && !isZ){
	dataSet = new RooDataSet("DataSet_prompt","DataSet_noprompt",*wspace.set("TreeSet"),RooFit::Import(*data_prompt),RooFit::Cut(cut.c_str()));
	RooDataSet dataSet_noprompt("DataSet","DataSet",*wspace.set("TreeSet"),RooFit::Import(*data_noprompt),RooFit::Cut(cut.c_str()));
	dataSet->append(dataSet_noprompt);
	}else
	dataSet = new RooDataSet("DataSet","DataSet",*wspace.set("TreeSet"),RooFit::Import(*data),RooFit::Cut(cut.c_str()));
	std::cout << "debug 2 " << std::endl;
//	wspace.var("mll_ET")->setRange("low",xmin,5.16);
//	wspace.var("mll_ET")->setRange("high",5.42,xmax);

	//std::cout << "histo stdev" << templ->GetStdDev() << std::endl; 
	if (sigPDF == 0){
		// Double CB
		wspace.factory("mean[3.0969+00, 3.08+00, 3.12+00]");

		wspace.factory("sigma[3e-2, 2.e-3, 8.0e-1]");
		wspace.factory("alpha1[1,0.001 ,10]");
		wspace.factory("alpha2[-7,-10, -0.001]");
		//	wspace.factory("alpha2[-0.5,-1 , 1.0e+2]");
		wspace.factory("n1[10, 0.001, 20]");
		wspace.factory("n2[3, 0.001, 40]");
		//	wspace.factory("n2[1, 0, 10]");
		wspace.factory(("CBShape::CB1("+x+",mean,sigma,alpha1,n1)").c_str());
		wspace.factory(("CBShape::CB2("+x+",mean,sigma,alpha2,n2)").c_str());
	//	wspace.factory("Gaussian::gaus("+x+",mean,sigma)");
	}
	if (sigPDF == 1){
		// Voigtian
		wspace.factory("mean[91.1e+00, 89.5e+00, 92e+00]");
		wspace.factory("width[2.4952, 2.4952, 2.4952]");
		//	wspace.factory("mean[3.0969+00, 3.06e+00, 3.10e+00]");
		wspace.factory("sigma[3.477e-02, 1.477e-5,2.5]");
		//wspace.factory("Gaussian::gaus("+x+",mean,sigma)");
		wspace.factory(("Voigtian::Voigt("+x+",mean,width,sigma)").c_str());
	}
	if (bkgPDF == 0){
		// Polynomial
		wspace.factory("c0[10, 0,100]");
		wspace.factory("c1[1, 0, 50]");
		wspace.factory("c2[0.03, 0, 50]");
		wspace.factory("c3[0, 0, 50]");
		wspace.factory(("Bernstein::bkg("+x+",{c0,c1,c2})").c_str());
	}
/*	if (bkgPDF == 1){
		wspace.factory("c1[0.0, -100.0, 100.0]");
		wspace.factory("Polynomial::bkg(mll_ET,{c1})");
	}*/
	if (bkgPDF == 2){
		// Emll_ETponential
		wspace.factory("exp_alpha[-10.0, -100.0, -0.001]");
		wspace.factory(("Exponential::bkg("+x+",exp_alpha)").c_str());
	}
/*	if (bkgPDF == 3){
		// Polynomial
		wspace.factory("c0[1.0, -1.0, 10.0]");
		wspace.factory("c1[0.1, -5.0, 5.0]");
		wspace.factory("c2[-0.1, -5.0, 5.0]");
		wspace.factory("c3[0.1, -5, 5]");
		wspace.factory("c4[-0.1, -5, 5]");
		wspace.factory("Bernstein::bkg(mll_ET,{c0,c1,c2,c3,c4})");
	}*/

	wspace.defineSet("obs",x.c_str());
	//wspace.defineSet("obs1","mll_ET1");
	RooArgList list = RooArgList(*wspace.set("obs"));


		if(isZ){
		wspace.factory("SUM::signal(nsig*Voigt)");
		wspace.factory("SUM::model(nsig*Voigt)");
		}else{
		wspace.factory("SUM::sig(NormC1*CB1,NormC2*CB2)");
		wspace.factory("SUM::signal(nsig*sig)");
		wspace.factory("SUM::model(nsig*sig,nbkg*bkg)");
		}
		wspace.factory("SUM::background(nbkg*bkg)");


	RooAbsPdf* mod = wspace.pdf("model");
	RooAbsPdf* bkg = wspace.pdf("bkg");
	RooAbsPdf* sig = wspace.pdf("sig");
	RooAbsPdf* signal = wspace.pdf("signal");
	
	if(isZ) RooFitResult* mc_template = signal->fitTo(mcSet,RooFit::Save(),RooFit::Range(88,95),RooFit::SumW2Error(kFALSE),RooFit::Extended(kTRUE));
	else{ 
//	wspace.var("mll_ET")->setRange("rFit",2.8,3.4);
	RooFitResult* mc_template = signal->fitTo(mcSet,RooFit::YVar(*wspace.var("PUweight")),RooFit::SumW2Error(kFALSE),RooFit::Save()/*,RooFit::Extended(kTRUE)*/);
	}
//	wspace.var("mll_ET")->setRange("r",2.2,4);
	resu[0]=wspace.var("mean")->getVal();
	resu[1]=wspace.var("mean")->getError();
	resu[4]=wspace.var("sigma")->getVal();
	resu[5]=wspace.var("sigma")->getError();
	if(isZ){
	wspace.var("width")->setConstant(kTRUE);
	}else{
	wspace.var("n1")->setConstant(kTRUE);
	wspace.var("n2")->setConstant(kTRUE);
	wspace.var("NormC1")->setConstant(kTRUE);
	wspace.var("NormC2")->setConstant(kTRUE);
//	wspace.var("width")->setConstant(kTRUE);
	wspace.var("alpha1")->setConstant(kTRUE);
	wspace.var("alpha2")->setConstant(kTRUE);
	}
	std::cout << " ________________________________________________________________________________________DEBUG " << dataSet->sumEntries() << std::endl;
	if(isZ)wspace.var("mean")->setVal(91);
	if(isZ)RooFitResult* result = mod->fitTo(*dataSet,RooFit::Save(),RooFit::Range(88,95),RooFit::Extended(kTRUE));
	else RooFitResult* result = mod->fitTo(*dataSet,RooFit::Save(),RooFit::Extended(kTRUE));
	resu[2] = wspace.var("mean")->getVal();
	resu[3] = wspace.var("mean")->getError();
	resu[6]=wspace.var("sigma")->getVal();
	resu[7]=wspace.var("sigma")->getError();
//	wspace.var("sigma")->setVal(wspace.var("sigma")->getVal());
//	wspace.var("sigma")->setConstant(kTRUE);
	
	TCanvas* c1 =new TCanvas("MC_Fit", "fit", 800, 1200);
	c1->Divide(1,2,0,0);
	c1->cd(1);
	c1->SetGrid();
	gPad->SetLeftMargin(0.10);
	gPad->SetRightMargin(0.05);
	gPad->SetBottomMargin(0.11);
	//gpad->setrightmargin(0.05);
	RooPlot* mll_ETmcframe = wspace.var(x.c_str())->frame();
	mcSet.plotOn(mll_ETmcframe, RooFit::Name("template data"),RooFit::DataError(RooAbsData::SumW2));
	wspace.var("mean")->setVal(resu[0]);
	wspace.var("sigma")->setVal(resu[4]);
	if (isZ)signal->plotOn(mll_ETmcframe,RooFit::Name("signal"),RooFit::Range(88,95),RooFit::LineColor(2),RooFit::MoveToBack()); // this will show fit overlay on canvas
	else signal->plotOn(mll_ETmcframe,RooFit::Name("signal"),RooFit::Range(2.2,4),RooFit::LineColor(2),RooFit::MoveToBack()); // this will show fit overlay on canvas
	mll_ETmcframe->GetYaxis()->SetTitleOffset(0.9);
	mll_ETmcframe->GetYaxis()->SetTitleFont(42);
	mll_ETmcframe->GetYaxis()->SetTitleSize(0.05);
	mll_ETmcframe->GetYaxis()->SetLabelSize(0.065);
	mll_ETmcframe->GetYaxis()->SetLabelSize(0.04);
	mll_ETmcframe->GetYaxis()->SetLabelFont(42);
	mll_ETmcframe->GetXaxis()->SetTitleOffset(0.9);
	mll_ETmcframe->GetXaxis()->SetTitleFont(42);
	mll_ETmcframe->GetXaxis()->SetTitleSize(0.05);
	mll_ETmcframe->GetXaxis()->SetLabelSize(0.065);
	mll_ETmcframe->GetXaxis()->SetLabelSize(0.04);
	mll_ETmcframe->GetXaxis()->SetLabelFont(42);

	mll_ETmcframe->GetYaxis()->SetTitle("Events");
	mll_ETmcframe->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV/c^{2}]");
	//mll_ETframe.GetXaxis().SetTitle("m(e^{+}e^{-}) [GeV/c^{2}]")
	mll_ETmcframe->SetStats(0);
	mll_ETmcframe->SetMinimum(0);
	mll_ETmcframe->Draw();
	c1->cd(2);
	c1->SetGrid();
	gPad->SetLeftMargin(0.10);
	gPad->SetRightMargin(0.05);
	RooPlot* mll_ETDataframe = wspace.var(x.c_str())->frame();
	wspace.var("mean")->setVal(resu[2]);
	wspace.var("sigma")->setVal(resu[6]);
	dataSet->plotOn(mll_ETDataframe, RooFit::Name("Dataset"));
	if(isZ)mod->plotOn(mll_ETDataframe,RooFit::Name("model"),RooFit::Range(88,95),RooFit::LineColor(2)); // this will show fit overlay on canvas
	else mod->plotOn(mll_ETDataframe,RooFit::Name("model"),RooFit::LineColor(2)); // this will show fit overlay on canvas
//	if (isZ) mod->plotOn(mll_ETDataframe,RooFit::Name("combinatorial"),RooFit::Components("bkg"),RooFit::DrawOption("FL"),RooFit::Range(69,111),RooFit::LineStyle(kDashed),RooFit::LineColor(kTeal),RooFit::FillColor(kTeal),RooFit::FillStyle(3001));
	if(!isZ)mod->plotOn(mll_ETDataframe,RooFit::Name("combinatorial"),RooFit::Components("bkg"),RooFit::DrawOption("FL"),RooFit::LineStyle(kDashed),RooFit::LineColor(kTeal),RooFit::FillColor(kTeal),RooFit::FillStyle(3001));
	mll_ETDataframe->GetYaxis()->SetTitleOffset(0.9);
	mll_ETDataframe->GetYaxis()->SetTitleFont(42);
	mll_ETDataframe->GetYaxis()->SetTitleSize(0.05);
	mll_ETDataframe->GetYaxis()->SetLabelSize(0.065);
	mll_ETDataframe->GetYaxis()->SetLabelSize(0.04);
	mll_ETDataframe->GetYaxis()->SetLabelFont(42);
	mll_ETDataframe->GetXaxis()->SetTitleOffset(0.9);
	mll_ETDataframe->GetXaxis()->SetTitleFont(42);
	mll_ETDataframe->GetXaxis()->SetTitleSize(0.05);
	mll_ETDataframe->GetXaxis()->SetLabelSize(0.065);
	mll_ETDataframe->GetXaxis()->SetLabelSize(0.04);
	mll_ETDataframe->GetXaxis()->SetLabelFont(42);

	mll_ETDataframe->GetYaxis()->SetTitle("Events");
	mll_ETDataframe->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV/c^{2}]");
	//mll_ETframe.GetXaxis().SetTitle("m(e^{+}e^{-}) [GeV/c^{2}]")
	mll_ETDataframe->SetStats(0);
	mll_ETDataframe->SetMinimum(0);
	mll_ETDataframe->Draw();
	


	CMS_lumi(c1,5,1.45,25);
//	if (isZ)c1->SaveAs((PNGPATH+"//fit_Z"+filename+".png").c_str());
	c1->SaveAs((PNGPATH+"/"+lable+"/"+"fit_"+filename+".png").c_str());
//	wspace.var("sigma1")->setConstant(kTRUE);
//	wspace.var("mean1")->setConstant(kTRUE);
//wspace.var("alpha2")->setConstant(kTRUE);
		//wspace.var("mean")->setConstant(kTRUE);
//		wspace.var("n1")->setConstant(kTRUE);
	
//		wspace.var("n2")->setConstant(kTRUE);
//	signal->plotOn(mll_ETmcframe,RooFit::Name("signal"),RooFit::Normalization(wspace.var("nsig")->getVal(),RooAbsReal::NumEvent),RooFit::LineColor(2),RooFit::MoveToBack()); // this will show fit overlay on canvas
//	mll_ETmcframe->Draw();
        //if(res)	c1->SaveAs((PNGPATH+"/plots/fit_MC"+std::to_string((int)set)+etaRegion+"_res.png").c_str());
       // else	c1->SaveAs((PNGPATH+"fit_MC"+std::to_string((int)set)+etaRegion+"_nores.png").c_str());
		
//	std::cout << "______________________________________________check" << wspace.var("nsig")->getVal() << std::endl;

//	std::cout << "sigma " << wspace.var("sigma")->getVal() << std::endl;
/*
	RooDataHist data = RooDataHist("data","data",list,hist);
	std::cout << hist->GetEntries() << std::endl;
	//RooAddPdf mod("model", "s+b",RooArgList(*wspace.obj("model")));
	
	//if (set <= 3)wspace.var("sigma")->setRange(templ->GetStdDev()*0.1,templ->GetStdDev()*0.75);
	wspace.var("sigma")->setRange(wspace.var("sigma")->getVal()*0.1,wspace.var("sigma")->getVal()*10);
//	wspace.var("sigma")->setVal(hist->GetStdDev()*0.5);
	wspace.var("mean")->setRange(3.06,3.12);
	wspace.var("mean")->setVal(3.09);
	wspace.var("mll_ET")->setRange("peak",2.6,3.5);
	RooFitResult* result;
	std::cout << " ____ "<< wspace.var("mean")->getVal()-3*wspace.var("sigma")->getVal()  << " " <<wspace.var("mean")->getVal()+3*wspace.var("sigma")->getVal() << std::endl;
	if (!res)result = mod->fitTo(data,RooFit::Save(),RooFit::Emll_ETtended(kTRUE) ,RooFit::Range("peak"));
	else result = mod->fitTo(data,RooFit::Save(),RooFit::Emll_ETtended(kTRUE));
	std::cout << " _____________________________________________________________________________________________________DEBUG "  << std::endl;
	result->Print(); 
	if(wspace.var("nsig")->getError()>0.8*wspace.var("nsig")->getVal()){

	std::cout << "_______________________________________________________Signal yield error too large, repeating fit with new paramaters" << std::endl; 
	result = mod->fitTo(data,RooFit::Save(),RooFit::Emll_ETtended(kTRUE)/* ,RooFit::Range("low","high"));
	}
	std::cout << "______________________________________________check" << wspace.var("nsig")->getVal() << std::endl;
	result->Print(); 
	std::cout << "sigma " << wspace.var("sigma")->getVal() << std::endl;
	wspace.var("mll_ET")->setRange("window",wspace.var("mean")->getVal()-3*wspace.var("sigma")->getVal(),wspace.var("mean")->getVal()+3*wspace.var("sigma")->getVal());
	RooAbsReal* fracBkgRange = bkg->createIntegral(*wspace.set("obs"),*wspace.set("obs"),"window") ;
	std::cout << "____________________________________________________________________________________________________debug"  << fracBkgRange<< std::endl;
	RooAbsReal* fracSigRange;
	fracSigRange = sig->createIntegral(*wspace.set("obs"),*wspace.set("obs"),"window") ;
//	else  fracSigRange = sig->createIntegral(*wspace.set("obs"),*wspace.set("obs"),"window") ;
	double  nbkgWindow = wspace.var("nbkg")->getVal() * fracBkgRange->getVal();
	double  nSigWindow = wspace.var("nsig")->getVal() * fracSigRange->getVal();
	// 	print(nbkg.getVal(), fracBkgRange.getVal());
	std::cout << "Number of signals:" << nSigWindow << ", Number of background: " << nbkgWindow << ", S/sqrt(S+B): " << nSigWindow/sqrt(nSigWindow + nbkgWindow)<< std::endl;
	TCanvas* c2 =new TCanvas("fig_binnedFit", "fit", 800, 600);
	c2->SetGrid();
	c2->cd();
	gPad->SetLeftMargin(0.10);
	gPad->SetRightMargin(0.05);
	RooPlot* mll_ETframe = wspace.var("x")->frame();
	gStyle->SetPalette(kPastel);


	if (!res){
		data.plotOn(mll_ETframe, RooFit::Name("data"));
		mod->plotOn(mll_ETframe,RooFit::Name("combinatorial"),RooFit::Components("bkg")/*,RooFit::Normalization(wspace.var("nbkg")->getVal(),RooAbsReal::NumEvent),RooFit::DrawOption("L"),RooFit::Range(xmin,xmax),RooFit::LineStyle(kDashed),RooFit::LineColor(kAzure+7));
		mod->plotOn(mll_ETframe,RooFit::Name("comb"),RooFit::Components("bkg")/*,RooFit::Normalization(wspace.var("nbkg")->getVal(),RooAbsReal::NumEvent),RooFit::DrawOption("FL"),RooFit::Range(xmin,xmax),RooFit::LineStyle(kDashed),RooFit::LineColor(kAzure+7),RooFit::FillColor(kAzure-3),RooFit::FillStyle(3001));
//		mod->plotOn(mll_ETframe,RooFit::Name("comb"),RooFit::Components("bkg")/*,RooFit::Normalization(wspace.var("nbkg")->getVal(),RooAbsReal::NumEvent),RooFit::Range(5.40,xmax),RooFit::DrawOption("FL"),RooFit::LineStyle(kDashed),RooFit::LineColor(kAzure+7),RooFit::FillColor(kAzure-3),RooFit::FillStyle(3001));
		mod->plotOn(mll_ETframe,RooFit::Name("model"),RooFit::LineColor(2)/*,RooFit::Normalization(wspace.var("nbkg")->getVal()+wspace.var("nmisreco")->getVal(),RooAbsReal::NumEvent)); // this will show fit overlay on canvas
	}else{

		data.plotOn(mll_ETframe, RooFit::Name("data"));
		mod->plotOn(mll_ETframe,RooFit::Name("model"),RooFit::LineColor(2)); // this will show fit overlay on canvas
		//mod ->plotOn(mll_ETframe,RooFit::Name("model"),RooFit::Components("sig"),RooFit::DrawOption("FL"),RooFit::FillColor(9),RooFit::FillStyle(3004),RooFit::LineStyle(6),RooFit::LineColor(9)) ;
		mod->plotOn(mll_ETframe,RooFit::Name("combinatorial"),RooFit::Components("bkg"),RooFit::DrawOption("L"),RooFit::FillColor(kMagenta),RooFit::FillStyle(3001),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
//		mod->plotOn(mll_ETframe,RooFit::Name("misreco"),RooFit::Components("misr"),RooFit::DrawOption("FL"),RooFit::AddTo("combinatorial",1,1),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3),RooFit::FillColor(kOrange-3),RooFit::FillStyle(3001));
		mod->plotOn(mll_ETframe,RooFit::Name("comb"),RooFit::Components("bkg"),RooFit::DrawOption("FL"),RooFit::FillColor(kAzure+7),RooFit::FillStyle(3001),RooFit::LineStyle(kDashed),RooFit::LineColor(kAzure+7));
		mod ->plotOn(mll_ETframe,RooFit::Name("Sig"),RooFit::Components("sig"),RooFit::DrawOption("FL"),RooFit::FillColor(kRed-3),RooFit::FillStyle(3004),RooFit::LineColor(kRed-3)) ;
		std::cout << "________________________________________________________" << wspace.var("nsig")->getVal() << std::endl;
	}
	//	mod->plotOn(mll_ETframe,RooFit::VisualizeError(result), RooFit::FillColor(kOrange), RooFit::MoveToBack())// # this will show fit overlay on canvas
	//pdf->plotOn(frame,Normalization(0.5,RooAbsReal::Relative)) ;
	
	mll_ETframe->GetYaxis()->SetTitleOffset(0.9);
	mll_ETframe->GetYaxis()->SetTitleFont(42);
	mll_ETframe->GetYaxis()->SetTitleSize(0.05);
	mll_ETframe->GetYaxis()->SetLabelSize(0.065);
	mll_ETframe->GetYaxis()->SetLabelSize(0.04);
	mll_ETframe->GetYaxis()->SetLabelFont(42);
	mll_ETframe->GetXaxis()->SetTitleOffset(0.9);
	mll_ETframe->GetXaxis()->SetTitleFont(42);
	mll_ETframe->GetXaxis()->SetTitleSize(0.05);
	mll_ETframe->GetXaxis()->SetLabelSize(0.065);
	mll_ETframe->GetXaxis()->SetLabelSize(0.04);
	mll_ETframe->GetXaxis()->SetLabelFont(42);

	mll_ETframe->GetYaxis()->SetTitle("Events");
	mll_ETframe->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV/c^{2}]");
	//mll_ETframe.GetXaxis().SetTitle("m(e^{+}e^{-}) [GeV/c^{2}]")
	mll_ETframe->SetStats(0);
	mll_ETframe->SetMinimum(0);
	mll_ETframe->Draw();
	mB_res->SetLineColor(kOrange+8);
	mB_res->SetLineWidth(3);
	mB_res->Draw("samehist");

	TLegend* l;
	l = new TLegend(0.67,0.75,0.95,0.90);
	//else  l = new TLegend(0.67,0.15,0.95,0.50);
	//	l->SetTemll_ETtFont(72);
	l->SetTemll_ETtSize(0.03);
	l->AddEntry(mll_ETframe->findObject("data"),"Data","lpe");
	l->AddEntry(mll_ETframe->findObject("comb"),"Combinatorial bkg","l");
	l->AddEntry(mll_ETframe->findObject("Sig"),"Signal","l");

	char v1[10],v2[10],v3[10];
	sprintf(v1,"%.2f",wspace.var("nsig")->getVal());
	sprintf(v2,"%.2f", nbkgWindow);
	sprintf(v3,"%.2f", nores);
	if (res) l->AddEntry("" ,("S = "+std::string(v1)).c_str(),"");
	else {
		l->AddEntry("" ,("S_pred = "+std::string(v3)).c_str(),"");
		//	l->AddEntry("" ,("S_fit = "+std::string(v1)).c_str(),"");


	}
	l->AddEntry("",("B = "+std::string(v2)).c_str(),"");
	char v[10];
	if(!res)sprintf(v,"%.2f",nores/sqrt(nores + nbkgWindow));
	if(res)sprintf(v,"%.2f",nSigWindow/sqrt(nSigWindow + nbkgWindow));
	l->AddEntry("" ,("S/#sqrt{S+B} = "+std::string(v)).c_str(),"");
	l->Draw();
	CMS_lumi(c2,5,0,25);
	if (res)c2->SaveAs((PNGPATH+"/plots/fit_peak"+std::to_string((int)set)+etaRegion+".png").c_str());
	else c2->SaveAs((PNGPATH+"fit_nores"+std::to_string((int)set)+".png").c_str());
	resu[2] = wspace.var("mean")->getVal();
	resu[3] = wspace.var("mean")->getError();
	resu[6]=wspace.var("sigma")->getVal();
	resu[7]=wspace.var("sigma")->getError();
	std::cout << " S/sqrt(S+B): " << nSigWindow/sqrt(nSigWindow + nbkgWindow)<< std::endl;
	std::pair<double,double> pair = std::make_pair(nSigWindow,nbkgWindow);
	std::cout << " S " << pair.first<< std::endl;
	std::cout << " S/sqrt(S+B): " << pair.second<< std::endl;
//	if(!res){double err = PoissonError(wspace);
	//std::cout <<  "_______________________________________eval" << err << std::endl;}
	delete mB_res;
	delete mB_nores;
	delete templ;*/

}


void setStyle() {


	// set the TStyle
	TStyle* style = new TStyle("DrawBaseStyle", "");
	style->SetCanvasColor(0);
	style->SetPadColor(0);
	style->SetFrameFillColor(0);
	style->SetStatColor(0);
	style->SetOptStat(0);
	style->SetOptFit(0000);
	style->SetTitleFillColor(0);
	style->SetCanvasBorderMode(0);
	style->SetPadBorderMode(0);
	style->SetFrameBorderMode(0);
	style->SetPadBottomMargin(0.12);
	style->SetPadLeftMargin(0.12);
	style->cd();
	// For the canvas:
	style->SetCanvasBorderMode(0);
	style->SetCanvasColor(kWhite);
	style->SetCanvasDefH(600); //Height of canvas
	style->SetCanvasDefW(600); //Width of canvas
	style->SetCanvasDefX(0); //POsition on screen
	style->SetCanvasDefY(0);
	// For the Pad:
	style->SetPadBorderMode(0);
	style->SetPadColor(kWhite);
	style->SetPadGridX(false);
	style->SetPadGridY(false);
	style->SetGridColor(0);
	style->SetGridStyle(3);
	style->SetGridWidth(1);
	// For the frame:
	style->SetFrameBorderMode(0);
	style->SetFrameBorderSize(1);
	style->SetFrameFillColor(0);
	style->SetFrameFillStyle(0);
	style->SetFrameLineColor(1);
	style->SetFrameLineStyle(1);
	style->SetFrameLineWidth(1);
	// Margins:
	style->SetPadTopMargin(0.10);
	style->SetPadBottomMargin(0.14);//0.13);
	style->SetPadLeftMargin(0.16);//0.16);
	style->SetPadRightMargin(0.1);//0.02);
	// For the Global title:
	style->SetOptTitle(0);
	style->SetTitleFont(42);
	style->SetTitleColor(1);
	style->SetTitleTextColor(1);
	style->SetTitleFillColor(10);
	style->SetTitleFontSize(0.05);
	// For the axis titles:
	style->SetTitleColor(1, "XYZ");
	style->SetTitleFont(42, "XYZ");
	style->SetTitleSize(0.05, "XYZ");
	style->SetTitleXOffset(1.15);//0.9);
	style->SetTitleYOffset(1.5); // => 1.15 if exponents
	// For the axis labels:
	style->SetLabelColor(1, "XYZ");
	style->SetLabelFont(42, "XYZ");
	style->SetLabelOffset(0.007, "XYZ");
	style->SetLabelSize(0.045, "XYZ");
	// For the axis:
	style->SetAxisColor(1, "XYZ");
	style->SetStripDecimals(kTRUE);
	style->SetTickLength(0.03, "XYZ");
	style->SetNdivisions(510, "XYZ");
	style->SetPadTickX(1); // To get tick marks on the opposite side of the frame
	style->SetPadTickY(1);
	// for histograms:
	style->SetHistLineColor(1);
	// for the pallete
	Double_t stops[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red  [5] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[5] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue [5] = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(5, stops, red, green, blue, 100);
	style->SetNumberContours(100);

	style->cd();

}


void DoubleSlicer(int bin,float min, float max, bool isZ, bool full){

	int i,k;
	float x[2][bin],mean[2][bin],RMS[2][bin];
	float peak_diff[2][bin],peak_diff_sigma[2][bin],pt[2][bin],pt_err[2][bin],sigma_ratio[2][bin],sigma_ratio_unc[2][bin];

	std::string barrel[2];
	std::string lable[2];
	std::string GT[2];
	std::ofstream outfile[2][2];
	barrel[0] = "BarrelCentral";
	barrel[1] = "BarrelEdge";
	lable[0] = "";
	lable[1] = "JPsiNoPrompt";
	GT[0] = " GT: 102X_dataRun2_Prompt_v14";
	GT[1] = " GT: 102X_dataRun2_v11";
	double resu[8];
	
	int nPtBin = 4;
	double min_edge[4]= {10,40,50,75};
	double max_edge[4]= {40,50,75,150};
	int EtaBin;
	if (isZ) EtaBin = 2;
	else EtaBin = 1;
	int color[2]={kRed-7, kCyan-3};

	for(k=0;k<EtaBin;k++){

	if(isZ){
	outfile[k][0].open(("linearity_Z_"+barrel[k]+".txt").c_str());
	outfile[k][1].open(("NormSigma_Z_"+barrel[k]+".txt").c_str());
	}else{
	outfile[k][0].open(("linearity_JPsi_"+barrel[k]+".txt").c_str());
	outfile[k][1].open(("NormSigma_JPsi_"+barrel[k]+".txt").c_str());
	}
	for(i=0; i<bin;i++){
	//	if(isZ && k==0)x[k][i] = min_edge[i];	
		x[k][i]=min +(max-min)/bin*i;
		char range[15]="";
		int n;
		n=sprintf(range,"%.3f;%.3f",x[k][i],x[k][i]+(max-min)/bin);
		std::string name;
//		if (isZ) name= "PtMin"+std::to_string((int)x[k][i])+"PtMax"+std::to_string((int)(max_edge[i]))+".root";
		name= "PtMin"+std::to_string((int)x[k][i])+"PtMax"+std::to_string((int)(x[k][i]+(max-min)/bin))+barrel[k]+".root";
		if (isZ)fit(name,lable[k],1,2,1,full,resu);
		else fit(name,lable[k],0,0,0,full,resu);

		double m_pdg;
		if (isZ) m_pdg = 91.1876;
		else m_pdg = MASS_JPSI;
		peak_diff[k][i]=(resu[2]-resu[0])/m_pdg;
		
		std::cout << "check: " << peak_diff[i] << std::endl;
		peak_diff_sigma[k][i]= sqrt(pow(resu[1],2)+pow(resu[3],2))/m_pdg;
		sigma_ratio[k][i] = resu[6]/resu[4];
		sigma_ratio_unc[k][i] = sqrt(pow(resu[7]/resu[4],2)+pow(resu[6]*resu[5]/(resu[4]*resu[4]),2));
	//	if (isZ)pt[k][i]= min_edge[i]+(max_edge[i]-min_edge[i])/2;
	//	else
		  pt[k][i]=x[k][i]+(max-min)/(2*bin);
	//	if (isZ)
	//	pt_err[k][i] = (max_edge[i]-min_edge[i])/2;
	//	else
		 pt_err[k][i] = (max-min)/(2*bin);
		//SavePlot (xaxis, proj[k][i],(PLOTPATH+filename+std::string(range)).c_str() ,false);
		outfile[k][0] << pt[k][i] << " "  << peak_diff[k][i] << " " << pt_err[k][i] <<  " " << peak_diff_sigma[k][i] << std::endl;
		outfile[k][1] << pt[k][i] << " "  << sigma_ratio[k][i] << " " << pt_err[k][i] <<  " " << sigma_ratio_unc[k][i] << std::endl;
	}

	outfile[k][0].close();
	outfile[k][1].close();
	}

	
	setStyle();
	TGraphErrors* lin[2];
	TGraphErrors* sigma[2];
	//std::cout << "check: " << peak_diff[0] << " " << peak_diff[1] << " " <<peak_diff[2] << " " <<peak_diff[3] << " " <<peak_diff[4]<< " " <<peak_diff[5]<< " " << MASS_JPSI << std::endl;
	lin[0] = new TGraphErrors(bin,pt[0],peak_diff[0],pt_err[0],peak_diff_sigma[0]);
	lin[1] = new TGraphErrors(bin,pt[1],peak_diff[1],pt_err[1],peak_diff_sigma[1]);
	TH2D* plotter;
	if(isZ)plotter = new TH2D("plotter","plotter",10,0,170,10,-0.01,0.01);
	else plotter = new TH2D("plotter","plotter",10,0,20,10,-0.01,0.01);
	TCanvas* canvas = new TCanvas("","",800,600);
        TLegend* l;  
   	if(!full)  l = new TLegend(0.3,0.9,0.80,0.75);
         else  l = new TLegend(0.3,0.88,0.60,0.8);
    	  l->SetFillStyle(0);
    	  l->SetLineStyle(0);
    	  l->SetBorderSize(0);
	  plotter->GetXaxis()->SetTitle("pt(Gev/c)");
	  plotter->GetYaxis()->SetTitle("(m^{peak}_{data}-m^{peak}_{MC})/m_{PDG}");
	  plotter->Draw();
	  for(k=0;k<2;k++){
	 	 lin[k]->SetMarkerStyle(8);
	 	 lin[k]->SetMarkerSize(1.5);
		 lin[k]->SetLineWidth(2);
		 lin[k]->SetMarkerColor(color[k]);
		 lin[k]->SetLineColor(color[k]);
	 	 lin[k]->Draw("sameP");
		 if (!full)l->AddEntry(lin[k],"Z ","lp");

	}
	if (full)l->AddEntry(lin[1],"Bparking 40 fb-1","lp");
	l->Draw(); 
	if (isZ) canvas->SaveAs("linearity_Z.pdf");
	else if (full) canvas->SaveAs("linearity_JPsiTot.pdf");
	else  canvas->SaveAs("linearity_JPsi.pdf");


	sigma[0] = new TGraphErrors(bin,pt[0],sigma_ratio[0],pt_err[0],sigma_ratio_unc[0]);
	sigma[1] = new TGraphErrors(bin,pt[1],sigma_ratio[1],pt_err[1],sigma_ratio_unc[1]);
	TH2D* plotter1;
	if(isZ)plotter1 = new TH2D("plotter1","plotter1",10,0,170,10,0,2);
	else plotter1 = new TH2D("plotter1","plotter1",10,0,20,10,0,2);
	  TCanvas* canva1 = new TCanvas("1","1",800,600);
          TLegend* l1 = new TLegend();
    	  l1->SetFillStyle(0);
    	  l1->SetLineStyle(0);
    	  l1->SetBorderSize(0);
	  plotter1->GetXaxis()->SetTitle("pt(Gev/c)");
	  plotter1->GetYaxis()->SetTitle("(#sigma_{data}/#sigma_{MC}");
	  plotter1->Draw();
	  for(k=0;k<2;k++){
	 	sigma[k]->SetMarkerStyle(8);
	 	sigma[k]->SetMarkerSize(1.5);
		sigma[k]->SetLineWidth(1);
		sigma[k]->SetMarkerColor(color[k]);
		sigma[k]->SetLineColor(color[k]);
	    	sigma[k]->Draw("sameP");
		l1->AddEntry(sigma[k],"Z","lp");

	}
	l1->Draw();
	if (isZ) canva1->SaveAs("sigma_Z.pdf");
	else if (full) canva1->SaveAs("sigma_JPsiTot.pdf");
	else  canva1->SaveAs("sigma_JPsi.pdf");
	  

}
