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
#include "TString.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TTree.h"


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



void fit(std::string filename, int  sigPDF,int bkgPDF, bool res, int set,std::string etaRegion){

	std::cout << "in fit________________________ " << std::endl;
	RooWorkspace wspace("w");
	gStyle->SetOptFit(0000);
	gROOT->SetBatch(true);
	gROOT->SetStyle("Plain");
	gStyle->SetGridStyle(3);
	gStyle->SetOptStat(000000);
	gStyle->SetOptTitle(0);

	std::string PNGPATH = "/eos/home-r/ratramon/www/";
	bool blind=true ;
	double eff_res, eff_nores,eff_misreco;
	double xmin, xmax;


	double resu[8];
	
	//mll_ETmin = hist->GetXaxis()->GetXmin()/*+2*hist->GetXaxis()->GetBinWidth(1)*/;
//	mll_ETmax = hist->GetXaxis()->GetXmax();//-2*hist->GetXaxis()->GetBinWidth(1);
//	std::cout << "bin width " << hist->GetXamll_ETis()->GetBinWidth(1) << std::endl;
//	TH1D* mB_res = new TH1D("sig_norm","sig_norm",(int)((mll_ETmax-xmin)/hist->GetXaxis()->GetBinWidth(1)), xmin, xmax); 
//	TH1D* mB_nores = new TH1D("sig_norm","sig_norm",(int)((mll_ETmax-xmin)/hist->GetXaxis()->GetBinWidth(1)), xmin, xmax); 

//	std::cout << "limits" << mll_ETmin << "   " << xmax << std::endl;

	//wspace.factory(("mll_ET[5.0,"+std::to_string(xmin)+","+std::to_string(xmax)+"]").c_str());

	TFile* DataIn = TFile::Open(("../data/mllTree_"+filename).c_str());
	TTree* data;
	DataIn->GetObject("tree",data);
	TFile* mcIn = TFile::Open(("../data/mllTreeMC_"+filename).c_str());
	TTree* mc;
	mcIn->GetObject("treeMC",mc);
	wspace.factory("nbkg[100000,0,1000000]");
	double lumi = 10.31;
	int c =set;

//	std::cout << templ->GetEntries() << std::endl;
	wspace.factory("nsig[10000,0,1000000]");
	wspace.factory("norm[100,20,1000000]");
	wspace.factory("NormC1[100,0,1000000]");
	wspace.factory("mll_ET[3,2.2,4]");
	wspace.factory("B_mass[5.3,5.2,5.4]");
	wspace.factory("NormC2[100,-0.000001,1000000]");
	wspace.factory("nmisreco[1000,0,1000000]");
	wspace.factory("nores[10,0,1000000]");
	wspace.factory("nVoigt[10,0.00,1000000]");
//	wspace.var("mll_ET");
//	wspace.var("mll_ET1");
	wspace.var("mll_ET");
	wspace.var("B_mass");
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
	wspace.defineSet("TreeSet","mll_ET,B_mass");
	RooDataSet mcSet("mcSet","mcSet",*wspace.set("TreeSet"),RooFit::Import(*mc),RooFit::Cut("B_mass<5.4 && B_mass>5.2"));
	RooDataSet dataSet("DataSet","DataSet",*wspace.set("TreeSet"),RooFit::Import(*data),RooFit::Cut("B_mass<5.4 && B_mass>5.2"));

//	wspace.var("mll_ET")->setRange("low",xmin,5.16);
//	wspace.var("mll_ET")->setRange("high",5.42,xmax);

	//std::cout << "histo stdev" << templ->GetStdDev() << std::endl; 
	if (sigPDF == 0){
		// Voigtian
		wspace.factory("width[1.000e-02, 9.000e-03, 1.000e-01]");
		wspace.factory("mean[3.02+00, 3.01+00, 3.12e+00]");
		wspace.factory("mean1[2.8+00, 2.55e+00, 3.1200e+00]");

		wspace.factory("sigma1[2.6e-1, 1.e-2, 3.01]");
	//	wspace.factory("mean[5.279e+00, 5.27e+00, 5.30e+00]");
	//	wspace.factory(("mean["+std::to_string(templ->GetXamll_ETis()->GetBinCenter(templ->GetMaximumBin()))+","+std::to_string(templ->GetXaxis()->GetBinCenter(templ->GetMaximumBin())-.5*templ->GetStdDev())+", "+std::to_string(templ->GetXaxis()->GetBinCenter(templ->GetMaximumBin())+-5*templ->GetStdDev())+"]").c_str());
	//	wspace.factory(("sigma["+std::to_string(templ->GetStdDev())+","+std::to_string(0.1*templ->GetStdDev())+","+std::to_string(5*templ->GetStdDev())+"]").c_str());
		wspace.factory("sigma[2e-2, 1.e-5, 5e-1]");
	wspace.factory("Voigtian::Voigt(mll_ET,mean,width,sigma)");
		wspace.factory("Gaussian::gaus(mll_ET,mean,sigma)");
		wspace.factory("Gaussian::gaus1(mll_ET,mean1,sigma1)");
	}
	if (sigPDF == 1){
		// Gaussian
		wspace.factory("mean[5.2418e+00, 5.20e+00, 5.35e+00]");
		//	wspace.factory("mean[3.0969+00, 3.06e+00, 3.10e+00]");
		wspace.factory("sigma[3.477e-02, 1.477e-02, 7.477e-02]");
			wspace.factory("Gaussian::sig(mll_ET,mean,sigma)");
	}
	if (sigPDF == 0){
		// Crystal-ball
		//	wspace.factory("mean[5.279e+00, 5.26e+00, 5.29e+00]");
		//	wspace.factory("mean[3.0969+00, 3.06e+00, 3.10e+00]");
		//	wspace.factory("sigma[3.1858e-02, 1.e-3, 7.e-2]");
		//		wspace.factory("sA[0, -1.e-3, 5.e-2]");
		wspace.factory("alpha1[4,0.001 ,100]");
		wspace.factory("alpha2[-7,-100, -0.001]");
		//	wspace.factory("alpha2[-0.5,-1 , 1.0e+2]");
		wspace.factory("n1[10, 0.001, 20]");
		wspace.factory("n2[3, 0.001, 40]");
		//	wspace.factory("n2[1, 0, 10]");
		wspace.factory("CBShape::CB1(mll_ET,mean,sigma,alpha1,n1)");
		wspace.factory("CBShape::CB2(mll_ET,mean,sigma,alpha2,n2)");
	}
	if (bkgPDF == 0){
		// Polynomial
		wspace.factory("c0[1, 0,10]");
		wspace.factory("c1[0, 0, 5]");
		wspace.factory("c2[0.03, 0, 5]");
		wspace.factory("c3[0, 0, 5]");
		wspace.factory("Bernstein::bkg(mll_ET,{c0,c1,c2})");
	}
	if (bkgPDF == 1){
		wspace.factory("c1[0.0, -100.0, 100.0]");
		wspace.factory("Polynomial::bkg(mll_ET,{c1})");
	}
	if (bkgPDF == 2){
		// Emll_ETponential
		wspace.factory("emll_ETp_alpha[-10.0, -100.0, -1]");
		wspace.factory("Emll_ETponential::bkg(x,exp_alpha)");
	}
	if (bkgPDF == 3){
		// Polynomial
		wspace.factory("c0[1.0, -1.0, 10.0]");
		wspace.factory("c1[0.1, -5.0, 5.0]");
		wspace.factory("c2[-0.1, -5.0, 5.0]");
		wspace.factory("c3[0.1, -5, 5]");
		wspace.factory("c4[-0.1, -5, 5]");
		wspace.factory("Bernstein::bkg(mll_ET,{c0,c1,c2,c3,c4})");
	}

	wspace.defineSet("obs","mll_ET");
	//wspace.defineSet("obs1","mll_ET1");
	RooArgList list = RooArgList(*wspace.set("obs"));


/*	if (res)
	{*/
		wspace.factory("SUM::sig(NormC1*CB1,NormC2*CB2)");
		wspace.factory("SUM::signal(nsig*sig)");
	//	wspace.factory("SUM::sig(nVoigt*Voigt)");
		wspace.factory("SUM::model(nsig*sig,nbkg*bkg)");
	//	wspace.factory("SUM::signal()");
		wspace.factory("SUM::background(nbkg*bkg)");
/*	}
	else{ 
		//wspace.factory("SUM::sig(pippo*CB1,mario*CB2)");
//		wspace.factory("mean[5.279, 5.20e+00, 5.35e+00]");
//		wspace.factory("sigma[7.477e-02, 7.477e-02, 7.477e-02]");
	//	wspace.factory("Gaussian::signal(mll_ET,mean,sigma)");
		//wspace.factory("SUM::model(nbkg*bkg, nmisreco*misr)");
		wspace.factory("SUM::sig(pippo*CB1,mario*CB2)");
		wspace.factory("SUM::model(nsig*sig,nbkg*bkg");
		wspace.factory("SUM::signal(nsig*sig)");
//wspace.factory("SUM::signal(nsig*sig)");
	}*/
	RooAbsPdf* signal= wspace.pdf("signal");
	RooAbsPdf* mod = wspace.pdf("model");
	RooAbsPdf* bkg = wspace.pdf("bkg");
	RooAbsPdf* background = wspace.pdf("background");
	RooAbsPdf* sig = wspace.pdf("sig");
	//bkg = wspace.pdf("bkg");
	std::cout << "sigma " << wspace.var("sigma")->getVal() << std::endl;

	RooFitResult* template_par = signal->fitTo(mcSet,RooFit::Save()/*,RooFit::Extended(kTRUE)*/);
	wspace.var("n1")->setConstant(kTRUE);
	wspace.var("n2")->setConstant(kTRUE);
	wspace.var("NormC1")->setConstant(kTRUE);
	wspace.var("NormC2")->setConstant(kTRUE);
//	wspace.var("width")->setConstant(kTRUE);
	wspace.var("alpha1")->setConstant(kTRUE);
	wspace.var("alpha2")->setConstant(kTRUE);
	std::cout << " _____________________________________________________________________________________________________DEBUG "  << std::endl;
	resu[0]=wspace.var("mean")->getVal();
	resu[1]=wspace.var("mean")->getError();
	resu[4]=wspace.var("sigma")->getVal();
	resu[5]=wspace.var("sigma")->getError();
	RooFitResult* result = mod->fitTo(dataSet,RooFit::Save(),RooFit::Extended(kTRUE));
	resu[2] = wspace.var("mean")->getVal();
	resu[3] = wspace.var("mean")->getError();
	resu[6]=wspace.var("sigma")->getVal();
	resu[7]=wspace.var("sigma")->getError();
	double data_mean = wspace.var("mean")->getVal();
	double data_sigma = wspace.var("sigma")->getVal();
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
	RooPlot* mll_ETmcframe = wspace.var("mll_ET")->frame();
	mcSet.plotOn(mll_ETmcframe, RooFit::Name("template data"));
	wspace.var("mean")->setVal(resu[0]);
	wspace.var("sigma")->setVal(resu[4]);
	signal->plotOn(mll_ETmcframe,RooFit::Name("signal"),RooFit::LineColor(2),RooFit::MoveToBack()); // this will show fit overlay on canvas
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
	RooPlot* mll_ETDataframe = wspace.var("mll_ET")->frame();
	wspace.var("mean")->setVal(resu[2]);
	wspace.var("sigma")->setVal(resu[6]);
	dataSet.plotOn(mll_ETDataframe, RooFit::Name("Dataset"));
	mod->plotOn(mll_ETDataframe,RooFit::Name("model"),RooFit::LineColor(2)); // this will show fit overlay on canvas
	mod->plotOn(mll_ETDataframe,RooFit::Name("combinatorial"),RooFit::Components("bkg"),RooFit::DrawOption("FL"),RooFit::LineStyle(kDashed),RooFit::LineColor(kTeal),RooFit::FillColor(kTeal),RooFit::FillStyle(3001));
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
	c1->SaveAs((PNGPATH+"/plots/fit_MCData"+filename+".png").c_str());
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
