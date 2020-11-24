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


#define MASS_JPSI 3.0969

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
			latex.DrawLatex(posX_*(5), posY_, cmsText);
//			latex.DrawLatex(posX_, posY_-posY_*0.49, cmsText);
			if( writeExtraText ) 
			{
				latex.SetTextFont(extraTextFont);
				latex.SetTextAlign(align_);
				latex.SetTextSize(extraTextSize*t);
			//	latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, cmsText);
		//		latex.DrawLatex(posX_, posY_-0.03*posY_, MCText);
				latex.DrawLatex(posX_*5, posY_*0.97, extraText);
				latex.SetTextFont(lumiTextFont);
				latex.DrawLatex(posX_*5, posY_*0.94,lumiText+" Bparking");
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

void fit(std::string filename, std::string lable, std::string cut_block, int  sigPDF,int bkgPDF, bool isZ, int set,double * resu){

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
	TFile* DataIn;
	TFile* mcIn;
	if(isZ) DataIn = TFile::Open(("../data/"+filename).c_str());
	else    DataIn = TFile::Open(("../data/"+filename).c_str());
	TTree* data;
	DataIn->GetObject("tree",data);
//	if(isZ) mcIn = TFile::Open(("../data/MCZ_fullStat/mllTreeMC_"+filename).c_str());
//	else    mcIn = TFile::Open(("../data/MCJPsi/mllTreeMC_"+filename).c_str());
//	TTree* mc;
//	mcIn->GetObject("treeMC",mc);
	if (lable == "") lable ="JPsiPrompt";
	std::string PNGPATH = "/eos/home-r/ratramon/www/plots";
        gSystem->Exec(("mkdir "+PNGPATH+"/"+lable).c_str());		
        gSystem->Exec(("cp  "+PNGPATH+"/index.php " + PNGPATH+"/"+lable+"/.").c_str());		
	wspace.factory("nbkg[1000,100,100000000]");
	double lumi = 10.31;
	int c =set;
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
	wspace.factory("nsig[10000,0,1000000]");
	wspace.factory("luminosityBlock[100,0,2500]");
	wspace.factory("run[100,0,325500]");
	wspace.factory("norm[100,20,1000000]");
	wspace.factory("NormC1[1000,0,1000000]");
	wspace.factory("B_mass[5.3,5.2,5.4]");
	wspace.factory("NormC2[1000,0.00,1000000]");
	wspace.factory("nmisreco[1000,0,1000000]");
	wspace.factory("nores[10,0,1000000]");
	wspace.factory("nVoigt[10,0.00,1000000]");
	wspace.factory("PUweight[0,-1.01,1000000]");
	wspace.var("mll_EE");
	wspace.var("luminosityBlock");
	wspace.var("run");
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
	if (!isZ)x_ext = x+",B_mass,PUweight,luminosityBlock,run";
	else x_ext = x+",PUweight,luminosityBlock,run";
	wspace.defineSet("TreeSet",(x_ext).c_str());
//	wspace.defineSet("TreeSet","mll_ET,B_mass,luminosityBlock");
	std::string cut;
	if(isZ){
	cut = "("+cut_block+")";
	}else{
	
	cut = "("+cut_block+") && (B_mass<5.37 && B_mass>5.23) ";
	}
	//RooDataSet mcSet("mcSet","mcSet",*wspace.set("TreeSet"),RooFit::Import(*mc),RooFit::Cut((cut+" && PUweight!=-1 ").c_str()),RooFit::WeightVar("PUweight"));
	RooDataSet dataSet("DataSet","DataSet",*wspace.set("TreeSet"),RooFit::Import(*data),RooFit::Cut(cut.c_str()));
	std::cout << "debug 2 " << dataSet.sumEntries()<< std::endl;
//	wspace.var("mll_ET")->setRange("low",xmin,5.16);
//	wspace.var("mll_ET")->setRange("high",5.42,xmax);

	//std::cout << "histo stdev" << templ->GetStdDev() << std::endl; 
	if (sigPDF == 0){
		// Double CB
		wspace.factory("mean[3.12+00, 3.094+00, 3.17+00]");

		wspace.factory("sigma[4e-2, 3.e-2, 2e-1]");
		wspace.factory("alpha1[1,0.01 ,10]");
		wspace.factory("alpha2[-7,-10, -0]");
		//	wspace.factory("alpha2[-0.5,-1 , 1.0e+2]");
		wspace.factory("n1[10, 0.001, 50]");
		wspace.factory("n2[3, 0.001, 100]");
		//	wspace.factory("n2[1, 0, 10]");
		wspace.factory(("CBShape::CB1("+x+",mean,sigma,alpha1,n1)").c_str());
		wspace.factory(("CBShape::CB2("+x+",mean,sigma,alpha2,n2)").c_str());
	//	wspace.factory(("Gaussian::gaus("+x+",mean,sigma)").c_str());
	}
	if (sigPDF == 1){
		// Voigtian
		wspace.factory("mean[91.1e+00, 90.0e+00, 92e+00]");
		wspace.factory("width[2.4952, 2.4952, 2.4952]");
		//	wspace.factory("mean[3.0969+00, 3.06e+00, 3.10e+00]");
		wspace.factory("sigma[3.477e-02, 1.477e-5,1.7]");
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
		wspace.factory("SUM::sig(CB1,CB2)");
		wspace.factory("SUM::signal(nsig*sig)");
		wspace.factory("SUM::model(nsig*sig,nbkg*bkg)");
		}
		wspace.factory("SUM::background(nbkg*bkg)");


	RooAbsPdf* mod = wspace.pdf("model");
	RooAbsPdf* bkg = wspace.pdf("bkg");
	RooAbsPdf* sig = wspace.pdf("sig");
	RooAbsPdf* signal = wspace.pdf("signal");
	
//	if(isZ) RooFitResult* mc_template = signal->fitTo(mcSet,RooFit::Save(),RooFit::Range(88,96),RooFit::Extended(kTRUE));
//	else RooFitResult* mc_template = signal->fitTo(mcSet,RooFit::Save(),RooFit::Extended(kTRUE));
	resu[0]=wspace.var("mean")->getVal();
	resu[1]=wspace.var("mean")->getError();
	resu[4]=wspace.var("sigma")->getVal();
	resu[5]=wspace.var("sigma")->getError();
	if(isZ){
	wspace.var("width")->setConstant(kTRUE);
	}else{
//	wspace.var("n1")->setConstant(kTRUE);
//	wspace.var("n2")->setConstant(kTRUE);
//	wspace.var("NormC1")->setConstant(kTRUE);
//	wspace.var("NormC2")->setConstant(kTRUE);
//	wspace.var("width")->setConstant(kTRUE);
//	wspace.var("alpha1")->setConstant(kTRUE);
//	wspace.var("alpha2")->setConstant(kTRUE);
	}
	std::cout << " _____________________________________________________________________________________________________DEBUG "  << std::endl;
	if(isZ)RooFitResult* result = mod->fitTo(dataSet,RooFit::Save(),RooFit::Range(88,96),RooFit::Extended(kTRUE));
	else RooFitResult* result = mod->fitTo(dataSet,RooFit::Save(),RooFit::Extended(kTRUE));
	resu[2] = wspace.var("mean")->getVal();
	resu[3] = wspace.var("mean")->getError();
	resu[6]=wspace.var("sigma")->getVal();
	resu[7]=wspace.var("sigma")->getError();
//	wspace.var("sigma")->setVal(wspace.var("sigma")->getVal());
//	wspace.var("sigma")->setConstant(kTRUE);
	
	TCanvas* c1 =new TCanvas("MC_Fit", "fit", 800, 600);
/*	c1->Divide(1,2,0,0);
	c1->cd(1);
	c1->SetGrid();
	gPad->SetLeftMargin(0.10);
	gPad->SetRightMargin(0.05);
	gPad->SetBottomMargin(0.11);
	//gpad->setrightmargin(0.05);
	RooPlot* mll_ETmcframe = wspace.var(x.c_str())->frame();
	mcSet.plotOn(mll_ETmcframe, RooFit::Name("template data"));
	wspace.var("mean")->setVal(resu[0]);
	wspace.var("sigma")->setVal(resu[4]);
	if (isZ)signal->plotOn(mll_ETmcframe,RooFit::Name("signal"),RooFit::Range(70,110),RooFit::LineColor(2),RooFit::MoveToBack()); // this will show fit overlay on canvas
	else signal->plotOn(mll_ETmcframe,RooFit::Name("signal"),RooFit::LineColor(2),RooFit::MoveToBack()); // this will show fit overlay on canvas
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
	c1->cd(2);*/
	c1->SetGrid();
	gPad->SetLeftMargin(0.10);
	gPad->SetRightMargin(0.05);
	RooPlot* mll_ETDataframe = wspace.var(x.c_str())->frame();
	wspace.var("mean")->setVal(resu[2]);
	wspace.var("sigma")->setVal(resu[6]);
	dataSet.plotOn(mll_ETDataframe, RooFit::Name("Dataset"));
	if(isZ)mod->plotOn(mll_ETDataframe,RooFit::Name("model"),RooFit::Range(69,111),RooFit::LineColor(2)); // this will show fit overlay on canvas
	else  mod->plotOn(mll_ETDataframe,RooFit::Name("model"),RooFit::LineColor(2)); // this will show fit overlay on canvas
	if (isZ) mod->plotOn(mll_ETDataframe,RooFit::Name("combinatorial"),RooFit::Components("bkg"),RooFit::DrawOption("FL"),RooFit::Range(69,111),RooFit::LineStyle(kDashed),RooFit::LineColor(kTeal),RooFit::FillColor(kTeal),RooFit::FillStyle(3001));
	else mod->plotOn(mll_ETDataframe,RooFit::Name("combinatorial"),RooFit::Components("bkg"),RooFit::DrawOption("FL"),RooFit::LineStyle(kDashed),RooFit::LineColor(kTeal),RooFit::FillColor(kTeal),RooFit::FillStyle(3001));
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
	


	CMS_lumi(c1,5,1,25);
	if (isZ)c1->SaveAs((PNGPATH+"/"+lable+"/fit_Z"+filename+cut_block+".png").c_str());
	else c1->SaveAs((PNGPATH+"/"+lable+"/fit_Jpsi"+filename+cut_block+".png").c_str());




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
void DoubleSlicer(int bin,float min, float max, bool isZ,std::string lable){

	int i,k;
	float x[2][bin],mean[2][bin],RMS[2][bin];
	float peak_diff[2][bin],peak_diff_sigma[2][bin],pt[2][bin],pt_err[2][bin],sigma_ratio[2][bin],sigma_ratio_unc[2][bin];

	std::string barrel[2];
	barrel[0] = "BarrelCentral";
	barrel[1] = "BarrelEdge";
	double resu[8];
	
	int nPtBin = 4;

//	for(k=0;k<2;k++){
	k= 0;
	for(i=0; i<bin;i++){
		x[k][i]=min +(max-min)/bin*i;
		char range[15]="";
		int n;
		n=sprintf(range,"run>%.3f && run<%.3f",x[k][i],x[k][i]+(max-min)/bin);
		std::string name;
//		if (isZ) name= "PtMin"+std::to_string((int)x[k][i])+"PtMax"+std::to_string((int)(max_edge[i]))+".root";
		if(isZ)name= "mllTree_ZFullStatPtMin40PtMax50BarrelCentral.root";
		else name= "mllTreeFullJPsi.root";
		if (isZ)fit(name,lable,range,1,2,1,i,resu);
		else fit(name,lable,range,0,0,0,i,resu);

		double m_pdg;
		if (isZ) m_pdg = 91.1876;
		else m_pdg = MASS_JPSI;
		peak_diff[k][i]=resu[2];
		
		std::cout << "check: " << peak_diff[i] << std::endl;
		peak_diff_sigma[k][i]= resu[3];
		sigma_ratio[k][i] = resu[6];
		sigma_ratio_unc[k][i] =resu[7];
		  pt[k][i]=x[k][i]+(max-min)/(2*bin);
		 pt_err[k][i] = (max-min)/(2*bin);
		//SavePlot (xaxis, proj[k][i],(PLOTPATH+filename+std::string(range)).c_str() ,false);

	}
//	}

	
	setStyle();
	//std::cout << "check: " << peak_diff[0] << " " << peak_diff[1] << " " <<peak_diff[2] << " " <<peak_diff[3] << " " <<peak_diff[4]<< " " <<peak_diff[5]<< " " << MASS_JPSI << std::endl;
	TGraphErrors* linCentral = new TGraphErrors(bin,pt[0],peak_diff[0],pt_err[0],peak_diff_sigma[0]);
	TGraphErrors* linEdge = new TGraphErrors(bin,pt[1],peak_diff[1],pt_err[1],peak_diff_sigma[1]);
	TH2D* plotter;
	if(isZ)plotter = new TH2D("plotter","plotter",10,320500,325500,10,85,95);
	else plotter = new TH2D("plotter","plotter",10,320500,325500,10,3.05,3.15);
	  TCanvas* canvas = new TCanvas("","",800,600);
	  plotter->GetXaxis()->SetTitle("run");
	  plotter->GetYaxis()->SetTitle("m^{peak}");
	  linCentral->SetMarkerStyle(8);
	  linCentral->SetMarkerSize(1.5);
	  linCentral->SetLineWidth(2);
	  linCentral->SetMarkerColor(kRed-7);
	  linCentral->SetLineColor(kRed-7);
//	  linEdge->SetMarkerStyle(23);
//	  linEdge->SetMarkerSize(2);
//	  linEdge->SetLineWidth(2);
//	  linEdge->SetMarkerColor(kRed-7);
//	  linEdge->SetLineColor(kRed-7);
	  //linearity->GetYaxis()->SetRangeUser(0,1);
	  plotter->Draw();
	  linCentral->Draw("sameP");
	 // linEdge->Draw("sameP");
	 if(isZ) canvas->SaveAs("Run_PeakZ.pdf");
	 else canvas->SaveAs("Run_PeakJPsi.pdf");
	TGraphErrors* sigmaCentral = new TGraphErrors(bin,pt[0],sigma_ratio[0],pt_err[0],sigma_ratio_unc[0]);
	TGraphErrors* sigmaEdge = new TGraphErrors(bin,pt[1],sigma_ratio[1],pt_err[1],sigma_ratio_unc[1]);
	TH2D* plotter1;
	if(isZ)plotter1 = new TH2D("plotter","plotter",10,320500,325500,10,0.6,3);
	else plotter1 = new TH2D("plotter","plotter",10,320500,325500,10,0.6,1.6);
	  TCanvas* canvas_sig = new TCanvas("sigma","sigma",800,600);
	  sigmaCentral->GetHistogram()->SetTitle("run");
	  sigmaCentral->GetHistogram()->SetTitle("#sigma_{data}");
	  sigmaCentral->SetMarkerStyle(8);
	  sigmaCentral->SetLineWidth(2);
	  sigmaCentral->SetMarkerColor(kRed-7);
	  sigmaCentral->SetLineColor(kRed-7);
//	  sigmaEdge->SetMarkerStyle(8);
//	  sigmaEdge->SetLineWidth(2);
//	  sigmaEdge->SetMarkerColor(kRed-7);
//	  sigmaEdge->SetLineColor(kRed-7);
	  //linearity->GetYaxis()->SetRangeUser(0,1);
//	  plotter1->Draw();
	  sigmaCentral->Draw("AP");
	 // sigmaEdge->Draw("sameP");
	 if(isZ) canvas_sig->SaveAs("Run_sigmaZ.pdf");
	 else canvas_sig->SaveAs("Run_sigmaJPsi.pdf");
	  delete plotter; 
	  delete canvas;
	  delete linCentral;

}
