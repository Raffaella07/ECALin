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
	style->SetPadLeftMargin(0.18);//0.16);
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
	style->SetTitleYOffset(1.7); // => 1.15 if exponents
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



TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold
float lumiTextFont   = 42;  // default is helvetic-bold

bool writeExtraText = false;
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
	latex.SetTextSize(lumiTextSize*t*0.01);    
//	latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

//	outOfFrame = false;
	if( outOfFrame )
	{
		latex.SetTextFont(cmsTextFont);
		latex.SetTextAlign(11); 
		latex.SetTextSize(cmsTextSize*t*0.04);    
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
//			latex.DrawLatex(posX_, posY_-posY_*0.49, cmsText);
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
//	outOfFrame = false;

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
int plotter(){


	setStyle();
	double pl_min[2];
	double pl_max[2];

	TGraphErrors* plots[3][2];
	TCanvas* canvas[2];
	TH2D* plotter[2];
	
	int color[3] =  {kOrange+8,kAzure+7,kAzure+7};
	int marker[3] =  {8,8,26};
	std::string Linfiles[3]={"linearity_JPsi_BarrelCentral.txt","linearity_Z_BarrelCentral.txt","linearity_Z_BarrelEdge.txt"};
	std::string Sigfiles[3]={"NormSigma_JPsi_BarrelCentral.txt","NormSigma_Z_BarrelCentral.txt","NormSigma_Z_BarrelEdge.txt"};
	std::string Lentries[3]={"J/#Psi |#eta| < 1","Z |#eta| < 1","Z 1 < |#eta| < 1.5"};
	std::string lable[2]={"(m_{data}^{fit}-m_{MC}^{fit})/m_{PDG}","#frac{#sigma_{data}}{#sigma_{MC}}"};
	int i,j;
	
	pl_min[0] =-0.01;
	pl_min[1] = 0;
	pl_max[0] = 0.01;
	pl_max[1] = 2;
	
	TLegend* l[2];
	
	for(j=0;j<2;j++){

        if (j ==0)l[j] = new TLegend(0.6,0.15,0.9,0.32);
        else l[j] = new TLegend(0.2,0.15,0.5,0.32);
    	l[j]->SetFillStyle(0);
    	l[j]->SetLineStyle(0);
    	l[j]->SetBorderSize(0);
	canvas[j] = new TCanvas(("linearity"+std::to_string(j)).c_str(),("linearity"+std::to_string(j)).c_str(),600,550);
	plotter[j] = new TH2D(("linPlotter"+std::to_string(j)).c_str(),("linPlotter"+std::to_string(j)).c_str(),10,0,80,10,pl_min[j],pl_max[j]);
	
	plotter[j]->GetXaxis()->SetTitle("p_{T} (GeV) ");
	plotter[j]->GetYaxis()->SetTitle(lable[j].c_str());
	plotter[j]->Draw();
		for(i=0;i<3;i++){
		
		if(j==0)plots[i][j]= new TGraphErrors(Linfiles[i].c_str(),"%lg %lg %lg %lg");
		else plots[i][j]= new TGraphErrors(Sigfiles[i].c_str(),"%lg %lg %lg %lg");
		plots[i][j]->SetLineWidth(1);
		plots[i][j]->SetLineColor(color[i]);
		plots[i][j]->SetMarkerColor(color[i]);
		plots[i][j]->SetMarkerStyle(marker[i]);
		plots[i][j]->SetMarkerSize(1.3);
		plots[i][j]->Draw("sameP");
		l[j]->AddEntry(plots[i][j],Lentries[i].c_str(),"lp");

		}	
	l[j]->Draw();
	CMS_lumi(canvas[j],5,0,25);
	if (j==0)canvas[j]->SaveAs("Linearity_full.pdf");
	else canvas[j]->SaveAs("Sigma_full.pdf");


}
	return 1;

}	
