#include "../interface/BPark_fitUtils.h"

#define BR_NORES  5.50E-07
#define BR_RES  0.0000618

TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold

bool writeExtraText = false;
TString extraText   = "Preliminary";
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.6;
float lumiTextOffset   = 0.2;
float cmsTextSize      = 0.75;
float cmsTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
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
	latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

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
		if( drawLogo )
		{
			posX_ =   l + 0.045*(1-l-r)*W/H;
			posY_ = 1-t - 0.045*(1-t-b);
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
			if( writeExtraText ) 
			{
				latex.SetTextFont(extraTextFont);
				latex.SetTextAlign(align_);
				latex.SetTextSize(extraTextSize*t);
				latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
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



void fit(TH1D* templ,TH1D* hist,int  sigPDF,int bkgPDF, bool res, int set, double* resu,std::string etaRegion){

	std::cout << "in fit________________________ " << std::endl;
	RooWorkspace wspace("w");
	gStyle->SetOptFit(0000);
	gROOT->SetBatch(true);
	gROOT->SetStyle("Plain");
	gStyle->SetGridStyle(3);
	gStyle->SetOptStat(000000);
	gStyle->SetOptTitle(0);

	std::string PNGPATH = "/eos/home-r/ratramon/www/";
	std::cout << "in fit________3_______________ " << std::endl;
	bool blind=true ;
	double eff_res, eff_nores,eff_misreco;
	double xmin, xmax;
	xmin = hist->GetXaxis()->GetXmin()/*+2*hist->GetXaxis()->GetBinWidth(1)*/;
	xmax = hist->GetXaxis()->GetXmax();//-2*hist->GetXaxis()->GetBinWidth(1);
	std::cout << "bin width " << hist->GetXaxis()->GetBinWidth(1) << std::endl;
	TH1D* mB_res = new TH1D("sig_norm","sig_norm",(int)((xmax-xmin)/hist->GetXaxis()->GetBinWidth(1)), xmin, xmax); 
	TH1D* mB_nores = new TH1D("sig_norm","sig_norm",(int)((xmax-xmin)/hist->GetXaxis()->GetBinWidth(1)), xmin, xmax); 

	std::cout << "limits" << xmin << "   " << xmax << std::endl;

	wspace.factory(("x[5.0,"+std::to_string(xmin)+","+std::to_string(xmax)+"]").c_str());


	wspace.factory("nbkg[100000,0,1000000]");
	double lumi = 10.31;
	int c =set;

	std::cout << templ->GetEntries() << std::endl;
	wspace.factory("nsig[100,0,1000000]");
	wspace.factory("norm[100,20,1000000]");
	wspace.factory("NormC1[100,0,1000000]");
	wspace.factory("x1[5.2,4.5,6]");
	wspace.factory("NormC2[100,0,1000000]");
	wspace.factory("nmisreco[1000,0,1000000]");
	wspace.factory("nores[10,0,1000000]");
	wspace.factory("nVoigt[10,0.00,1000000]");
	wspace.var("x");
	wspace.var("x1");
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
	//RooDoubleCBShape DoubleCB("DoubleCBShape","DoubleCBShape",*wspace.var("x"),*wspace.var("mean"),*wspace.var("sigma"),*wspace.var("alpha1"),*wspace.var("n1"),*wspace.var("alpha2"),*wspace.var("n2"));
	//	wspace.import(DoubleCB);		
	//wspace.var("npred");

	wspace.var("x")->setRange("low",xmin,5.16);
	wspace.var("x")->setRange("high",5.42,xmax);

	std::cout << "histo stdev" << templ->GetStdDev() << std::endl; 
	if (sigPDF == 0){
		// Voigtian
		wspace.factory("width[1.000e-02, 9.000e-03, 1.000e-01]");
		wspace.factory("mean[3.09+00, 3.075+00, 3.1200e+00]");
		wspace.factory("mean1[2.8+00, 2.55e+00, 3.1200e+00]");

		wspace.factory("sigma1[2.2e-1, 1.e-2, 3.01]");
	//	wspace.factory("mean[5.279e+00, 5.27e+00, 5.30e+00]");
	//	wspace.factory(("mean["+std::to_string(templ->GetXaxis()->GetBinCenter(templ->GetMaximumBin()))+","+std::to_string(templ->GetXaxis()->GetBinCenter(templ->GetMaximumBin())-.5*templ->GetStdDev())+", "+std::to_string(templ->GetXaxis()->GetBinCenter(templ->GetMaximumBin())+-5*templ->GetStdDev())+"]").c_str());
	//	wspace.factory(("sigma["+std::to_string(templ->GetStdDev())+","+std::to_string(0.1*templ->GetStdDev())+","+std::to_string(5*templ->GetStdDev())+"]").c_str());
		wspace.factory("sigma[2.2e-2, 1.e-3, 9e-2]");
	wspace.factory("Voigtian::Voigt(x,mean,width,sigma)");
		wspace.factory("Gaussian::gaus(x,mean,sigma)");
		wspace.factory("Gaussian::gaus1(x,mean1,sigma1)");
	}
	if (sigPDF == 1){
		// Gaussian
		wspace.factory("mean[5.2418e+00, 5.20e+00, 5.35e+00]");
		//	wspace.factory("mean[3.0969+00, 3.06e+00, 3.10e+00]");
		wspace.factory("sigma[3.477e-02, 1.477e-02, 7.477e-02]");
			wspace.factory("Gaussian::sig(x,mean,sigma)");
	}
	if (sigPDF == 0){
		// Crystal-ball
		//	wspace.factory("mean[5.279e+00, 5.26e+00, 5.29e+00]");
		//	wspace.factory("mean[3.0969+00, 3.06e+00, 3.10e+00]");
		//	wspace.factory("sigma[3.1858e-02, 1.e-3, 7.e-2]");
		//		wspace.factory("sA[0, -1.e-3, 5.e-2]");
		wspace.factory("alpha1[1,0.001 ,300]");
		wspace.factory("alpha2[-1,-100, -0.001]");
		//	wspace.factory("alpha2[-0.5,-1 , 1.0e+2]");
		wspace.factory("n1[100, 0.001, 200]");
		wspace.factory("n2[3, 0.001, 350]");
		//	wspace.factory("n2[1, 0, 10]");
		wspace.factory("CBShape::CB1(x,mean,sigma,alpha1,n1)");
		wspace.factory("CBShape::CB2(x,mean,sigma,alpha2,n2)");
	}
	if (bkgPDF == 0){
		// Polynomial
		wspace.factory("c0[1, -10.0, 10.0]");
		wspace.factory("c1[-0.01, -10.0, 10.0]");
		wspace.factory("c2[0.4, -10.0, 5.0]");
		wspace.factory("c3[-0.1, -10.0, 10.0]");
		wspace.factory("Polynomial::bkg(x,{c0,c1,c2,c3})");
	}
	if (bkgPDF == 1){
		wspace.factory("c1[0.0, -100.0, 100.0]");
		wspace.factory("Polynomial::bkg(x,{c1})");
	}
	if (bkgPDF == 2){
		// Exponential
		wspace.factory("exp_alpha[-10.0, -100.0, -1]");
		wspace.factory("Exponential::bkg(x,exp_alpha)");
	}
	if (bkgPDF == 3){
		// Polynomial
		wspace.factory("c0[1.0, -1.0, 10.0]");
		wspace.factory("c1[0.1, -5.0, 5.0]");
		wspace.factory("c2[-0.1, -5.0, 5.0]");
		wspace.factory("c3[0.1, -5, 5]");
		wspace.factory("c4[-0.1, -5, 5]");
		wspace.factory("Polynomial::bkg(x,{c0,c1,c2,c3,c4})");
	}

	wspace.defineSet("obs","x");
	//wspace.defineSet("obs1","x1");
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
	//	wspace.factory("Gaussian::signal(x,mean,sigma)");
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

	if (set < 3)wspace.var("sigma")->setRange(templ->GetStdDev()*0.1,templ->GetStdDev());
        else wspace.var("sigma")->setRange(templ->GetStdDev()*0.2,templ->GetStdDev());
	if (set < 3)wspace.var("sigma")->setVal(templ->GetStdDev()*0.6);
        else wspace.var("sigma")->setVal(templ->GetStdDev()*0.6);
	wspace.var("mean")->setRange(templ->GetXaxis()->GetBinCenter(templ->GetMaximumBin())-0.02*templ->GetXaxis()->GetBinCenter(templ->GetMaximumBin()),templ->GetXaxis()->GetBinCenter(templ->GetMaximumBin())+0.02*templ->GetXaxis()->GetBinCenter(templ->GetMaximumBin()));
	wspace.var("mean")->setVal(templ->GetXaxis()->GetBinCenter(templ->GetMaximumBin()));
	
	//RooKeysPdf misreco1("mis","mis", *wspace.var("x"),d,RooKeysPdf::MirrorBoth,1.5);
	//	RooKeysPdf misreco2("misNo","misNo", *wspace.var("x1"),d,RooKeysPdf::NoMirror);
	RooDataHist mc_template = RooDataHist("template fit","template fit",list,templ);
	std::cout << " _____________________________________________________________________________________________________DEBUG "  << std::endl;
	RooFitResult* template_par = signal->fitTo(mc_template,RooFit::Save()/*,RooFit::Extended(kTRUE)*/);
	std::cout << " _____________________________________________________________________________________________________DEBUG "  << std::endl;
	
//	wspace.var("sigma")->setVal(wspace.var("sigma")->getVal());
//	wspace.var("sigma")->setConstant(kTRUE);
	//template_par = signal->fitTo(mc_template,RooFit::Save()/*,RooFit::Extended(kTRUE)*/);
	resu[0]=wspace.var("mean")->getVal();
	resu[1]=wspace.var("mean")->getError();
	resu[4]=wspace.var("sigma")->getVal();
	resu[5]=wspace.var("sigma")->getError();
	
	TCanvas* c1 =new TCanvas("MC_Fit", "fit", 800, 600);
	RooPlot* xmcframe = wspace.var("x")->frame();
	mc_template.plotOn(xmcframe, RooFit::Name("template fit"));
//	wspace.var("sigma1")->setConstant(kTRUE);
//	wspace.var("mean1")->setConstant(kTRUE);
	wspace.var("n1")->setConstant(kTRUE);
	wspace.var("n2")->setConstant(kTRUE);
	wspace.var("NormC1")->setConstant(kTRUE);
	wspace.var("NormC2")->setConstant(kTRUE);
//	wspace.var("width")->setConstant(kTRUE);
	wspace.var("alpha1")->setConstant(kTRUE);
	wspace.var("alpha2")->setConstant(kTRUE);
//wspace.var("alpha2")->setConstant(kTRUE);
		//wspace.var("mean")->setConstant(kTRUE);
//		wspace.var("n1")->setConstant(kTRUE);
	
//		wspace.var("n2")->setConstant(kTRUE);
	signal->plotOn(xmcframe,RooFit::Name("signal"),RooFit::Normalization(wspace.var("nsig")->getVal(),RooAbsReal::NumEvent),RooFit::LineColor(2),RooFit::MoveToBack()); // this will show fit overlay on canvas
	xmcframe->Draw();
        if(res)	c1->SaveAs((PNGPATH+"/plots/fit_MC"+std::to_string((int)set)+etaRegion+"_res.png").c_str());
        else	c1->SaveAs((PNGPATH+"fit_MC"+std::to_string((int)set)+etaRegion+"_nores.png").c_str());
		
//	wspace.var("sigma")->setVal(8e-2);
	double nores;
	std::cout << "______________________________________________check" << wspace.var("nsig")->getVal() << std::endl;

	std::cout << "sigma " << wspace.var("sigma")->getVal() << std::endl;

	RooDataHist data = RooDataHist("data","data",list,hist);
	std::cout << hist->GetEntries() << std::endl;
	//RooAddPdf mod("model", "s+b",RooArgList(*wspace.obj("model")));
	
	//if (set <= 3)wspace.var("sigma")->setRange(templ->GetStdDev()*0.1,templ->GetStdDev()*0.75);
	wspace.var("sigma")->setRange(wspace.var("sigma")->getVal()*0.1,wspace.var("sigma")->getVal()*10);
//	wspace.var("sigma")->setVal(hist->GetStdDev()*0.5);
	wspace.var("mean")->setRange(3.06,3.12);
	wspace.var("mean")->setVal(3.09);
	wspace.var("x")->setRange("peak",2.6,3.5);
	RooFitResult* result;
	std::cout << " ____ "<< wspace.var("mean")->getVal()-3*wspace.var("sigma")->getVal()  << " " <<wspace.var("mean")->getVal()+3*wspace.var("sigma")->getVal() << std::endl;
	if (!res)result = mod->fitTo(data,RooFit::Save(),RooFit::Extended(kTRUE) ,RooFit::Range("peak"));
	else result = mod->fitTo(data,RooFit::Save(),RooFit::Extended(kTRUE));
	std::cout << " _____________________________________________________________________________________________________DEBUG "  << std::endl;
	result->Print(); 
	if(wspace.var("nsig")->getError()>0.8*wspace.var("nsig")->getVal()){

	std::cout << "_______________________________________________________Signal yield error too large, repeating fit with new paramaters" << std::endl; 
	result = mod->fitTo(data,RooFit::Save(),RooFit::Extended(kTRUE)/* ,RooFit::Range("low","high")*/);
	}
	std::cout << "______________________________________________check" << wspace.var("nsig")->getVal() << std::endl;
	result->Print(); 
	std::cout << "sigma " << wspace.var("sigma")->getVal() << std::endl;
	wspace.var("x")->setRange("window",wspace.var("mean")->getVal()-3*wspace.var("sigma")->getVal(),wspace.var("mean")->getVal()+3*wspace.var("sigma")->getVal());
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
	RooPlot* xframe = wspace.var("x")->frame();
	gStyle->SetPalette(kPastel);


	if (!res){
		data.plotOn(xframe, RooFit::Name("data"));
		mod->plotOn(xframe,RooFit::Name("combinatorial"),RooFit::Components("bkg")/*,RooFit::Normalization(wspace.var("nbkg")->getVal(),RooAbsReal::NumEvent)*/,RooFit::DrawOption("L"),RooFit::Range(xmin,xmax),RooFit::LineStyle(kDashed),RooFit::LineColor(kAzure+7));
		mod->plotOn(xframe,RooFit::Name("comb"),RooFit::Components("bkg")/*,RooFit::Normalization(wspace.var("nbkg")->getVal(),RooAbsReal::NumEvent)*/,RooFit::DrawOption("FL"),RooFit::Range(xmin,xmax),RooFit::LineStyle(kDashed),RooFit::LineColor(kAzure+7),RooFit::FillColor(kAzure-3),RooFit::FillStyle(3001));
//		mod->plotOn(xframe,RooFit::Name("comb"),RooFit::Components("bkg")/*,RooFit::Normalization(wspace.var("nbkg")->getVal(),RooAbsReal::NumEvent)*/,RooFit::Range(5.40,xmax),RooFit::DrawOption("FL"),RooFit::LineStyle(kDashed),RooFit::LineColor(kAzure+7),RooFit::FillColor(kAzure-3),RooFit::FillStyle(3001));
		mod->plotOn(xframe,RooFit::Name("model"),RooFit::LineColor(2)/*,RooFit::Normalization(wspace.var("nbkg")->getVal()+wspace.var("nmisreco")->getVal(),RooAbsReal::NumEvent)*/); // this will show fit overlay on canvas
	}else{

		data.plotOn(xframe, RooFit::Name("data"));
		mod->plotOn(xframe,RooFit::Name("model"),RooFit::LineColor(2)); // this will show fit overlay on canvas
		//mod ->plotOn(xframe,RooFit::Name("model"),RooFit::Components("sig"),RooFit::DrawOption("FL"),RooFit::FillColor(9),RooFit::FillStyle(3004),RooFit::LineStyle(6),RooFit::LineColor(9)) ;
		mod->plotOn(xframe,RooFit::Name("combinatorial"),RooFit::Components("bkg"),RooFit::DrawOption("L"),RooFit::FillColor(kMagenta),RooFit::FillStyle(3001),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
//		mod->plotOn(xframe,RooFit::Name("misreco"),RooFit::Components("misr"),RooFit::DrawOption("FL"),RooFit::AddTo("combinatorial",1,1),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3),RooFit::FillColor(kOrange-3),RooFit::FillStyle(3001));
		mod->plotOn(xframe,RooFit::Name("comb"),RooFit::Components("bkg"),RooFit::DrawOption("FL"),RooFit::FillColor(kAzure+7),RooFit::FillStyle(3001),RooFit::LineStyle(kDashed),RooFit::LineColor(kAzure+7));
		mod ->plotOn(xframe,RooFit::Name("Sig"),RooFit::Components("sig"),RooFit::DrawOption("FL"),RooFit::FillColor(kRed-3),RooFit::FillStyle(3004),RooFit::LineColor(kRed-3)) ;
		std::cout << "________________________________________________________" << wspace.var("nsig")->getVal() << std::endl;
	}
	//	mod->plotOn(xframe,RooFit::VisualizeError(result), RooFit::FillColor(kOrange), RooFit::MoveToBack())// # this will show fit overlay on canvas
	//pdf->plotOn(frame,Normalization(0.5,RooAbsReal::Relative)) ;
	
	xframe->GetYaxis()->SetTitleOffset(0.9);
	xframe->GetYaxis()->SetTitleFont(42);
	xframe->GetYaxis()->SetTitleSize(0.05);
	xframe->GetYaxis()->SetLabelSize(0.065);
	xframe->GetYaxis()->SetLabelSize(0.04);
	xframe->GetYaxis()->SetLabelFont(42);
	xframe->GetXaxis()->SetTitleOffset(0.9);
	xframe->GetXaxis()->SetTitleFont(42);
	xframe->GetXaxis()->SetTitleSize(0.05);
	xframe->GetXaxis()->SetLabelSize(0.065);
	xframe->GetXaxis()->SetLabelSize(0.04);
	xframe->GetXaxis()->SetLabelFont(42);

	xframe->GetYaxis()->SetTitle("Events");
	xframe->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV/c^{2}]");
	//xframe.GetXaxis().SetTitle("m(e^{+}e^{-}) [GeV/c^{2}]")
	xframe->SetStats(0);
	xframe->SetMinimum(0);
	xframe->Draw();
	mB_res->SetLineColor(kOrange+8);
	mB_res->SetLineWidth(3);
	mB_res->Draw("samehist");

	TLegend* l;
	l = new TLegend(0.67,0.75,0.95,0.90);
	//else  l = new TLegend(0.67,0.15,0.95,0.50);
	//	l->SetTextFont(72);
	l->SetTextSize(0.03);
	l->AddEntry(xframe->findObject("data"),"Data","lpe");
	l->AddEntry(xframe->findObject("comb"),"Combinatorial bkg","l");
	l->AddEntry(xframe->findObject("Sig"),"Signal","l");

	char v1[10],v2[10],v3[10];
	sprintf(v1,"%.2f",wspace.var("nsig")->getVal());
	sprintf(v2,"%.2f", nbkgWindow);
	sprintf(v3,"%.2f", nores);
/*	if (res) l->AddEntry("" ,("S = "+std::string(v1)).c_str(),"");
	else {
		l->AddEntry("" ,("S_pred = "+std::string(v3)).c_str(),"");
		//	l->AddEntry("" ,("S_fit = "+std::string(v1)).c_str(),"");


	}
	l->AddEntry("",("B = "+std::string(v2)).c_str(),"");
	char v[10];
	if(!res)sprintf(v,"%.2f",nores/sqrt(nores + nbkgWindow));
	if(res)sprintf(v,"%.2f",nSigWindow/sqrt(nSigWindow + nbkgWindow));
	l->AddEntry("" ,("S/#sqrt{S+B} = "+std::string(v)).c_str(),"");*/
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
	delete templ;

}

double PoissonError(RooWorkspace ws){

	std::string PNGPATH = "/eos/home-r/ratramon/www/";
	TH1D* predicted_b = new TH1D("pred","pred",43,4.71,6);
	TH1D* predicted_s = new TH1D("pred_sig","pred_sig",43,4.71,6);
	TH1D* predicted = new TH1D("pred_t","pred_t",43,4.71,6);
	TRandom* rdm = new TRandom();
	int seed;
	seed = time(0);
	rdm->SetSeed(seed);
	int i;
	for (i=1;i<=predicted_s->GetNbinsX();i++){
	ws.var("x")->setVal(predicted_s->GetXaxis()->GetBinCenter(i));
	double mean_b = ws.pdf("bkg")->getVal(ws.set("obs"));
	double mean_s = ws.pdf("sig")->getVal(ws.set("obs"));
	std::cout << "eval" << mean_s <<std::endl;
	predicted_b->SetBinContent(i,mean_b);
	predicted_b->SetBinError(i,sqrt(mean_b));
	predicted_s->SetBinContent(i,mean_s);
	predicted_s->SetBinError(i,sqrt(mean_s));

	}
	
	predicted_b->Scale(ws.var("nbkg")->getVal()/predicted_b->Integral());
	predicted_s->Scale(ws.var("nsig")->getVal()/predicted_s->Integral());

	std::cout << "nbkg" << ws.var("nbkg")->getVal() <<std::endl;
	std::cout << "nsig" << ws.var("nsig")->getVal() <<std::endl;
	for (i=1;i<=predicted->GetNbinsX();i++){
	
	predicted->SetBinContent(i,rdm->Poisson(predicted_b->GetBinContent(i)+predicted_s->GetBinContent(i)));
	predicted->SetBinError(i,sqrt(rdm->Poisson(predicted_b->GetBinContent(i)+predicted_s->GetBinContent(i)))); // check for right formulation of errors 

	}
	
	std::cout << "nsig" << ws.var("nsig")->getVal() <<std::endl;
/*	TCanvas* c = new TCanvas("poisson","poisson",800,600);
	predicted->SetMarkerStyle(8);
	predicted->GetXaxis()->SetTitle("m_{B}(GeV)");
	predicted->Draw();*/
	RooArgList list = RooArgList(*ws.set("obs"));
	RooAbsPdf* signal= ws.pdf("signal");
	RooAbsPdf* mod = ws.pdf("model");
	RooAbsPdf* bkg = ws.pdf("bkg");
	RooAbsPdf* background = ws.pdf("background");
	RooAbsPdf* sig = ws.pdf("sig");

	RooDataHist data = RooDataHist("data","data",list,predicted);
	ws.var("nsig")->setConstant(kFALSE);
	RooFitResult* result = mod->fitTo(data,RooFit::Save(),RooFit::Extended(kTRUE)/* ,RooFit::Range("low","high")*/);
	
	TCanvas* c2 =new TCanvas("fig_binnedFit", "fit", 800, 600);
	c2->SetGrid();
	c2->cd();
	gPad->SetLeftMargin(0.10);
	gPad->SetRightMargin(0.05);
	RooPlot* xframe = ws.var("x")->frame();
		data.plotOn(xframe, RooFit::Name("data"));
		mod->plotOn(xframe,RooFit::Name("model"),RooFit::LineColor(2)); // this will show fit overlay on canvas
		//mod ->plotOn(xframe,RooFit::Name("model"),RooFit::Components("sig"),RooFit::DrawOption("FL"),RooFit::FillColor(9),RooFit::FillStyle(3004),RooFit::LineStyle(6),RooFit::LineColor(9)) ;
		mod->plotOn(xframe,RooFit::Name("combinatorial"),RooFit::Components("bkg"),RooFit::DrawOption("L"),RooFit::FillColor(kMagenta),RooFit::FillStyle(3001),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
		mod->plotOn(xframe,RooFit::Name("misreco"),RooFit::Components("misr"),RooFit::DrawOption("FL"),RooFit::AddTo("combinatorial",1,1),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange-3),RooFit::FillColor(kOrange-3),RooFit::FillStyle(3001));
		mod->plotOn(xframe,RooFit::Name("comb"),RooFit::Components("bkg"),RooFit::DrawOption("FL"),RooFit::FillColor(kAzure+7),RooFit::FillStyle(3001),RooFit::LineStyle(kDashed),RooFit::LineColor(kAzure+7));
		mod ->plotOn(xframe,RooFit::Name("Sig"),RooFit::Components("sig"),RooFit::DrawOption("FL"),RooFit::FillColor(kRed-3),RooFit::FillStyle(3004),RooFit::LineColor(kRed-3)) ;
	//	std::cout << "________________________________________________________" << wspace.var("nsig")->getVal() << std::endl;
	xframe->GetYaxis()->SetTitleOffset(0.9);
	xframe->GetYaxis()->SetTitleFont(42);
	xframe->GetYaxis()->SetTitleSize(0.05);
	xframe->GetYaxis()->SetLabelSize(0.065);
	xframe->GetYaxis()->SetLabelSize(0.04);
	xframe->GetYaxis()->SetLabelFont(42);
	xframe->GetXaxis()->SetTitleOffset(0.9);
	xframe->GetXaxis()->SetTitleFont(42);
	xframe->GetXaxis()->SetTitleSize(0.05);
	xframe->GetXaxis()->SetLabelSize(0.065);
	xframe->GetXaxis()->SetLabelSize(0.04);
	xframe->GetXaxis()->SetLabelFont(42);

	xframe->GetYaxis()->SetTitle("Events");
	xframe->GetXaxis()->SetTitle("m(K^{+}e^{+}e^{-}) [GeV/c^{2}]");
	//xframe.GetXaxis().SetTitle("m(e^{+}e^{-}) [GeV/c^{2}]")
	xframe->SetStats(0);
	xframe->SetMinimum(0);
	xframe->Draw();

	TLegend* l;
	l = new TLegend(0.67,0.55,0.95,0.90);
	//else  l = new TLegend(0.67,0.15,0.95,0.50);
	//	l->SetTextFont(72);
	l->SetTextSize(0.03);
		l->AddEntry(xframe->findObject("data"),"Data","lpe");
//	else	l->AddEntry(xframe->findObject("data_blind"),"Data","lpe");
	l->AddEntry(xframe->findObject("comb"),"Combinatorial bkg","l");
	l->AddEntry(xframe->findObject("misreco"),"Misreco bkg ","l");
	l->AddEntry(xframe->findObject("Sig"),"Signal","l");
/*
	char v1[10],v2[10],v3[10];
	sprintf(v1,"%.2f",wspace.var("nsig")->getVal());
	sprintf(v2,"%.2f", nbkgWindow);
	sprintf(v3,"%.2f", nores);
	if (res) l->AddEntry("" ,("S = "+std::string(v1)).c_str(),"");
	else {
		l->AddEntry("" ,("S_pred = "+std::string(v3)).c_str(),"");
		//	l->AddEntry("" ,("S_fit = "+std::string(v1)).c_str(),"");


	}
//	l->AddEntry("",("B = "+std::string(v2)).c_str(),"");
//	char v[10];
	if(!res)sprintf(v,"%.2f",nores/sqrt(nores + nbkgWindow));
	if(res)sprintf(v,"%.2f",nSigWindow/sqrt(nSigWindow + nbkgWindow));
	l->AddEntry("" ,("S/#sqrt{S+B} = "+std::string(v)).c_str(),"");*/
	l->Draw();
	CMS_lumi(c2,5,0,37);
	c2->SaveAs((PNGPATH+"/newElectronPF/Poisson.png").c_str());
	return 34.;
}



