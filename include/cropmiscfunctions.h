#ifndef CROPMISCFUNCTIONS_HH
#define CROPMISCFUNCTIONS_HH
#include "stdio.h"
#include <cmath>
#include "string"
#include <TMath.h>
#include <TString.h>
#include <stdlib.h>
#include <iostream>
#include <TH1.h>
#include <TColor.h>
#include <TStyle.h>
using std::cout;
using std::endl;

inline TString toString(Double_t value){
	char pretty[20];
	sprintf(pretty,"%g",value);
	return TString(pretty);
}

inline TString prettyPrint(Double_t value){
	char pretty[20];
	TString prettyString;
	sprintf (pretty, "%4.3g",value);
	prettyString = pretty;
	return prettyString;
}

inline void binomEff(Double_t num, Double_t denom,Double_t *eff=0, Double_t *err=0){
	if(denom==0){
		*eff = 1.0;
		*err = 0.0;
	}else{
		*eff = num / denom;
		*err = sqrt(fabs(*eff*(1.0-*eff))/denom);
	}
}

inline void quadEff(Double_t num, Double_t numerr, Double_t denom, Double_t denomerr, Double_t *eff=0, Double_t *err=0){
	*eff = num / denom;
	*err = *eff*sqrt((numerr / num)*(numerr / num) + (denomerr / denom)*( denomerr / denom));
}

inline void quadProd(Double_t a, Double_t d_a, Double_t b, Double_t d_b, Double_t *ab, Double_t *d_ab){
	*ab = a*b;
	*d_ab = *ab*sqrt(( d_a/a)*(d_a/a) + (d_b/b)*(d_a/a));
}

inline void quadDiff(Double_t a, Double_t d_a, Double_t b, Double_t d_b, Double_t *aminusb=0, Double_t *d_aminusb=0){
	*aminusb = a - b;
	*d_aminusb = sqrt( d_a*d_a + d_b*d_b);
}

inline void quadSum(Double_t a, Double_t d_a, Double_t b, Double_t d_b, Double_t *aplusb=0, Double_t *d_aplusb=0){
	*aplusb = a + b;
	*d_aplusb = sqrt(d_a*d_a + d_b*d_b);
}

inline void quadPow(Double_t a, Double_t d_a, Double_t tothe, Double_t *atothe=0, Double_t *d_atothe=0){
	*atothe = pow(a,tothe);
	*d_atothe = (*atothe * tothe * d_a)/ a;
}

inline void genericInfo(TString tool){
        TString SimpleToolsVersion = "2.0q";
	cout << tool << ": Part of the SimpleTools Package v" << SimpleToolsVersion << endl;
	cout << "Author: Conor Fitzpatrick conor.fitzpatrick@cern.ch" << endl;
}

inline void weightInfo(){
	cout << "Weightfile syntax: " << endl;
	cout << "<S/B (signal/background)>  <filepath> <ntuplepath> <weight> <preproc. cut> <legend title> <fillstyle> <fill color>" << endl;
	cout << "weight can be a double or a formula, eg a per-event weight" << endl;
	cout << "preproc cut is standard root TCut syntax" << endl;
	cout << "for fillstyle and color see the root TAttFill documentation" << endl;
}
inline void varInfo(){
	cout << "Varfile syntax: " << endl;
	cout << "<var>  <min>   <max>   <bins>  <units> <log/lin>" << endl;
	cout << "Units can be left as an empty pair of quotes (\"\")" << endl;
}

inline void cutInfo(){
	cout << "Cutfile syntax: " << endl;
	cout << "<cut> <min> <max> <res> " << endl;
	cout << "cut must end in either \"<\" or \">\" " << endl;
}

inline void stackerInfo(){
	genericInfo("stacker");
	cout << "Makes pretty stacked or unstacked histogram pages for inclusion in publications" << endl;
	cout << "Usage: " << endl;
	cout << "stacker <Weightfile> <varfile> <plotsY> <plotsX> <resY> <resX> <outputdir> <stackargs>" << endl;
	cout << "plotsY, plotsX define the number of plots on each page in X and Y eg: 2, 2 means 4 plots per page, 2 in X 2 in Y" << endl;
	cout << "resY, resX define the page resolution eg: 2048 1536 means 2048x1536 pixels" << endl;
	cout << "outputdir will contain the plots in svg, pdf, png, eps format as well as outputdir/plots.root containing all plots" << endl;
	cout << "stackargs can be any TH1/THStack draw argument eg: nostack e1 hist etc. Combinations separated by commas" << endl;
	cout << "" << endl;
	weightInfo();

	cout << "" << endl;
	varInfo();
	cout << "" << endl;

}

inline void stackergeninfo(){
	genericInfo("stackergen");
	cout << "Generates var files compatible with stacker from an input weightfile" << endl;
	cout << "Usage: " << endl;
	cout << "stackergen <inputweightfile> <outputvarfile> <bins> <log/lin (x)> <log/lin (y)>" << endl;
	cout << "specify whether the var file should use log x/y axes using \"log\" or \"lin\" twice" << endl;
	cout << "" << endl;
	weightInfo();

	cout << "" << endl;
	varInfo();
	cout << "" << endl;


}

inline void varstocutsInfo(){
	genericInfo("varstocuts");
	cout << "Turns a stacker varfile into a crop cutfile" << endl;
	cout << "Usage: " << endl;
	cout << "varstocuts <inputweightfile> <inputvarfile> <outputcutfile>" << endl;
	cout << "An input weightfile is needed in order to determine if the cut should be > or <" << endl;

	cout << "" << endl;
	weightInfo();

	cout << "" << endl;
	varInfo();
	cout << "" << endl;
	cutInfo();
	cout << "" << endl;

}

inline void mergevarsInfo(){
	genericInfo("mergevars");
	cout << "merges variables of type var1_X var2_X into min(var1_X,var2_X) and max(var1_X,var2_X)" << endl;
	cout << "Usage: " << endl;
	cout << "mergevars <inputvarfile> <outputvarfile> <var1> <var2> {var3}" << endl;
	 cout << "" << endl;
	 varInfo();
	 cout << "" << endl;


}

inline void updatedatastoreInfo(){
	genericInfo("updatedatastore");
	cout << "Rewrites old (simpletooks v1.0j or lower) crop and stacker weightfiles into the new weightfile format" << endl;
		cout << "Usage: " << endl;
	cout << "updatedatastore <inputweightfile> <outputweightfile>" << endl;
		cout << "" << endl;
	weightInfo();
	cout << "" << endl;
}

inline void corrInfo(){
	genericInfo("corr");
	cout << "Creates 2D histograms and correlation matrix for two varfiles" << endl;
	cout << "Usage: " << endl;
	cout << "corr <inputweightfile> <inputvarfile 1> <inputvargile 2> <plotsY> <plotsX> <resY> <resX> <outputdir> <2dargs> <maxcorr>" << endl;
		cout << "To get the full correlation matrix for a set of variables use the same var file as 1 and 2" << endl;
	cout << "Beware though that the number of plots will be the number of vars in 1 X the number of vars in 2" << endl;
	cout << "plotsY, plotsX define the number of plots on each page in X and Y eg: 2, 2 means 4 plots per page, 2 in X 2 in Y" << endl;
	cout << "resY, resX define the page resolution eg: 2048 1536 means 2048x1536 pixels" << endl;
	cout << "outputdir will contain the plots in svg, pdf, png, eps format as well as outputdir/plots.root containing all plots" << endl;
	cout << "2dargs can be any TH2 draw argument eg: COLZ, LEGO etc. Combinations separated by commas" << endl;
	cout << "maxcorr denotes the maximum level of correlation you're willing to tolerate. Anything above this level gets removed from uniquevars.txt in the output folder" << endl;
	cout << "" << endl;
	weightInfo();
	cout << "" << endl;
	varInfo();
	cout << "" << endl;

}

inline void effInfo(){
	genericInfo("eff");
	cout << "Creates binned efficiency distributions from two weightfiles" << endl;
	cout << "Usage: " << endl;
	cout << "eff <inputweightfile numerator> <inputweightfile denominator> <inputvarfile> <plotsY> <plotsX> <resY> <resX> <outputdir> <effargs>" << endl;
	cout << "plotsY, plotsX define the number of plots on each page in X and Y eg: 2, 2 means 4 plots per page, 2 in X 2 in Y" << endl;
	cout << "resY, resX define the page resolution eg: 2048 1536 means 2048x1536 pixels" << endl;
	cout << "outputdir will contain the plots in svg, pdf, png, eps format as well as outputdir/plots.root containing all plots" << endl;
	cout << "effargs can be any TH1 draw argument eg: e1, hist, etc.  Combinations separated by commas" << endl;
	cout << "" << endl;
	weightInfo();
	cout << "" << endl;
	varInfo();
	cout << "" << endl;

}

inline void sepperInfo(){
	genericInfo("sepper");
	cout << "Sorts a varfile by separation power, removing any vars with 0 power" << endl;
	cout << "Usage: " << endl;
	cout << "sepper <inputweightfile> <inputvarfile> <outputvarfile>" << endl;
	cout << "The output varfile can then be passed to varstocuts" << endl;
	cout << "" << endl;
	weightInfo();
	cout << "" << endl;
	varInfo();
	cout << "" << endl;

}


inline void cropBanner(){
	cout << "  ________  ____  ___    "<< endl;
	cout << " / ___/ _ \\/ __ \\/ _ \\   "<<endl;
	cout << "/ /__/ , _/ /_/ / ___/   "<<endl;
	cout << "\\___/_/|_|\\____/_/ "<< endl;
	cout << ""  << endl;
	genericInfo("crop");
}



inline void bwdivInfo(){
	genericInfo("bwdiv");
	cout << "Optimises a set of L0 thresholds" << endl;
	cout << "Usage: " << endl;
	cout << "bwdiv <WeightFile> <cutfile> <max. bkg. eff>" << endl;
	cout << "" << endl;
	weightInfo();
	cout << "" << endl;
	cutInfo();
	cout << "" << endl;

}




inline void cropInfo(){
	cout << "Optimises a set of rectangular cuts" << endl;
	cout << "Usage: " << endl;
	cout << "crop <WeightFile> <cutfile> <steps> <OrderMethod> <output.root> [-b]" << endl;
	cout << "-b specifies batchmode. No plots are produced on-screen, only in output.root" << endl;
	cout << "<steps> are the number of iterations crop should at most attempt to obtain convergency." << endl;
	cout << "OrderMethod can be one of S, B, SB, R, O:" << endl;
	cout << "       S orders based on the cut that is least efficient on signal" << endl;
	cout << "       B orders based on the cut that is most efficient on background" << endl;
	cout << "       SB is the product of the above" << endl;
	cout << "       R orders randomly" << endl;
	cout << "       O uses the original ordering in the cutfile at each step" << endl;
	cout << "" << endl;
	weightInfo();
	cout << "" << endl;
	cutInfo();
	cout << "" << endl;

}



inline Double_t GetSeparation( const TH1D& S, const TH1D& B ){
	// compute "separation" defined as
	// <s2> = (1/2) Int_-oo..+oo { (S^2(x) - B^2(x))/(S(x) + B(x)) dx }
	Double_t separation = 0;
	// sanity checks
	// signal and background histograms must have same number of bins and 
	// same limits
	if ((S.GetNbinsX() != B.GetNbinsX()) || (S.GetNbinsX() <= 0)) {
		cout << "<GetSeparation> signal and background"
			<< " histograms have different number of bins: "
			<< S.GetNbinsX() << " : " << B.GetNbinsX() << endl;
	}

	if (S.GetXaxis()->GetXmin() != B.GetXaxis()->GetXmin() ||
			S.GetXaxis()->GetXmax() != B.GetXaxis()->GetXmax() ||
			S.GetXaxis()->GetXmax() <= S.GetXaxis()->GetXmin()) {
		cout << S.GetXaxis()->GetXmin() << " " << B.GetXaxis()->GetXmin()
			<< " " << S.GetXaxis()->GetXmax() << " " << B.GetXaxis()->GetXmax() 
			<< " " << S.GetXaxis()->GetXmax() << " " << S.GetXaxis()->GetXmin() << endl;
		cout << "<GetSeparation> signal and background"
			<< " histograms have different or invalid dimensions:" << endl;
	}

	Int_t    nbins = S.GetNbinsX(); 
	Double_t nS    = S.Integral( 0, nbins+1, "width" ); // include under/overflow bins
	Double_t nB    = B.Integral( 0, nbins+1, "width" );

	if (nS > 0 && nB > 0) {
		// include under/overflow bins
		for (Int_t ibin=0; ibin<=nbins+1; ibin++) {
			Double_t s = S.GetBinContent( ibin )/nS;
			Double_t b = B.GetBinContent( ibin )/nB;
			// separation

			if (s + b > 0) separation += 0.5*(s - b)*(s - b)/(s + b)*S.GetXaxis()->GetBinWidth(ibin);
			//separation += 0.5*(s - b)*(s - b)/(s + b)*S.GetXaxis()->GetBinWidth(ibin);

			//   cout << nS << "    " << nB << "    " << separation << endl; 
		}
	}
	else {
		cout << "<GetSeparation> histograms with zero entries: "
			<< nS << " : " << nB << " cannot compute separation"
			<< endl;
		separation = 0;
	}
	//cout << separation << endl;
	return separation;
}
class point1D{
	public:
		UInt_t x;
		Double_t y;

};

inline bool pointgt(point1D p1, point1D p2){return p1.y > p2.y;}

inline void set_plot_style()
{
    const Int_t NRGBs = 3;
    const Int_t NCont = 255;
	Double_t stops[NRGBs]  = { 0.00, 0.5, 1.0 };
	Double_t red[NRGBs] = { 0.00, 0.00, 1.00 };
  	Double_t green[NRGBs]   = { 0.00, 1.00, 0.00 };
	Double_t blue[NRGBs] = { 1.00, 0.00, 0.00 };


//    const Int_t NRGBs = 5;
//    const Int_t NCont = 255;
//	Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
// 	Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
//	Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
//	Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };

    
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
#endif
