#ifndef CROPcutspace_HH
#define CROPcutspace_HH
#include <TMath.h>
#include <TString.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TArrow.h>
#include <TGraphErrors.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
using std::cout;
using std::endl;
class cropcutspace {
	public:

		cropcutspace(UInt_t, Double_t, Double_t);
		cropcutspace(Double_t, Double_t, Double_t);
		~cropcutspace();
		void addPoint(UInt_t *, Double_t *, Double_t *, Double_t *, Double_t *, Double_t *, Double_t *,Double_t *,Double_t *,Double_t *,Double_t *);
		UInt_t setMaxFoMStep(bool);

		void plotAll(TString, TCanvas*)const;
		void plotSeff(TString)const;
		void plotBrej(TString)const;
		void plotFoM(TString)const;
		void plotRoC(TString)const;
	inline	Double_t getVal(UInt_t n)const{return Vals[n];}
	inline	Double_t getOptimalVal()const{return Vals[maxFoMStep];}
	inline	UInt_t getOptimalStep() const{return maxFoMStep;}

	inline	Double_t getOptimalSeff() const{return Seff[maxFoMStep];}
	inline	Double_t getOptimalBrej() const{return Brej[maxFoMStep];}
	inline	Double_t getOptimaldSeff() const{return d_Seff[maxFoMStep];}
	inline	Double_t getOptimaldBrej() const{return d_Brej[maxFoMStep];}

	inline	Double_t getMaxFoM() const{return FoM[maxFoMStep];}
	inline	Double_t getdMaxFoM() const{return d_FoM[maxFoMStep];}

	inline	UInt_t getSteps() const {return steps;}


	private:
		void setVals(Double_t);
		Double_t* S;
		Double_t* d_S;
		Double_t* B;
		Double_t* d_B;
		Double_t* Seff;
		Double_t* d_Seff;
		Double_t* Brej;
		Double_t* d_Brej;
		Double_t* FoM;
		Double_t* d_FoM;
		Double_t* Vals;
		Double_t* d_Vals;
		UInt_t maxFoMStep;
		UInt_t steps;
		Double_t res;

};

#endif
