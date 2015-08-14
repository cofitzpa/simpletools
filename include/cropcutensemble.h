#ifndef CROPcutensemble_HH
#define CROPcutensemble_HH
#include "cropmiscfunctions.h"

#include <stdlib.h>
#include <iostream>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TNtuple.h>
#include <fstream>
#include "stdio.h"
#include "string"
#include "Riostream.h"
#include <cctype>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <TH1.h>
#include "cropmiscfunctions.h"
#include "cropcutspace.h"
#include "cropdatastore.h"
#include "cropvarensemble.h"
#include "TRandom3.h"

using boost::lexical_cast;
using boost::bad_lexical_cast;
using std::cout;
using std::endl;
using std::vector;
const char cropcutensembledelimiters[] = " \t\n;";
const TString cropcutensembleCstr = "#", cropcutensembleAnd = "&&", cropcutensembleGreaterThan = ">",  cropcutensembleLessThan = "<", cropcutensembleTab = "\t";

class cropvarensemble;
class cropdatastore;
class cropcutensemble {
	public:
		UInt_t NCutVars;
		cropcutensemble();
		~cropcutensemble(){cout << "deleting cutensemble" << endl;}
		cropcutensemble(TString, TRandom3 *);
		cropcutensemble(TString, TString, TRandom3 *);
		cropcutensemble(cropvarensemble *, cropdatastore *);
		inline	void addCutSpacePoint(UInt_t n, UInt_t step, Double_t * _s, Double_t * _d_s, Double_t * _b, Double_t * _d_b, Double_t * _seff, Double_t * _d_seff,Double_t * _brej,Double_t * _d_brej,Double_t *_fom,Double_t * _d_fom){CutSpaces[CutOrder[n]].addPoint(&step,_s,_d_s,_b,_d_b,_seff,_d_seff,_brej,_d_brej,_fom,_d_fom);}
		void addCut(TString,TString, UInt_t);
		void addCut(TString, Double_t, Double_t, Double_t);
		void addCut(TString, Double_t, Double_t, UInt_t);
		void print() const;
		void printOneLine(UInt_t) const;
		bool OrderAscendingSeff();
		bool OrderDescendingBrej();
		bool OrderAscendingSeffBeff();
		bool OrderRandom();
		bool OrderOriginal();
		inline	TString getName() const{return name;}
		void setName(TString _name){name = _name;}
		inline	void getCutVar(UInt_t n, TString *cutOut) const{*cutOut = CutVars[CutOrder[n]];}
		inline	void getCut(UInt_t n, Double_t val, TString *cutOut) const{
			*cutOut = CutVars[CutOrder[n]];
			cutOut->Append(toString(val));
		}
		inline	void getCut(UInt_t n, UInt_t s, TString *cutOut) const{
			*cutOut = CutVars[CutOrder[n]];
			cutOut->Append(toString(CutSpaces[CutOrder[n]].getVal(s)));
		}
		inline	void getOptimalCut(UInt_t n, TString *cutOut) const{
			*cutOut = CutVars[CutOrder[n]];
			cutOut->Append(toString(CutSpaces[CutOrder[n]].getOptimalVal()));
		}
		inline	Double_t getCutMinVal(UInt_t n) const{return CutMinVals[CutOrder[n]];}
		inline	Double_t getCutMaxVal(UInt_t n) const{return CutMaxVals[CutOrder[n]];}
		inline	Double_t getCutResolution(UInt_t n)const{return CutResolutions[CutOrder[n]];}
		inline	UInt_t getCutSteps(UInt_t n)const{return CutSpaces[CutOrder[n]].getSteps();}
		inline	UInt_t getOptimalStep(UInt_t n)const{return CutSpaces[CutOrder[n]].getOptimalStep();}
		inline	UInt_t setOptimalStep(UInt_t n){return CutSpaces[CutOrder[n]].setMaxFoMStep(isAGreaterThanCut(CutOrder[n]));}
		void getOptimalEnsemble(UInt_t n, TString *)const;
		inline	Double_t getMaxFoM(UInt_t n)const{return CutSpaces[CutOrder[n]].getMaxFoM();}
		inline	Double_t getdMaxFoM(UInt_t n)const{return CutSpaces[CutOrder[n]].getdMaxFoM();}
		bool isAGreaterThanCut(UInt_t) const;
		inline	void plotCut(UInt_t n, TCanvas * c){CutSpaces[CutOrder[n]].plotAll(CutVars[CutOrder[n]],c);}
		void writeAllPlots(TCanvas * c);
		TString toLine(UInt_t) const;
		void writeToFile(TString) const;

	private:
		TRandom3 *rng;
		TString name;
		vector<TString> CutVars;
		vector<Double_t> CutMinVals;
		vector<Double_t> CutMaxVals;
		vector<Double_t> CutResolutions;
		vector<cropcutspace> CutSpaces;
		vector<UInt_t> CutOrder;
};

#endif
