#ifndef cropvarensemble_HH
#define cropvarensemble_HH
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
#include "TRandom3.h"
#include <float.h>
#include <vector>
#include <algorithm>
#include "cropdatastore.h"
#include "cropcutensemble.h"
#include "cropdataset.h"
#include "TPaveStats.h"
#include "THStack.h"
using boost::lexical_cast;
using boost::bad_lexical_cast;
using std::cout;
using std::endl;
using std::vector;
const char cropvarensembledelimiters[] = " \t\n;";
const TString cropvarensembleCstr = "#", cropvarensembleAnd = "&&", cropvarensembleTab = "\t";

inline bool tstringcmp(TString c1, TString c2){return c1 == c2;}
inline bool tstringlt(TString c1, TString c2){return c1.Hash() < c2.Hash();}

class cropcutensemble;
class cropdataset;
class cropdatastore;
class cropvarensemble {
	public:
		UInt_t Nvars;
		cropvarensemble();
		~cropvarensemble();
		cropvarensemble(TString);	
		cropvarensemble(TString, TString);

		cropvarensemble(cropcutensemble *, bool);
		cropvarensemble(cropdataset *, bool, bool, UInt_t);
		cropvarensemble(cropdatastore *, bool, bool, UInt_t);
		cropvarensemble(cropdatastore *, TString, bool, bool, UInt_t);
		vector <TString> VarsFromDataSet(cropdataset *);
		vector <TString> VarsFromList(TString);
		void addVar(TString,TString, UInt_t);
		void addVar(TString, TString, Double_t, Double_t, Double_t, bool, bool);
		void addVar(TString, TString, Double_t, Double_t, UInt_t, bool, bool);
		void print() const;
		void printOneLine(UInt_t) const;
		inline TString getName() const{return name;}
		inline void setName(TString _name){name = _name;}
		inline TString getVar(UInt_t n)const{return vars[varOrder[n]];}
		//inline	void getVar(UInt_t n, TString *varOut) const{*varOut = vars[varOrder[n]];}
		inline	Double_t getVarMinVal(UInt_t n) const{return mins[varOrder[n]];}
		inline	Double_t getVarMaxVal(UInt_t n) const{return maxs[varOrder[n]];}
		inline	Double_t getVarResolution(UInt_t n)const{return ress[varOrder[n]];}
		inline UInt_t getVarBins(UInt_t n)const{return bins[varOrder[n]];}
		inline void getUnits(UInt_t n, TString *unit)const{*unit = units[varOrder[n]];}
		inline bool isLogScale(UInt_t n)const{return logScale[varOrder[n]];}
		inline bool isUseless(UInt_t n)const{return useless[varOrder[n]];}
		inline void makeUseless(UInt_t n){useless[n] = true;}
		TString toLine(UInt_t n)const;
		void writeToFile(TString) const;
		THStack *getStack(UInt_t n, cropdatastore *data);
		TH1D *getHisto(UInt_t n, cropdatastore *data);
		TH1D *getSignalHisto(UInt_t n, cropdatastore *data);
		TH1D *getBackgroundHisto(UInt_t n, cropdatastore *data);
		UInt_t findVar(TString * _var)const;
		void sortbySepPower(cropdatastore *_data);


		void mergeVars(TString v1, TString v2, TString v3, TString v4, TString v5, TString v6);
		void mergeVars(TString v1, TString v2, TString v3, TString v4, TString v5);
		void mergeVars(TString v1, TString v2, TString v3, TString v4);
		void mergeVars(TString v1, TString v2, TString v3);
		void mergeVars(TString v1, TString v2);

	private:
		TString name;
		vector<TString> units;
		vector<TString> vars;
		vector<Double_t> mins;
		vector<Double_t> maxs;
		vector<Double_t> ress;
		vector<UInt_t> bins;
		vector<UInt_t> varOrder;
		vector<bool> logScale;
		vector<bool> useless;
};

#endif
