#ifndef CROPDATASET_HH
#define CROPDATASET_HH

#include <TEntryList.h>
#include <stdlib.h>
#include <iostream>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <fstream>
#include "stdio.h"
#include "string"
#include "Riostream.h"
#include <cctype>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <TH1.h>
#include <TH2.h>
#include "cropmiscfunctions.h"
#include "cropvarensemble.h"

using boost::lexical_cast;
using boost::bad_lexical_cast;
using std::cout;
using std::endl;
const char cropdatasetdelimiters[] = " \t\n;";
const TString cropdatasetBstr = "B",cropdatasetSstr = "S",cropdatasetCstr = "#", cropdatasetAnd = "&&", cropdatasetTab ="\t";

class cropvarensemble;

class cropdataset {
	public:
		void init(TString, TString, TString, TString, TString, UInt_t, UInt_t, bool);
		void stackerinit(TString, TString, TString, TString, TString, UInt_t, UInt_t);
		void cropinit(TString,TString, TString, bool, TString, TString);

		cropdataset();
		~cropdataset(){cout << "deleting dataset" << endl;}
		cropdataset(TString, TString, TString, TString, TString, UInt_t, UInt_t, bool);
		cropdataset(TString, TString, TString, TString, TString, UInt_t, UInt_t);
		cropdataset(TString,TString, TString, bool, TString, TString);
		cropdataset(TString, TString, UInt_t);
		void getWeightedEntries(TString *, Double_t *, Double_t *) const;
		void print() const;
		void printOneLine() const;
		void printOneLine(TString *, TString *) const;
		TString getName() const;
		void setName(TString);
		bool isSignal() const{return signal;}
		inline Double_t getTotalEntries() const{return totalEntries;}
		inline 	Double_t getProcEntries() const{return procEntries;}
		inline 	Double_t getProcWeightedEntries() const{return procWeightedEntries;}
		inline 	Double_t getProcWeightedEntriesError() const{return d_procWeightedEntries;}
		inline 	Double_t getTotalWeightedEntries() const{return totalWeightedEntries;}
		inline 	Double_t getTotalWeightedEntriesError() const{return d_totalWeightedEntries;}
		inline UInt_t getColor() const{return color;}
		inline UInt_t getFill() const{return fill;}
		void getEfficiency(TString *, TString *, Double_t *, Double_t *)const;


		TH2D * getHisto(cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const;
		TH2D * getHisto(TString cut, cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const;
		TH2D * getHisto(TString, UInt_t, Double_t, Double_t,TString, UInt_t, Double_t, Double_t)const;
		TH2D * getHisto(TString, TString, UInt_t, Double_t, Double_t,TString, UInt_t, Double_t, Double_t)const;
		TH1D * getHisto(cropvarensemble *vars, UInt_t n)const;
		TH1D * getHisto(TString cut, cropvarensemble *vars, UInt_t n)const;
		TH1D * getHisto(TString, TString, UInt_t, Double_t, Double_t)const;
		TH1D * getHisto(TString, UInt_t, Double_t, Double_t)const;
		inline 	TTree * getTree()const{return procNtuple;}
		inline 	TString getWeightVar()const{return perEventWeightVar;}
		TString toLine()const;
		void setCut(TString *cut);
		void resetCut();
	private:
		TEntryList *tmpelist;
		TEntryList * elist;
		TString name, filePath, ntuplePath, tmpFilePath;
		bool signal;
		TString preprocCutVar, perEventWeightVar;
		TTree *originalNtuple, *procNtuple;
		TFile *originalTFile, *procTFile;
		bool preprocCut;
		UInt_t fill;
		UInt_t color;
		Double_t totalEntries, procEntries;
		Double_t totalWeightedEntries, d_totalWeightedEntries;
		Double_t procWeightedEntries, d_procWeightedEntries;
};

#endif
