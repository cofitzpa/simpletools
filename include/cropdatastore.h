#ifndef CROPDATASTORE_HH
#define CROPDATASTORE_HH

#include "cropdataset.h"
#include <TH1.h>
#include <TH2.h>
#include "THStack.h"
#include "TPaveStats.h"
#include <stdlib.h>
#include "cropvarensemble.h"
using std::cout;
using std::endl;
using std::vector;
class cropdataset;
class cropvarensemble;
class cropdatastore {
	public:
		cropdatastore();
		cropdatastore(TString);
		cropdatastore(TString, TString);
		~cropdatastore();
		void addDataset(cropdataset*);
		void print() const;
		void print(TString, TString) const;
		void printOneLine(bool) const;
		void printOneLine(bool, TString *, TString *) const;
		cropdataset* getSignalDataset(UInt_t) const;
		cropdataset* getBackgroundDataset(UInt_t) const;
		cropdataset* getDataset(TString) const;
		inline UInt_t getNSignalDatasets() const{return NSignalDatasets;}
		inline UInt_t getNBackgroundDatasets() const{return NBackgroundDatasets;}
		void getWeightedSignalEntries(TString *, Double_t *, Double_t *) const;
		void getWeightedBackgroundEntries(TString *, Double_t *, Double_t *) const;
		inline Double_t getTotalSignalEntries() const{return totalSignalEntries;}
		inline Double_t getProcSignalEntries() const{return procSignalEntries;}
		inline Double_t getProcSignalWeightedEntries() const{return procSignalWeightedEntries;}
		inline Double_t getProcSignalWeightedEntriesError() const{return d_procSignalWeightedEntries;}
		inline Double_t getTotalSignalWeightedEntries() const{return totalSignalWeightedEntries;}
		inline Double_t getTotalSignalWeightedEntriesError() const{return d_totalSignalWeightedEntries;}
		inline Double_t getTotalBackgroundEntries() const{return totalBackgroundEntries;}
		inline Double_t getProcBackgroundEntries() const{return procBackgroundEntries;}
		inline Double_t getProcBackgroundWeightedEntries() const{return procBackgroundWeightedEntries;}
		inline Double_t getProcBackgroundWeightedEntriesError() const{return d_procBackgroundWeightedEntries;}
		inline Double_t getTotalBackgroundWeightedEntries() const{return totalBackgroundWeightedEntries;}
		inline Double_t getTotalBackgroundWeightedEntriesError() const{return d_totalBackgroundWeightedEntries;}
		void getSignalEfficiency(TString *, TString *, Double_t*,Double_t*)const;
		void getBackgroundEfficiency(TString *, TString *, Double_t*,Double_t*)const;
		void setName(TString _name){name = _name;}
		TString getName() const {return name;}

		TH2D * getSignalHisto(cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const;
		TH2D * getSignalHisto(TString cut, cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const;
		TH2D * getSignalHisto(TString, UInt_t, Double_t, Double_t,TString, UInt_t, Double_t, Double_t)const;
		TH2D * getSignalHisto(TString, TString, UInt_t, Double_t, Double_t,TString, UInt_t, Double_t, Double_t)const;
		TH2D * getBackgroundHisto(cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const;
		TH2D * getBackgroundHisto(TString cut, cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const;
		TH2D * getBackgroundHisto(TString, UInt_t, Double_t, Double_t,TString, UInt_t, Double_t, Double_t)const;
		TH2D * getBackgroundHisto(TString, TString, UInt_t, Double_t, Double_t,TString, UInt_t, Double_t, Double_t)const;
		TH2D * getHisto(cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const;
		TH2D * getHisto(TString cut, cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const;
		TH2D * getHisto(TString, UInt_t, Double_t, Double_t,TString, UInt_t, Double_t, Double_t)const;
		TH2D * getHisto(TString, TString, UInt_t, Double_t, Double_t,TString, UInt_t, Double_t, Double_t)const;

		TH1D* getHisto(cropvarensemble *_data, UInt_t n)const;
		TH1D* getHisto(TString cut, cropvarensemble *_data, UInt_t n)const;
		TH1D* getBackgroundHisto(cropvarensemble *_data, UInt_t n)const;
		TH1D* getBackgroundHisto(TString cut, cropvarensemble *_data, UInt_t n)const;
		TH1D* getSignalHisto(cropvarensemble *_data, UInt_t n)const;
		TH1D* getSignalHisto(TString cut, cropvarensemble *_data, UInt_t n)const;
		TH1D* getSignalHisto(TString varString, UInt_t bins, Double_t min, Double_t max)const;
		TH1D* getBackgroundHisto(TString varString, UInt_t bins, Double_t min, Double_t max)const;
		TH1D* getSignalHisto(TString cut, TString varString, UInt_t bins, Double_t min, Double_t max)const;
		TH1D* getBackgroundHisto(TString cut, TString varString, UInt_t bins, Double_t min, Double_t max)const;
		TH1D* getHisto(TString cut, TString varString, UInt_t bins, Double_t min, Double_t max)const;
		TH1D* getHisto(TString varString, UInt_t bins, Double_t min, Double_t max)const;

		THStack * getStack(cropvarensemble *_data, UInt_t n)const;
		THStack * getStack(TString cut, cropvarensemble *_data, UInt_t n)const;
		THStack * getStack(TString cut, TString varString, UInt_t bins, Double_t min, Double_t max)const;
		THStack * getStack(TString varString, UInt_t bins, Double_t min, Double_t max)const;
		void initStats() const;
		void finalStats(TString *, TString *) const;
		void writeToFile(TString)const;
		void setCut(TString *);
		void resetCut();

	private:
		TString name;
		UInt_t NSignalDatasets, NBackgroundDatasets;
		vector<cropdataset*> SignalDatasets, BackgroundDatasets;
		Double_t totalSignalEntries, procSignalEntries;
		Double_t totalBackgroundEntries, procBackgroundEntries;
		Double_t totalSignalWeightedEntries, d_totalSignalWeightedEntries;
		Double_t totalBackgroundWeightedEntries, d_totalBackgroundWeightedEntries;
		Double_t procSignalWeightedEntries, d_procSignalWeightedEntries;
		Double_t procBackgroundWeightedEntries, d_procBackgroundWeightedEntries;

};
#endif
