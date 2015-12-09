/* crop: Part of the simpletools package
 * (c) Conor Fitzpatrick, 2008
 *
 * If you find this program useful in whole or in part 
 * please cite this paper: 
 *
 * Feel free to send bugreports, feature requests, patches etc to:
 * conor.fitzpatrick@cern.ch
 *
 */


#include "TFunctor.h"
#include "cropdatastore.h"
#include "cropcutensemble.h"
#include <iostream>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TStyle.h>
#include <TRandom3.h>
//GUI
#include <TApplication.h>
#include <TGFrame.h>
#include<TRootEmbeddedCanvas.h>


//TMVA
#include "TMVA/Interval.h"
//#include "TMVA/LogInterval.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/FitterBase.h"
#include "TMVA/GeneticFitter.h"




//http://www.akiti.ca/mulmatvec2.html
using std::cout;
using std::endl;
using namespace boost;


Double_t ratelimit;
Double_t *maxeffs;

class initfcn : public TMVA::IFitterTarget {
	public: 
		//default constructor
		initfcn(cropdatastore* _dstore, cropdataset* _dset, cropcutensemble * _thresholds){
			dstore = _dstore;
			dset = _dset;
			thresholds = _thresholds;
			datasettotal = dset->getTotalWeightedEntries();
			backgroundtotal = dstore->getTotalBackgroundWeightedEntries();
		}

		double ratelimiter(TString *cut){
			Double_t B = 0;
			Double_t d_B = 0;
			dstore->getWeightedBackgroundEntries(cut,&B,&d_B);
			Double_t eff = B/backgroundtotal;
			if(eff < ratelimit){
				return 0.01*eff;
			}else{
				return (0.01*eff)+10000*(eff-ratelimit)*(eff-ratelimit);
			}
		//		return 10000*(eff-ratelimit)*(eff-ratelimit);
		}

		double EstimatorFunction( std::vector<double> &par){
			TString cut;
			thresholds->buildorcut(par,&cut);
			Double_t S=0;
			Double_t d_S=0;
			dset->getWeightedEntries(&cut,&S,&d_S);
			double f=ratelimiter(&cut);
			f+= ((1.0-(S/datasettotal))*(1.0-(S/datasettotal)));
			return f;
		}

	private:
		cropdatastore* dstore;
		cropdataset* dset;
		cropcutensemble* thresholds;
		Double_t datasettotal, backgroundtotal;
};


class optfcn : public TMVA::IFitterTarget {
	public: 
		//default constructor
		optfcn(cropdatastore* _dstore, cropcutensemble * _thresholds){
			dstore = _dstore;
			thresholds = _thresholds;
			backgroundtotal = dstore->getTotalBackgroundWeightedEntries();
		}

		double ratelimiter(TString *cut){
			Double_t B = 0;
			Double_t d_B = 0;
			dstore->getWeightedBackgroundEntries(cut,&B,&d_B);
			Double_t eff = B/backgroundtotal;
			if(eff < ratelimit){
				return 0.01*eff;
			}else{
				return (0.01*eff)+10000*(eff-ratelimit)*(eff-ratelimit);
			}
		//		return 10000*(eff-ratelimit)*(eff-ratelimit);
		}

		double EstimatorFunction(std::vector<double> &par){
			TString cut;
			thresholds->buildorcut(par,&cut);
			Double_t S=0;
			Double_t d_S=0;
			Double_t Stot = 0;
			double f=ratelimiter(&cut);
			for(UInt_t k = 0; k<dstore->getNSignalDatasets(); k++){
					dstore->getSignalDataset(k)->getWeightedEntries(&cut,&S,&d_S);
					Stot=dstore->getSignalDataset(k)->getTotalWeightedEntries();
				f+=  dstore->getSignalDataset(k)->getSpecial()*(1.0-((S/Stot)/maxeffs[k]))*(1.0-((S/Stot)/maxeffs[k]));
				return f;
			}
		}

	private:
		cropdatastore* dstore;
		cropcutensemble* thresholds;
		Double_t datasettotal, backgroundtotal;
};


int main(int argc, char *argv[]) {

	if(argc != 4){	
		bwdivInfo();
		exit(EXIT_FAILURE);

	}

	std::string weightListName = argv[1];
	std::string thresholdsName = argv[2];
	ratelimit = atof(argv[3]);

	cout << "INFO: Parsing weightfile..." << endl;
	cropdatastore *datastore = new cropdatastore("test",weightListName);
	datastore->initStats();



	cout << "INFO: thresholds..." << endl;
	TRandom3 rndGen(0);
	cropcutensemble* thresholds = new cropcutensemble("thresholds", thresholdsName, &rndGen);
	cout << endl;
	TStopwatch timer;
	timer.Start();

	std::vector<TMVA::Interval*> parameterRanges;
	for(int i=0; i<thresholds->NCutVars; i++){
		parameterRanges.push_back(new TMVA::Interval(thresholds->getCutMinVal(i),thresholds->getCutMaxVal(i),thresholds->getCutResolution(i)));
//		cout << parameterRanges[i]->GetStepSize() << endl;
	}
	maxeffs = new Double_t [datastore->getNSignalDatasets()];
	std::vector<Double_t> pars(parameterRanges.size());



	//TString opt="PopSize=40:Steps=30:Cycles=3:ConvCrit=0.01:SaveBestCycle=5";
	TString opt="PopSize=20:Steps=30:Cycles=3:ConvCrit=0.01:SaveBestCycle=5";
	TMVA::IFitterTarget *ffit; 
	TMVA::FitterBase* fitter;

	TString bestcut, thiscut;
	Double_t S, d_S, Seff, d_Seff;
	Double_t B, d_B, Beff, d_Beff;


	for(UInt_t k = 0; k<datastore->getNSignalDatasets(); k++){
		cropdataset* dataset = datastore->getSignalDataset(k);
		// Now ready for minimization step
		ffit = new initfcn(datastore,dataset,thresholds);
		fitter = new TMVA::GeneticFitter( *ffit,"FitterGA",parameterRanges, opt);
		fitter->Run(pars);
		thresholds->buildorcut(pars, &bestcut);
		cout << bestcut << endl;
		dataset->getEfficiency(new TString(""),&bestcut, &Seff,&d_Seff);
		datastore->getBackgroundEfficiency(new TString(""),&bestcut,&Beff,&d_Beff);
		cout << endl;
		cout << "=========================== Max Eff =========================" << endl;
		cout << dataset->getName() << " MAX EFF: $" << prettyPrint(Seff)<<"\\pm"<<prettyPrint(d_Seff)<< "$ BACKGROUND EFF: $" << prettyPrint(Beff)<<"\\pm"<<prettyPrint(d_Beff) << "$"<< endl;
		maxeffs[k]=Seff;



		cout << "=========================== THRESHOLDS at max eff =========================" << endl;
		for(int i=0; i<thresholds->NCutVars; i++){
			thresholds->getCutVar(i,&bestcut);
			bestcut+=prettyPrint(pars[i]);
			datastore->getBackgroundEfficiency(new TString(""),&bestcut,&Beff,&d_Beff);
			dataset->getEfficiency(new TString(""),&bestcut, &Seff,&d_Seff);
			cout << bestcut << "	$" << prettyPrint(Seff)<<"\\pm"<<prettyPrint(d_Seff) << "$ 	&	$" <<  prettyPrint(Beff)<<"\\pm"<<prettyPrint(d_Beff) << "$" << endl;
		}
		cout << "===========================================================================" << endl;
	}


	ffit = new optfcn(datastore,thresholds);
	fitter = new TMVA::GeneticFitter( *ffit,"FitterGA",parameterRanges, opt);
	fitter->Run(pars);
	thresholds->buildorcut(pars, &bestcut);
	cout << bestcut << endl;


	cout << endl;
	cout << "=========================== THRESHOLDS ====================================" << endl;

	for(int i=0; i<thresholds->NCutVars; i++){
		thresholds->getCutVar(i,&bestcut);
		cout << bestcut << "	" << prettyPrint(pars[i]) << endl;
	}

	cout << "=========================== EFFs ==========================================" << endl;
	cout << endl;

	/*
	   cout << "Thresh.				";
	   for(UInt_t k = 0; k<datastore->getNSignalDatasets(); k++){ 
	   cout << datastore->getSignalDataset(k)->getName() << "	";
	   }
	   cout << endl;
	   for(int i=0; i<thresholds->NCutVars; i++){
	   thresholds->getCutVar(i,&thiscut);
	   thiscut+=prettyPrint(pars[i]);
	   cout << thiscut << "		"; 
	   for(UInt_t k = 0; k<datastore->getNSignalDatasets(); k++){
	   cropdataset* dataset = datastore->getSignalDataset(k);
	   dataset->getEfficiency(new TString(""),&thiscut, &Seff,&d_Seff);
	   cout << prettyPrint(Seff) << "+/-" << prettyPrint(d_Seff)<< "		";
	   }
	   cout << endl;
	   }
	   */

	cout << "Cut	&	";
	for(int i=0; i<thresholds->NCutVars; i++){
		thresholds->getCutVar(i,&thiscut);
		thiscut+=prettyPrint(pars[i]);
		cout << thiscut << "     &       "; 
	}
	cout << "\\\\" << endl;

	for(UInt_t k = 0; k<datastore->getNSignalDatasets(); k++){
		cropdataset* dataset = datastore->getSignalDataset(k);
		cout << datastore->getSignalDataset(k)->getName() << "   &    ";
		for(int i=0; i<thresholds->NCutVars; i++){
			thresholds->getCutVar(i,&thiscut);
			thiscut+=prettyPrint(pars[i]);
			dataset->getEfficiency(new TString(""),&thiscut, &Seff,&d_Seff);
			cout << "$" << prettyPrint(Seff) << "\\pm" << prettyPrint(d_Seff)<< "$      &    ";
		}
		cout << "\\\\" << endl;
	}
	for(UInt_t k = 0; k<datastore->getNBackgroundDatasets(); k++){
		cropdataset* dataset = datastore->getBackgroundDataset(k);
		cout << datastore->getBackgroundDataset(k)->getName() << "   &    ";
		for(int i=0; i<thresholds->NCutVars; i++){
			thresholds->getCutVar(i,&thiscut);
			thiscut+=prettyPrint(pars[i]);
			dataset->getEfficiency(new TString(""),&thiscut, &Seff,&d_Seff);
			cout << "$" << prettyPrint(Seff) << "\\pm" << prettyPrint(d_Seff)<< "$      &    ";
		}
		cout << "\\\\" << endl;
	}

	thresholds->buildorcut(pars, &bestcut);

	//datastore->getWeightedBackgroundEntries(&bestcut,&B,&d_B);
	datastore->getBackgroundEfficiency(new TString(""),&bestcut,&Beff,&d_Beff);
	for(UInt_t k = 0; k<datastore->getNSignalDatasets(); k++){
		cropdataset * dataset = datastore->getSignalDataset(k);
		dataset->getEfficiency(new TString(""), &bestcut,&Seff,&d_Seff);
		cout << dataset->getName() << " OPT EFF: $" << prettyPrint(Seff) << "\\pm" << prettyPrint(d_Seff) <<"$"<< endl;

	}
	cout << "=========================================================================" << endl;
	cout << "BACKGROUND MAX EFF: $" << prettyPrint(Beff) << "\\pm" << prettyPrint(d_Beff) <<"$"<< endl;  

	timer.Stop();
	timer.Print();
	exit(0);
}
