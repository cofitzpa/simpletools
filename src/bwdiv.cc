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
Double_t *maxeffs, *dmaxeffs, *finaleffs, *dfinaleffs; 



std::vector< std::vector<Double_t> > maxincleffs; 
std::vector< std::vector<Double_t> > dmaxincleffs;
std::vector< std::vector<Double_t> > finalincleffs;
std::vector< std::vector<Double_t> > dfinalincleffs;

class initfcn : public TMVA::IFitterTarget {
	public: 
		//default constructor
		initfcn(cropdatastore* _dstore, cropdataset* _dset, cropcutensemble * _thresholds){
			dstore = _dstore;
			dset = _dset;
			thresholds = _thresholds;
		}

		double ratelimiter(TString *cut){
			Double_t Beff = 0;
			Double_t d_Beff = 0;
			dstore->getBackgroundEfficiency(new TString(""),cut, &Beff,&d_Beff);
			if(Beff < ratelimit){
				return 0.01*Beff;
			}else{
				return (0.01*Beff)+10000*(Beff-ratelimit)*(Beff-ratelimit);
			}
		}

		double EstimatorFunction( std::vector<double> &par){
			TString cut;
			thresholds->buildorcut(par,&cut);
			Double_t Seff=0;
			Double_t d_Seff=0;
			dset->getEfficiency(new TString(""),&cut, &Seff,&d_Seff);
			double f=ratelimiter(&cut);
			f+= ((1.0-Seff)*(1.0-Seff));
			return f;
		}

	private:
		cropdatastore* dstore;
		cropdataset* dset;
		cropcutensemble* thresholds;
};


class optfcn : public TMVA::IFitterTarget {
	public: 
		//default constructor
		optfcn(cropdatastore* _dstore, cropcutensemble * _thresholds){
			dstore = _dstore;
			thresholds = _thresholds;
		}

		double ratelimiter(TString *cut){
			Double_t Beff = 0;
			Double_t d_Beff = 0;
			dstore->getBackgroundEfficiency(new TString(""),cut, &Beff,&d_Beff);
			if(Beff < ratelimit){
				return 0.01*Beff;
			}else{
				return (0.01*Beff)+10000*((Beff-ratelimit)*(Beff-ratelimit));
			}
		}

		double EstimatorFunction(std::vector<double> &par){
			TString cut;
			thresholds->buildorcut(par,&cut);
			Double_t Seff=0;
			Double_t d_Seff=0;
			Double_t Stot = 0;
			double f=ratelimiter(&cut);
			for(UInt_t k = 0; k<dstore->getNSignalDatasets(); k++){
				dstore->getSignalDataset(k)->getEfficiency(new TString(""),&cut, &Seff,&d_Seff);
				f+=  dstore->getSignalDataset(k)->getSpecial()*((1.0-(Seff/maxeffs[k]))*(1.0-(Seff/maxeffs[k])));
			}
				return f;
		}

	private:
		cropdatastore* dstore;
		cropcutensemble* thresholds;
};


int main(int argc, char *argv[]) {

	if(argc != 4){
		if(argc != 5){
			bwdivInfo();
			exit(EXIT_FAILURE);
		}
	}

	std::string weightListName = argv[1];
	std::string thresholdsName = argv[2];
	ratelimit = atof(argv[3]);
	TString outname = "output.C";
	if(argc==5){
		outname = argv[4];
	}

	TCanvas* c=new TCanvas("Effs","Effs",1024,768);
	c->SetBottomMargin(0.5);

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
	dmaxeffs = new Double_t [datastore->getNSignalDatasets()];
	finaleffs = new Double_t [datastore->getNSignalDatasets()];
	dfinaleffs = new Double_t [datastore->getNSignalDatasets()];



	maxincleffs.resize(datastore->getNSignalDatasets());
	dmaxincleffs.resize(datastore->getNSignalDatasets());
	finalincleffs.resize(datastore->getNSignalDatasets());
	dfinalincleffs.resize(datastore->getNSignalDatasets());

	std::vector<Double_t> pars(parameterRanges.size());



	//TString opt="PopSize=40:Steps=30:Cycles=3:ConvCrit=0.01:SaveBestCycle=5";
	TString opt="PopSize=20:Steps=30:Cycles=3:ConvCrit=0.01:SaveBestCycle=5";
	TMVA::IFitterTarget *ffit; 
	TMVA::FitterBase* fitter;

	TString bestcut, thiscut;
	Double_t S, d_S, Seff, d_Seff;
	Double_t B, d_B, Beff, d_Beff;


	for(UInt_t k = 0; k<datastore->getNSignalDatasets(); k++){
		maxincleffs[k].resize(thresholds->NCutVars);
		dmaxincleffs[k].resize(thresholds->NCutVars);
		finalincleffs[k].resize(thresholds->NCutVars);
		dfinalincleffs[k].resize(thresholds->NCutVars);
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
		cout << "CHI^2: " << ffit->EstimatorFunction(pars) << endl;
		cout << dataset->getName() << " MAX EFF: $" << prettyPrint(Seff)<<"\\pm"<<prettyPrint(d_Seff)<< "$ BACKGROUND EFF: $" << prettyPrint(Beff)<<"\\pm"<<prettyPrint(d_Beff) << "$"<< endl;

		maxeffs[k]=Seff;
		dmaxeffs[k]=d_Seff;


		cout << "=========================== THRESHOLDS at max eff =========================" << endl;
		for(int i=0; i<thresholds->NCutVars; i++){
			thresholds->getCutVar(i,&bestcut);
			bestcut+=prettyPrint(pars[i]);
			datastore->getBackgroundEfficiency(new TString(""),&bestcut,&Beff,&d_Beff);
			dataset->getEfficiency(new TString(""),&bestcut, &Seff,&d_Seff);

			maxincleffs[k][i] = Seff;
			dmaxincleffs[k][i] = d_Seff;
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
	cout << "============================== DONE =======================================" << endl;
	cout << "FINAL CHI^2: " << ffit->EstimatorFunction(pars) << endl;
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

		TH1D* inclmaxeffbars = new TH1D(dataset->getName()+"inclmaxEffbars","Efficiencies per line for "+dataset->getName(),thresholds->NCutVars,0,thresholds->NCutVars);
		//inclmaxeffbars->SetXTitle("Line");
		inclmaxeffbars->SetYTitle("Efficiency");
		inclmaxeffbars->SetStats(kFALSE);
		inclmaxeffbars->SetFillStyle(3345);

		inclmaxeffbars->SetLineWidth(2);
		inclmaxeffbars->SetLineColor(kRed);
		inclmaxeffbars->SetFillColor(kRed);

		TH1D* inclfinaleffbars = new TH1D(dataset->getName()+"inclfinalEffbars","Efficiencies per line for "+dataset->getName(),thresholds->NCutVars,0,thresholds->NCutVars);
		//inclfinaleffbars->SetXTitle("Line");
		inclfinaleffbars->SetYTitle("Efficiency");
		inclfinaleffbars->SetStats(kFALSE);
		inclfinaleffbars->SetFillStyle(3354);

		inclfinaleffbars->SetLineWidth(2);
		inclfinaleffbars->SetLineColor(kBlue);
		inclfinaleffbars->SetFillColor(kBlue);

		cout << dataset->getName() << "   &    ";
		for(int i=0; i<thresholds->NCutVars; i++){
			thresholds->getCutVar(i,&thiscut);
			thiscut+=prettyPrint(pars[i]);
			dataset->getEfficiency(new TString(""),&thiscut, &Seff,&d_Seff);
			cout << "$" << prettyPrint(Seff) << "\\pm" << prettyPrint(d_Seff)<< "$      &    ";

			finalincleffs[k][i] = Seff;
			dfinalincleffs[k][i] = d_Seff;



		}
		cout << "\\\\" << endl;

		for(int i=1; i<=thresholds->NCutVars; i++){
			inclmaxeffbars->SetBinContent(i,maxincleffs[k][i-1]);
			inclmaxeffbars->SetBinError(i,dmaxincleffs[k][i-1]);
			inclfinaleffbars->SetBinContent(i,finalincleffs[k][i-1]);
			inclfinaleffbars->SetBinError(i,dfinalincleffs[k][i-1]);
			thresholds->getCutVar(i-1,&thiscut);
	
			inclmaxeffbars->GetXaxis()->SetBinLabel(i,thiscut);
			inclfinaleffbars->GetXaxis()->SetBinLabel(i,thiscut);

		}

		inclmaxeffbars->GetXaxis()->LabelsOption("v");
		inclmaxeffbars->SetMinimum(0.0);
		inclmaxeffbars->Draw("HIST,E1");

		inclfinaleffbars->GetXaxis()->LabelsOption("v");
		inclfinaleffbars->SetMinimum(0.0);
		inclfinaleffbars->Draw("HIST,E1,SAME");

		c->Print(outname+"_"+dataset->getName()+".pdf");
		c->Print(outname+"_"+dataset->getName()+".C");


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


	c->SetBottomMargin(0.2);

	TH1D* maxeffbars = new TH1D("maxEffbars","Total Efficiencies per channel",datastore->getNSignalDatasets(),0,datastore->getNSignalDatasets());
	//maxeffbars->SetXTitle("Mode");
	maxeffbars->SetYTitle("Efficiency");
	maxeffbars->SetStats(kFALSE);
	maxeffbars->SetFillStyle(3345);
	maxeffbars->SetLineWidth(2);
	maxeffbars->SetLineColor(kRed);
	maxeffbars->SetFillColor(kRed);


	TH1D* finaleffbars = new TH1D("finalEffbars","Total Efficiencies per channel",datastore->getNSignalDatasets(),0,datastore->getNSignalDatasets());
	//finaleffbars->SetXTitle("Mode");
	finaleffbars->SetYTitle("Efficiency");
	finaleffbars->SetStats(kFALSE);
	finaleffbars->SetFillStyle(3354);
	finaleffbars->SetLineColor(kBlue);
	finaleffbars->SetLineWidth(2);
	finaleffbars->SetLineColor(kBlue);
	finaleffbars->SetFillColor(kBlue);


	thresholds->buildorcut(pars, &bestcut);
	cout << "Mode & $\\epsilon^{\\text{max}}$ & $\\epsilon^{\\text{final}}$ \\\\" << endl;
	datastore->getBackgroundEfficiency(new TString(""),&bestcut,&Beff,&d_Beff);
	for(UInt_t k = 0; k<datastore->getNSignalDatasets(); k++){
		cropdataset * dataset = datastore->getSignalDataset(k);
		dataset->getEfficiency(new TString(""), &bestcut,&Seff,&d_Seff);
		cout << dataset->getName() << " & $"<< prettyPrint(maxeffs[k]) << "\\pm"<< prettyPrint(dmaxeffs[k])<< "$ & $" << prettyPrint(Seff) << "\\pm" << prettyPrint(d_Seff) <<"$ \\\\"<< endl;
		finaleffs[k]=Seff;
		dfinaleffs[k]=d_Seff;



	}
	cout << "Ratelimit & $"<< prettyPrint(ratelimit) << "$ & $" << prettyPrint(Beff) << "\\pm" << prettyPrint(d_Beff) <<"$ \\\\"<< endl;  


	for(UInt_t k = 1; k<=datastore->getNSignalDatasets(); k++){
		maxeffbars->SetBinContent(k,maxeffs[k-1]);
		maxeffbars->SetBinError(k,dmaxeffs[k-1]);
		maxeffbars->GetXaxis()->SetBinLabel(k,datastore->getSignalDataset(k-1)->getName());
		finaleffbars->SetBinContent(k,finaleffs[k-1]);
		finaleffbars->SetBinError(k,dfinaleffs[k-1]);
		finaleffbars->GetXaxis()->SetBinLabel(k,datastore->getSignalDataset(k-1)->getName());
	}
	maxeffbars->GetXaxis()->LabelsOption("v");
	maxeffbars->GetYaxis()->LabelsOption("v");
	maxeffbars->SetMinimum(0.0);
	maxeffbars->Draw("HIST,E1");
	finaleffbars->GetXaxis()->LabelsOption("v");
	finaleffbars->GetYaxis()->LabelsOption("v");
	finaleffbars->SetMinimum(0.0);
	finaleffbars->Draw("HIST,E1,SAME");

	c->Print(outname+"_totals.pdf");
	c->Print(outname+"_totals.C");








	timer.Stop();
	timer.Print();
	exit(0);
}
