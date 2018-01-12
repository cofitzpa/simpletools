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
#include <TSpline.h>
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


Double_t ratelimit, hertz;
TSpline3 *deadtimespline;
Double_t *nodtmaxeffs, *dnodtmaxeffs,*maxeffs, *dmaxeffs, *finaleffs, *dfinaleffs;



std::vector< std::vector<Double_t> > maxincleffs;
std::vector< std::vector<Double_t> > dmaxincleffs;
std::vector< std::vector<Double_t> > finalincleffs;
std::vector< std::vector<Double_t> > dfinalincleffs;


double gradient(double rate){
return rate*0.0001;
}

//Determines 1.0 - the derandomiser deadtime
//deadtimespline is built from scanning the ODIN with a given filling scheme
//note: 1100 is hardcoded, assuming the target rate is for 1.1MHz.
double technicaldeadtime(double rate){
	if(((rate/ratelimit)*hertz)>2000.0){
	return 0.0;
	}
	double dts = deadtimespline->Eval((rate/ratelimit)*hertz);
	if(dts<1){
	return dts;
	}
	return 1.0;
}

//returns 1-deadtime for _physics_ deadtime
//ie: if your thresholds accept more than 1.1MHz
//your signal effs will be throttled by the returned amount
double physicsdeadtime(double rate){
	if(rate < ratelimit){
		return 1.0;
	}
	return ratelimit/rate;
//	return 0.0;
}


//returns the larger of the tech and phys deadtimes
double ratelimiter(double rate){
	return technicaldeadtime(rate)*physicsdeadtime(rate);
}

class initfcn : public TMVA::IFitterTarget {
	public:
		//default constructor
		initfcn(cropdatastore* _dstore, cropdataset* _dset, cropcutensemble * _thresholds){
			dstore = _dstore;
			dset = _dset;
			thresholds = _thresholds;
		}
		//Chi2 to minimise for a single dataset
		double EstimatorFunction( std::vector<double> &par){
			TString cut;
			thresholds->buildorcut(par,&cut);
			dset->getEfficiency(new TString(""),&cut, &Seff,&d_Seff);
			dstore->getBackgroundEfficiency(new TString(""),&cut, &Beff,&d_Beff);
			Seff = Seff * ratelimiter(Beff);
			return ((1.0-Seff)*(1.0-Seff)) + gradient(Beff);
		}

	private:
		cropdatastore* dstore;
		cropdataset* dset;
		cropcutensemble* thresholds;
		Double_t Seff, Beff;
		Double_t d_Seff, d_Beff;
};


class optfcn : public TMVA::IFitterTarget {
	public:
		//default constructor
		optfcn(cropdatastore* _dstore, cropcutensemble * _thresholds){
			dstore = _dstore;
			thresholds = _thresholds;
		}
		//Chi2 to minimise for the ensemble of datasets
		double EstimatorFunction(std::vector<double> &par){
			TString cut;
			thresholds->buildorcut(par,&cut);
			double f=0.0;
			double spec = 0.0;
			dstore->getBackgroundEfficiency(new TString(""),&cut, &Beff,&d_Beff);
			double rlim = ratelimiter(Beff);
			for(UInt_t k = 0; k<dstore->getNSignalDatasets(); k++){
				spec = dstore->getSignalDataset(k)->getSpecial();
				if(spec>0.0){
					dstore->getSignalDataset(k)->getEfficiency(new TString(""),&cut, &Seff,&d_Seff);
					Seff = Seff*rlim;
					f+=  spec*((1.0-(Seff/maxeffs[k]))*(1.0-(Seff/maxeffs[k])));
				}
			}
			f+=gradient(Beff);
			return f;
		}

	private:
		cropdatastore* dstore;
		cropcutensemble* thresholds;
		Double_t Seff, Beff;
		Double_t d_Seff, d_Beff;
};


int main(int argc, char *argv[]) {
	bwdivBanner();
	if(argc != 6){
		if(argc != 7){
		bwdivInfo();
		exit(EXIT_FAILURE);
		}
		hertz = atof(argv[6]);
	}else{
		hertz = 1100.;
	}

	std::string weightListName = argv[1];
	std::string thresholdsName = argv[2];
	ratelimit = atof(argv[3]);
	TFile *deadtimefile = TFile::Open(argv[4]);
	if(!deadtimefile){cout << "Couldn't find deadtime file" << endl;
		exit(0);
	}
	TGraph *deadtimegraph;
	deadtimefile->GetObject("deadtimegraph",deadtimegraph);
	if(!deadtimegraph){
		cout << "Couldn't find graph in file" << endl;
		exit(0);
	}
	deadtimespline = new TSpline3("deadtimespline",deadtimegraph);
	TString outname = argv[5];


	TCanvas* c=new TCanvas("Effs","Effs",1024,768);
	c->SetBottomMargin(0.5);

	cout << "optimising to target rate: " << hertz << " using deadtime " << argv[4] << " target retention " << ratelimit << endl;
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


	nodtmaxeffs = new Double_t [datastore->getNSignalDatasets()];
	dnodtmaxeffs = new Double_t [datastore->getNSignalDatasets()];
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
	TString opt="PopSize=100:Steps=30:Cycles=3:ConvCrit=0.01:SaveBestCycle=5";
	TMVA::IFitterTarget *ffit;
	TMVA::FitterBase* fitter;

	TString bestcut, thiscut;
	Double_t S, d_S, Seff, d_Seff;
	Double_t B, d_B, Beff, d_Beff;


	for(UInt_t k = 0; k<datastore->getNSignalDatasets(); k++){
		cropdataset* dataset = datastore->getSignalDataset(k);
		finalincleffs[k].resize(thresholds->NCutVars);
		dfinalincleffs[k].resize(thresholds->NCutVars);
		maxincleffs[k].resize(thresholds->NCutVars);
		dmaxincleffs[k].resize(thresholds->NCutVars);
		if (dataset->getSpecial()<=0.0){
			maxeffs[k]=0.0;
			dmaxeffs[k]=0.0;
			nodtmaxeffs[k]=0.0;
			dnodtmaxeffs[k]=0.0;
		}else{
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
			cout << dataset->getName() << " MAX EFF: $" << prettyPrint(Seff)<<"\\pm"<<prettyPrint(d_Seff)<< "$ BACKGROUND EFF: $" << prettyPrint(Beff)<<"\\pm"<<prettyPrint(d_Beff) << endl;
			cout << "PHYS DEADTIME: "<< prettyPrint(1.0 - physicsdeadtime(Beff)) << " TECH DEADTIME: " <<  prettyPrint(1.0 - technicaldeadtime(Beff)) << endl;

			maxeffs[k]=Seff*ratelimiter(Beff);
			dmaxeffs[k]=d_Seff*ratelimiter(Beff);
			nodtmaxeffs[k]=Seff;
			dnodtmaxeffs[k]=d_Seff;

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
		cout << dataset->getName() << " & $"<< prettyPrint(nodtmaxeffs[k]) << "\\pm"<< prettyPrint(dnodtmaxeffs[k])<< "$ & $" << prettyPrint(Seff) << "\\pm" << prettyPrint(d_Seff) <<"$ \\\\"<< endl;
		finaleffs[k]=Seff;
		dfinaleffs[k]=d_Seff;



	}

	cout << "Retention & $"<< prettyPrint(ratelimit) << "$ & $" << prettyPrint(Beff) << "\\pm" << prettyPrint(d_Beff) <<"$ \\\\"<< endl;
	cout << "Phys Deadtime & $"<< " 0.0 " << "$ & $" << prettyPrint(1.0 - physicsdeadtime(Beff)) <<"$ \\\\"<< endl;
	cout << "Tech Deadtime & $"<< " 0.0 " << "$ & $" << prettyPrint(1.0 - technicaldeadtime(Beff)) <<"$ \\\\"<< endl;

	for(UInt_t k = 1; k<=datastore->getNSignalDatasets(); k++){
		maxeffbars->SetBinContent(k,nodtmaxeffs[k-1]);
		maxeffbars->SetBinError(k,dnodtmaxeffs[k-1]);
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
