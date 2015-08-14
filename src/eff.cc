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
#include "cropdatastore.h"
#include "cropvarensemble.h"
#include <stdlib.h>
#include <iostream>
#include "stdio.h"
#include "TROOT.h"
#include <TFile.h>
#include <TSystem.h>
#include<TPaveStats.h>
#include <TStyle.h>

vector<TPaveStats*> getStats(THStack* hs);

using std::cout;
using std::endl;
cropdatastore * datastorenum, * datastoredenom;
cropvarensemble * varensemble;
int main(int argc, char *argv[]) {
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetStatColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetFillColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	//	gStyle->SetMarkerStyle(20);
	//	gStyle->SetMarkerSize(1);
	//	gStyle->SetHistLineWidth(1);
	//	gStyle->SetLineStyleString(2,"[12 12]");


	TCanvas *c = new TCanvas("null","null",0,0);
	if(argc != 10){
		effInfo();
		exit(EXIT_FAILURE);
	}


	TString weightListName1 = argv[1];
	TString weightListName2 = argv[2];

	TString varListName = argv[3];
	UInt_t ny = atoi(argv[4]);
	UInt_t nx = atoi(argv[5]);
	UInt_t resy = atoi(argv[6]);
	UInt_t resx = atoi(argv[7]);
	TString outdir = argv[8];
	TString stackargs = argv[9];

	UInt_t ppp = ny * nx;

	cout << "INFO: Parsing weightfile..." << endl;
	UInt_t color = 1;
	datastorenum = new cropdatastore("numerator",weightListName1);
	if(datastorenum->getNSignalDatasets()==1 && datastorenum->getNBackgroundDatasets()==0){
		datastorenum->setName(datastorenum->getSignalDataset(0)->getName());

		color = datastorenum->getSignalDataset(0)->getColor();
	}else{
		if(datastorenum->getNSignalDatasets()==0 && datastorenum->getNBackgroundDatasets()==1){
			datastorenum->setName(datastorenum->getBackgroundDataset(0)->getName());
			color = datastorenum->getBackgroundDataset(0)->getColor();
		}
	}
	datastoredenom = new cropdatastore("denominator",weightListName2);
	if(datastoredenom->getNSignalDatasets()==1 && datastoredenom->getNBackgroundDatasets()==0){datastoredenom->setName(datastoredenom->getSignalDataset(0)->getName());}else{
		if(datastoredenom->getNSignalDatasets()==0 && datastoredenom->getNBackgroundDatasets()==1){datastoredenom->setName(datastoredenom->getBackgroundDataset(0)->getName());}
	}
	varensemble = new cropvarensemble("test", varListName);

	gSystem->mkdir( outdir );
	TFile *outFile = new TFile(outdir+"/plots.root","RECREATE");

	UInt_t nVarsDone = 0, page = 0;
	TString units;
	while(nVarsDone<varensemble->Nvars){
		cout << "page =" << page << endl;
		outFile->cd();
		char *canvName = new char[256];
		sprintf(canvName,outdir+"_page_%04i",page);
		TCanvas *c = new TCanvas(canvName,canvName,resy,resx);
		c->SetFillColor(0);
		c->Divide(ny,nx);
		c->Update();

		for(UInt_t v = 0; v<ppp; v++){
			if(nVarsDone<varensemble->Nvars){
				c->cd(v+1);
				gPad->SetFillColor(0);
				TH1D* num = datastorenum->getHisto(varensemble,nVarsDone);
				TH1D* denom = datastoredenom->getHisto(varensemble,nVarsDone);
				num->Divide(num, denom, 1, 1, "B");
				delete denom;
				num->SetStats(false);
				num->SetName(varensemble->getVar(nVarsDone));
				num->Draw("E1");
				num->SetTitle(varensemble->getVar(nVarsDone));
				varensemble->getUnits(nVarsDone,&units);
				num->GetXaxis()->SetTitle(units);
				num->GetYaxis()->SetTitle("#varepsilon("+datastorenum->getName()+"/"+datastoredenom->getName()+")");
				num->Draw(stackargs);
				num->SetMinimum(0);
				if(num->GetMaximum()>2.0){
					num->SetMaximum(2.0);
				}else{

					num->SetMaximum(num->GetMaximum()*1.1);
				}
				num->SetLineColor(color);
				num->SetMarkerColor(color);
				num->SetMarkerStyle(20);
				num->SetLineWidth(2);

				if(varensemble->isLogScale(nVarsDone)){

					gPad->SetLogy();
					num->SetMinimum(1.0e-32);

				}
				num->Write();
				c->Update();
				cout << varensemble->toLine(nVarsDone) << endl;
				nVarsDone++;
			}
		}
		page++;
		c->Write();
		c->Print(outdir+"/"+canvName+".png");
		c->Print(outdir+"/"+canvName+".eps");
		c->Print(outdir+"/"+canvName+".pdf");
		c->Print(outdir+"/"+canvName+".svg");
		delete c;

	}
	outFile->Write();
	outFile->Close();

}
