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
#include <cmath>
#include "TROOT.h"
#include <TFile.h>
#include <TSystem.h>
#include<TPaveStats.h>
#include <TStyle.h>
#include<TLegend.h>

vector<TPaveStats*> getStats(THStack* hs){
	TCanvas *d = new TCanvas("2null","2null",0,0);
	vector<TPaveStats*> stats;
	TList *hists = hs->GetHists();
	TListIter next(hists);
	TObject *obj;
	while((obj=next())){
		TH1D* hist = (TH1D*)obj;
		hist->SetStats(true);
		hist->Draw();
		d->Update();
		TPaveStats* stat = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
		stat->SetTextColor(hist->GetFillColor());
		stats.push_back(stat);

	}
	delete d;
	return stats;

}

using std::cout;
using std::endl;
cropdatastore * datastore;
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
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(0.5);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLineStyleString(2,"[12 12]");


	TCanvas *c = new TCanvas("null","null",0,0);
	if(argc != 9){	
	stackerInfo();
	exit(EXIT_FAILURE);
	}

	TString weightListName = argv[1];
	TString varListName = argv[2];
	UInt_t ny = atoi(argv[3]);
	UInt_t nx = atoi(argv[4]);
	UInt_t resy = atoi(argv[5]);
	UInt_t resx = atoi(argv[6]);
	TString outdir = argv[7];
	TString stackargs = argv[8];

	UInt_t ppp = ny * nx;

	cout << "INFO: Parsing weightfile..." << endl;
	datastore = new cropdatastore("test",weightListName);
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
				TPad *ctmp = (TPad*)c->cd(v+1);
				gPad->SetFillColor(0);
				THStack *hs = datastore->getStack(varensemble,nVarsDone);
				vector <TPaveStats *> stats = getStats(hs);
				ctmp->cd();
				hs->Draw(stackargs);

				if(v==0 && stats.size()>5){
				TLegend *leg = (TLegend*)ctmp->BuildLegend(0.70,0.55,0.89,0.89);

				leg->SetLineColor(0);
				leg->SetFillStyle(0);
				}
				
				TString yTitle = "Candidates/";
				yTitle += varensemble->getVarResolution(nVarsDone);
				varensemble->getUnits(nVarsDone,&units);
				if(units !="\"\""){
				hs->GetXaxis()->SetTitle(units);
				yTitle += units;
				}

				hs->GetYaxis()->SetTitle(yTitle);
				if(varensemble->isLogScale(nVarsDone)){

					if(hs->GetMinimum()<0){
					hs->SetMinimum(hs->GetMaximum()*10E-5); //Set log scale to cover 5 orders of magnitude.
					}

					gPad->SetLogy();
				}
			if(stats.size()<6){
				Double_t x1 = 0.825;
				Double_t x2 = 0.975;
				Double_t y1 = 0.825;
				Double_t y2 = 0.975;
				for(UInt_t i=0; i<stats.size(); i++){
					stats[i]->SetX1NDC(x1); stats[i]->SetX2NDC(x2);
					stats[i]->SetY1NDC(y1); stats[i]->SetY2NDC(y2);
					y1 = y1 - 0.165;
					y2 = y2 - 0.165;
					stats[i]->Draw();
				}
			}
				ctmp->Update();
				cout << varensemble->toLine(nVarsDone) << endl;
				nVarsDone++;
			}
		}
		page++;
		outFile->cd();
		c->Write();
		c->Print(outdir+"/"+canvName+".C");
		c->Print(outdir+"/"+canvName+".eps");
		c->Print(outdir+"/"+canvName+".pdf");
		c->Print(outdir+"/"+canvName+".png");
		delete c;

	}
	outFile->Write();
	outFile->Close();

}


