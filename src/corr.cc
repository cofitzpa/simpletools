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


using std::cout;
using std::endl;
cropdatastore * datastore;
cropvarensemble * varensemble1, *varensemble2;
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
    set_plot_style();
    //	gStyle->SetMarkerStyle(20);
    //	gStyle->SetMarkerSize(1);
    //	gStyle->SetHistLineWidth(1);
    //	gStyle->SetLineStyleString(2,"[12 12]");


    TCanvas *c = new TCanvas("null","null",0,0);
    if(argc != 11){
	corrInfo();
	exit(EXIT_FAILURE);
    }

    TString weightListName = argv[1];
    TString varListName1 = argv[2];
    TString varListName2 = argv[3];
    UInt_t ny = atoi(argv[4]);
    UInt_t nx = atoi(argv[5]);
    UInt_t resy = atoi(argv[6]);
    UInt_t resx = atoi(argv[7]);
    TString outdir = argv[8];
    TString stackargs = argv[9];
    Double_t maxcorr = atof(argv[10]);

    UInt_t ppp = ny * nx;

    cout << "INFO: Parsing weightfile..." << endl;

    datastore = new cropdatastore("test",weightListName);

    varensemble1 = new cropvarensemble("test", varListName1);
    varensemble2 = new cropvarensemble("test", varListName2);

    gSystem->mkdir( outdir );

    UInt_t nVarsX = varensemble1->Nvars;
    UInt_t nVarsY = varensemble2->Nvars;

    TFile *outFile = TFile::Open(outdir+"/plots.root","RECREATE");
    TH2D *corrs = new TH2D("correlations","correlations",nVarsX+1,0,nVarsX+1,nVarsY+1,0,nVarsY+1);
    UInt_t nVarsDone = 0, page = 0;
    UInt_t nVarsDonepp = ppp;
    TString units, vX,vY;
    Double_t corr;
    TText *t;
    if(varListName1 == varListName2){
	vX = varensemble1->getVar(0);
	corrs->Fill(vX,vX,1.0);
	for(UInt_t nX = 1; nX<nVarsX; nX++){
	    if(!varensemble2->isUseless(nX)){
		vX = varensemble1->getVar(nX);
		corrs->Fill(vX,vX,1.0);
		for(UInt_t nY = 0; nY<nX; nY++){
		    if(!varensemble2->isUseless(nY)){
			vY = varensemble2->getVar(nY);
			if(nVarsDonepp<ppp){
			    c->cd(nVarsDonepp+1);
			    gPad->SetRightMargin ( 0.1 );
			    gPad->SetLeftMargin ( 0.1 );
			    gPad->SetTopMargin ( 0.1 );
			    gPad->SetBottomMargin ( 0.1 );
			    //gPad->SetLogz();
			    TH2D* Histo = datastore->getHisto(varensemble1,varensemble2,nX,nY);
			    corr = Histo->GetCorrelationFactor();
			    Histo->SetTitle("");
			    Histo->SetStats(false);
			    Histo->Draw(stackargs);
			    Histo->GetXaxis()->SetTitle(vX);
			    Histo->GetYaxis()->SetTitle(vY);
			    Histo->SetMinimum(0.0);
			    corrs->Fill(vX,vY,corr);
			    corrs->Fill(vY,vX,corr);
			    t = new TText( 0.1, 0.95, "Correlation: "+prettyPrint(corr));
			    t->SetNDC();
			    t->SetTextSize( 0.05 );
			    t->AppendPad();
			    c->Update();
			    cout << "pad = " << nVarsDonepp << "      Variable = " << vX<<" : "<<vY << "  Correlation = " << prettyPrint(corr);
			    if(fabs(corr)>maxcorr){
				varensemble2->makeUseless(nX);
				cout << "\t Removing First Variable";
			    }
			    cout << endl;
			    nVarsDonepp++;
			    nVarsDone++;
			}else{
			    if(nVarsDone>0){
				c->Write();
				c->Print(outdir+"/"+c->GetTitle()+".png");
				c->Print(outdir+"/"+c->GetTitle()+".eps");
				c->Print(outdir+"/"+c->GetTitle()+".pdf");
				c->Print(outdir+"/"+c->GetTitle()+".C");
				delete c;
			    }
			    char *canvName = new char[256];
			    sprintf(canvName,outdir+"_page_%04i",page);
			    c = new TCanvas(canvName,canvName, resy,resx);
			    c->SetFillColor(0);
			    c->Divide(ny,nx);
			    c->Update();
			    nVarsDonepp = 0;
			    page++;
			    nY--;
			}
		    }
		}
	    }
	}
    }else{
	for(UInt_t nX = 0; nX<nVarsX; nX++){
	    vX = varensemble1->getVar(nX);
	    for(UInt_t nY = 0; nY<nVarsY; nY++){
		vY = varensemble2->getVar(nY);
		if(vX != vY){
		    if(nVarsDonepp<ppp){
			c->cd(nVarsDonepp+1);
			gPad->SetRightMargin ( 0.1 );
			gPad->SetLeftMargin ( 0.1 );
			gPad->SetTopMargin ( 0.1 );
			gPad->SetBottomMargin ( 0.1 );
			//gPad->SetLogz();
			TH2D* Histo = datastore->getHisto(varensemble1,varensemble2,nX,nY);
			corr = Histo->GetCorrelationFactor();
			Histo->SetTitle("");
			Histo->SetStats(false);
			Histo->Draw(stackargs);
			Histo->GetXaxis()->SetTitle(vX);
			Histo->GetYaxis()->SetTitle(vY);
			corrs->Fill(vX,vY,corr);
			t = new TText( 0.1, 0.95, "Correlation: "+prettyPrint(corr));
			t->SetNDC();
			t->SetTextSize( 0.05 );
			t->AppendPad();
			c->Update();
			cout << "pad = " << nVarsDonepp << "      Variable = " << vX<<" : "<<vY << "  Correlation = " << prettyPrint(corr);
			if(fabs(corr)>maxcorr){
			    varensemble2->makeUseless(nY);
			    cout << "\t Removing Second Variable";
			}
			cout << endl;
			nVarsDonepp++;
			nVarsDone++;
		    }else{
			if(nVarsDone>0){
			    c->Write();
			    c->Print(outdir+"/"+c->GetTitle()+".png");
			    c->Print(outdir+"/"+c->GetTitle()+".eps");
			    c->Print(outdir+"/"+c->GetTitle()+".pdf");
			    c->Print(outdir+"/"+c->GetTitle()+".C");
			    delete c;
			}
			char *canvName = new char[256];
			sprintf(canvName,outdir+"_page_%04i",page);
			c = new TCanvas(canvName,canvName, resy,resx);
			c->SetFillColor(0);
			c->Divide(ny,nx);
			c->Update();
			nVarsDonepp = 0;
			page++;
			nY--;


		    }
		}else{
		    corrs->Fill(vX,vY,1.0);
		}
	    }
	}
    }
    c->Write();
    c->Print(outdir+"/"+c->GetTitle()+".png");
    c->Print(outdir+"/"+c->GetTitle()+".eps");
    c->Print(outdir+"/"+c->GetTitle()+".pdf");
    c->Print(outdir+"/"+c->GetTitle()+".C");
    delete c;
    corrs->LabelsDeflate("X");
    corrs->LabelsDeflate("Y");
    corrs->LabelsOption("u");

    TCanvas *d = new TCanvas("Correlations","Correlations", resy,resx);
    d->SetLeftMargin  ( 0.3 );
    d->SetBottomMargin( 0.3 );
    d->SetRightMargin ( 0.1 );
    d->SetTopMargin   ( 0.1 );
    corrs->SetStats(false);
    corrs->GetZaxis()->SetLabelSize( 0.03 );
    corrs->GetXaxis()->SetLabelSize( 0.02 );
    corrs->GetYaxis()->SetLabelSize( 0.03 );
    corrs->SetMaximum(1.0);
    corrs->SetMinimum(-1.0);
    corrs->Draw("COLZ");
    d->Update();
    d->Write();
    d->Print(outdir+"/Correlations.png");
    d->Print(outdir+"/Correlations.eps");
    d->Print(outdir+"/Correlations.pdf");
    d->Print(outdir+"/Correlations.C");
    if(varListName1 == varListName2){
	varensemble2->writeToFile(outdir+"/uniquevars.txt");
    }
    outFile->Write();
    outFile->Close();

}
