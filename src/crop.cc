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
#include "cropoptimisationengine.h"
#include <iostream>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TRandom3.h>
//GUI
#include <TApplication.h>
#include <TGFrame.h>
#include<TRootEmbeddedCanvas.h>


//http://www.akiti.ca/mulmatvec2.html
using std::cout;
using std::endl;
using namespace boost;
int main(int argc, char *argv[]) {
	cropBanner();

	Bool_t gui=true;
	if(argc != 6){	
		if(argc != 7){

				cropInfo();
				exit(EXIT_FAILURE);
		}else{

			TString b = "-b";
			if(argv[6] == b){ 
				gui = false;
			}else{
				cropInfo();
				exit(EXIT_FAILURE);
			}

		}
	}
	
	gROOT->SetStyle("Plain");
	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetStatColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetFillColor(0);
	gStyle->SetFrameFillColor(0);


	std::string weightListName = argv[1];
	std::string cutListName = argv[2];
	UInt_t maxIters = atoi(argv[3]);
	TString mode = argv[4];
	TString outname = argv[5];
	TApplication theApp("App", &argc, argv);
	TFile * outFile = new TFile(outname,"RECREATE");
	TCanvas *guiACanvas;
	if(gui){
		TGMainFrame * guiMain = new TGMainFrame(gClient->GetRoot(),1024,768);
		TRootEmbeddedCanvas * guiCanvas = new TRootEmbeddedCanvas("guicanvas",guiMain,1024,768);
		guiMain->AddFrame(guiCanvas,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,10,10,10,1));
		guiMain->SetWindowName("CROP");
		guiMain->MapSubwindows();
		guiMain->Resize(guiMain->GetDefaultSize());
		guiMain->MapWindow();
		guiACanvas = guiCanvas->GetCanvas();
		guiACanvas->SetDoubleBuffer(1);
	}else{
		guiACanvas = new TCanvas("CROP","CROP", 1024,768);
	}
	guiACanvas->cd();
	guiACanvas->Divide(2,2);

	cout << "INFO: Parsing weightfile..." << endl;
	cropdatastore *datastore = new cropdatastore("test",weightListName);
	datastore->initStats();

	cout << "INFO: Parsing cutfile..." << endl;
	TRandom3 rndGen(0);
	cropcutensemble *cutensemble = new cropcutensemble("test", cutListName, &rndGen);

	TSpecificFunctor<cropcutensemble> SBOrderMethod(cutensemble, &cropcutensemble::OrderAscendingSeffBeff);
	TSpecificFunctor<cropcutensemble> SOrderMethod(cutensemble, &cropcutensemble::OrderAscendingSeff);
	TSpecificFunctor<cropcutensemble> BOrderMethod(cutensemble, &cropcutensemble::OrderDescendingBrej);
	TSpecificFunctor<cropcutensemble> ROrderMethod(cutensemble, &cropcutensemble::OrderRandom);
	TSpecificFunctor<cropcutensemble> OOrderMethod(cutensemble, &cropcutensemble::OrderOriginal);
	TFunctor *OrderMethod = &SBOrderMethod;
	if(mode == "B"){OrderMethod = &BOrderMethod;

	cout << "Using Background rejection ordering" << endl;
	}else{
		if(mode == "R"){OrderMethod = &ROrderMethod;
		
	cout << "Using Random ordering" << endl;
		}else{
			if(mode == "S"){OrderMethod = &SOrderMethod;
	cout << "Using Signal efficiency ordering" << endl;
			}else{
				if(mode == "O"){OrderMethod = &OOrderMethod;
				cout << "Using default ordering" << endl;
				}else{
					if(mode == "SB"){
					cout << "Using Purity Ordering" << endl;
					}else{
						cout << "Ordering method not found. This should be one of S,B,O,R,SB" << endl;
						exit(1);
					}
				}
			}
		}
	}

	TStopwatch timer;
	timer.Start();


	outFile->cd();
	TDirectory *initDir = outFile->mkdir("Initial","Initial");
	initDir->cd();

	//Uncomment to use S/sqrt(S+B) FoM
	cropoptimisationengine *optimisationengine = new cropoptimisationengine(datastore, cutensemble, guiACanvas, &rndGen, *cropFoM_SoRSB, OrderMethod);
	//Uncomment to use Punzi FoM
	//cropoptimisationengine *optimisationengine = new cropoptimisationengine(datastore, cutensemble, guiACanvas, &rndGen, *cropFoM_Punzi, OrderMethod);

	cout << "INFO: Initialising: Finding initial working point..." << endl;
	optimisationengine->initialise();
	outFile->cd();
	outFile->Write();
	TDirectory *finalDir = outFile->mkdir("Optimised","Optimised");
	finalDir->cd();

	cout << "INFO: Starting Optimisation...." << endl;
	bool converged = optimisationengine->optimise(maxIters);
	if(converged){
		cout << "INFO: Optimisation Complete. Printing final cut values and efficiencies:" << endl;
	}else{
		cout << "WARNING: FAILED TO CONVERGE IN THE ALLOCATED NUMBER OF STEPS! RESULT MAY NOT BE OPTIMAL!" << endl;
		cout << "WARNING: Check for correlated cuts (oscillations in the evolution plot)," << endl;
		cout << "WARNING: check for invalid or unreasonable cut ranges" << endl;
		cout << "WARNING: or increase the max. number of steps" << endl;
	}
	outFile->cd();
	optimisationengine->writeEvolutionPlot();

	optimisationengine->finalStats();
	timer.Stop();
	timer.Print();
	outFile->Write();
	outFile->Close();
	exit(0);
}
