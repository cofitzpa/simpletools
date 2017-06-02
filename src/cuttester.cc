/* cuttester: Part of the simpletools package
 * (c) Conor Fitzpatrick, 2008
 *
 * If you find this program useful in whole or in part
 * please cite this paper:
 *
 * Feel free to send bugreports, feature requests, patches etc to:
 * conor.fitzpatrick@cern.ch
 *
 */


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <vector>
#include <cropdataset.h>
using std::cout;
using std::endl;
using std::vector;
using std::ifstream;
using std::ofstream;

char pretty[10];
UInt_t ncuts=0;
vector<TString> cuts;
const TString And = "&&";
TString allCuts = "";
TString weightVar = "";
TTree* inTree;

TString createCut(UInt_t exclude){ //creates a cut string of all optimal cuts except the one to optimise
	TString outString ="";
	for(UInt_t i =0; i<ncuts; i++){
		if(i!=exclude){
			outString.Append(cuts[i]);
			outString.Append(And);
		}
	}
	outString = outString.Strip(TString::kTrailing,*And);
	return outString;
}


int main(int argc, char *argv[]) {
	if(argc != 4 && argc !=5){
		cout << "cuttester:   	tests a set of cuts, giving inclusive and exclusive efficiencies"<< endl;
		cout << "author:        Conor Fitzpatrick, 2008"<< endl;
		cout << "Syntax: " << argv[0] << " <input.root> <path/to/ntuple> <cutlist.txt> [weight var]"<< endl;
		cout << "cutlist.txt is a list of cuts without spaces, eg: B_s_Mass<500"<< endl;
		return EXIT_FAILURE;
	}
	Bool_t usingWeights = false;
	if(argc == 5){
		usingWeights = true;
		weightVar = argv[4];
	}

	TString filepath = argv[1];
	TString ntuplepath = argv[2];
	std::string cutListName = argv[3];
	cout << "--------CUTTESTER - Conor Fitzpatrick, 2008 ----------" << endl;
	cout << "testing cuts in:	" << cutListName << endl;
	cout <<	"to ntuple:		" << ntuplepath 	<< endl;
	cout <<	"in file:		" << filepath 	<< endl;
	cout << "-------------------------------------------------------" << endl;

	cropdataset *dataset = new cropdataset("data",filepath,ntuplepath,"1",weightVar,1,1,true);
	dataset->print();
	ifstream fileStream;
	TString cut;
	std::string lineread;
	fileStream.open(cutListName.c_str());
	if (!fileStream){
		std::cout << "Error in opening cut file" << std::endl;
		return EXIT_FAILURE;
	}

	cout << "cut no. 	cut" << endl;
	cout << "-------------------------------------------------------" << endl;

	while(!cut.ReadLine(fileStream).eof()){
		if(cut.BeginsWith("#")){}else{	cuts.push_back(cut);

			allCuts.Append(cut);
			allCuts.Append(And);
			cout <<  ncuts << "	" << cut << endl;
			ncuts++;
		}
	}
	allCuts = allCuts.Strip(TString::kTrailing,*And);
	cout << "-------------------------------------------------------------------------------------------------------------------------------------" << endl;

	TStopwatch sw;
	sw.Start();
	cout << "cut no.		excl. eff.			incl. eff.			total eff.			cands. passing" << endl;


	cout << "-------------------------------------------------------------------------------------------------------------------------------------" << endl;
	Double_t excleff, dexcleff, incleff, dincleff, toteff, dtoteff, cands, dcands;
	TString cutssofar = "";
	TString noCut = "1";
	for(UInt_t i=0; i< ncuts; i++){
		TString thiscut = cuts[i];
		cutssofar += thiscut;
		TString allothercuts = createCut(i);
		dataset->getWeightedEntries(&cutssofar,&cands,&dcands);
		dataset->getEfficiency(&noCut,&thiscut,&excleff,&dexcleff);
		dataset->getEfficiency(&allothercuts,&thiscut,&incleff,&dincleff);
		dataset->getEfficiency(&noCut,&cutssofar,&toteff,&dtoteff);
		cutssofar += "&&";
		cout << i << "	&	$" << prettyPrint(excleff*100.0) << "\\pm" << prettyPrint(dexcleff*100.0) << "\\%$	&	$" << prettyPrint(incleff*100.0)<<"\\pm"<<prettyPrint(dincleff*100.0)<<"\\%$	&	$" <<prettyPrint(toteff*100.0) << "\\pm" << prettyPrint(dtoteff*100.0) << "\\%$	&       $" << prettyPrint(cands) <<"\\pm" << prettyPrint(dcands) << "$	\\\\" << endl;
	}
	cout << "-------------------------------------------------------------------------------------------------------------------------------------" << endl;
	sw.Stop();
	sw.Print();
	cout<< "done." << endl;

}
