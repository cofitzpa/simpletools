/* cuttester: Part of the simpletools package
 * (c) Conor Fitzpatrick, 2008
 *
 * If you find this program useful in whole or in part
 * please cite this paper: LHCb-INT-2009-029
 *
 * Feel free to send bugreports, feature requests, patches etc to:
 * conor.fitzpatrick@cern.ch
 *
 */


#include <stdlib.h>
#include <libgen.h>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <vector>
#include <cropdataset.h>

using std::cout;
using std::string;
using std::endl;
using std::vector;
using std::ifstream;
using std::ofstream;

char pretty[10];
UInt_t ncuts=0;
vector<TString> cuts;
const TString And = "&&";
const TString Or = "||";
TString allCuts = "";
TString weightVar = "";
TTree* inTree;

TString createCut(UInt_t exclude){ //creates a cut string of all optimal cuts except the one to optimise
	TString outString ="";
	for(UInt_t i =0; i<ncuts; i++){
		if(i!=exclude){
			outString.Append("(");
			outString.Append(cuts[i]);
			outString.Append(")");
			outString.Append(Or);
		}
	}
	outString = outString.Strip(TString::kTrailing,*Or);
	return outString;
}


int main(int argc, char *argv[]) {
	if(argc != 5 && argc !=6){
		cout << "thresholdeffs:   	tests a set of L0 thresholds, giving exclusive and total efficiencies"<< endl;
		cout << "author:        Conor Fitzpatrick, 2015"<< endl;
		cout << "Syntax: " << argv[0] << " <input.root> <path/to/ntuple> <thresholds.txt> <Deadtime> [weight var]"<< endl;
		cout << "thresholds.txt is a list of thresholds without spaces, eg: L0HadronMeV>2800"<< endl;
		cout << "Deadtime is the efficiency lost to deadtime, eg: 0 means no deadtime, 1 means 100\% deadtime"<< endl;
		cout << "It is assumed the deadtime is the same for all thresholds" << endl;
		return EXIT_FAILURE;
	}
	Bool_t usingWeights = false;
	if(argc == 6){
		usingWeights = true;
		weightVar = argv[5];
	}

	Double_t ineff = 1.0-atof(argv[4]);
	TString filepath = argv[1];
	TString ntuplepath = argv[2];
	std::string cutListName = argv[3];
	cout << "--------THRESHOLDEFFS - Conor Fitzpatrick, 2015 ----------" << endl;
	cout << "testing cuts in:	" << cutListName << endl;
	cout <<	"to ntuple:		" << ntuplepath 	<< endl;
	cout <<	"in file:		" << filepath 	<< endl;
	cout << "-------------------------------------------------------" << endl;

	string tuplename = basename(argv[1]);
	tuplename.erase(tuplename.find_last_of("."), string::npos);
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
			allCuts.Append(Or);
			cout <<  ncuts << "	" << cut << endl;
			ncuts++;
		}
	}
	allCuts = allCuts.Strip(TString::kTrailing,*Or);
	cout << "-------------------------------------------------------------------------------------------------------------------------------------" << endl;

	TStopwatch sw;
	sw.Start();
	cout << "cut no.		excl. eff." << endl;



	cout << "-------------------------------------------------------------------------------------------------------------------------------------" << endl;
	Double_t excleff, dexcleff, incleff, dincleff, toteff, dtoteff, cands, dcands;
	TString cutssofar = "";
	TString noCut = "1";
	for(UInt_t i=0; i< ncuts; i++){
		TString thiscut = cuts[i];
		dataset->getEfficiency(&noCut,&thiscut,&excleff,&dexcleff);
		excleff*=ineff;
		dexcleff*=ineff;
		cout << i << "	&	$" << prettyPrint(excleff*100.0) << "\\pm" << prettyPrint(dexcleff*100.0) << "\\%$  \\\\" << endl;
	}
	cout << "-------------------------------------------------------------------------------------------------------------------------------------" << endl;
		dataset->getEfficiency(&noCut,&allCuts,&toteff,&dtoteff);
		toteff*=ineff;
		dtoteff*=ineff;

	cout << "dataset                total eff" << endl;
	cout << "-------------------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << tuplename << "  &       $" << prettyPrint(toteff*100.0) << "\\pm" << prettyPrint(dtoteff*100.0) << "\\%$  \\\\" << endl;

	sw.Stop();
	sw.Print();
	cout<< "done." << endl;

}
