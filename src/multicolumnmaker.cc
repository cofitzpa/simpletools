/* multicolumnmaker: Part of the simpletools package
 * (c) Conor Fitzpatrick, 2008, Adam Morris, 2018
 *
 * If you find this program useful in whole or in part
 * please cite this paper:
 *
 * Feel free to send bugreports, feature requests, patches etc to:
 * conor.fitzpatrick@cern.ch
 *
 */
#include "cropcutensemble.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TFormula.h>
#include <TTreeFormula.h>
#include <TTreeFormulaManager.h>
using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
struct column {
	column(TTree* inTree, TTree* outTree, TString columnname, TString formulaname) {
		//Ensure branch already named this way isn't copied:
		inTree->SetBranchStatus(columnname, 0);
		// Create the formula object (ROOT won't let me copy or assign)
		formula = new TTreeFormula(columnname, formulaname, outTree);
		// Create new branch
		formbranch = outTree->Branch(columnname, &branchvalue);
	}
	column(const column& other) {
		formbranch = other.formbranch;
		formbranch->SetAddress(&branchvalue);
		formula = new TTreeFormula(other.formula->GetName(), other.formula->GetTitle(), other.formula->GetTree());
	}
	~column() {
		delete formula;
	}
	void fill() {
		branchvalue = formula->EvalInstance();
		formbranch->Fill();
	}
	Double_t branchvalue;
	TTreeFormula *formula;
	TBranch *formbranch;
};
std::vector<column> make_columns(TTree* inTree, TTree* outTree, TString filename) {
	std::vector<column> columns;
	ifstream fileStream;
	fileStream.open(filename);
	if (!fileStream){
		cout << "FATAL: Error opening file" << filename << endl;
		exit(EXIT_FAILURE);
	}
	TString line;
	while(!line.ReadLine(fileStream).eof()){
		cout << "Parsing " << line << "\n";
		if(!line.BeginsWith("#")){
			TObjArray *lineitems = line.Tokenize(cropcutensembledelimiters);
			columns.push_back(column(inTree, outTree, ((TObjString*)lineitems->At(0))->GetString(), ((TObjString*)lineitems->At(1))->GetString()));
		}
	}
	return columns;
}
int main(int argc, char *argv[]) {
	if(argc != 5 ){
		cout << "multicolumnmaker:   	creates a new column in an ntuple based on a formula" << endl;
		cout << "author:     	   Conor Fitzpatrick, 2008"<< endl;
		cout << "Syntax: " << argv[0] << " <input.root> <path/to/ntuple> <formula file> <output.root>"<< endl;
		return EXIT_FAILURE;
	}

	TString inname = argv[1];
	TString tpath = argv[2];
	TString fname = argv[3];
	TString soutname = argv[4];
	TStopwatch sw;

	cout << "--------MULTICOLUMNMAKER - Conor Fitzpatrick, 2008 ----------" << endl;
	cout << "adding colums using formulae in:		" << fname 	<< endl;
	cout <<	"to ntuple:		" << tpath 	<< endl;
	cout <<	"in file:		" << inname 	<< endl;
	cout << "output file:		" << soutname 	<< endl;

	cout << "-------------------------------------------------------" << endl;

	TFile* in = TFile::Open( inname );

	TTree* inTree = (TTree*)in->Get(tpath);
	UInt_t total = (UInt_t)inTree->GetEntries();

	TString slash = "/";
	TString name = tpath;
	tpath.Resize(std::max(tpath.First(slash),0));

	TFile* sout = TFile::Open(soutname,"RECREATE");
	if(tpath!=name){
	sout->mkdir(tpath);
	sout->cd(tpath);
	}else{
	sout->cd();
	}

	cout << "copying ntuple" << endl; sw.Start();
	TTree *soutTree = inTree->CloneTree(-1);
	cout << "creating new columns" << endl; sw.Start();
	std::vector<column> columns = make_columns(inTree, soutTree, fname);
	int k=0;
	int pc =0;
	for(UInt_t l=0; l<total; l++){
		soutTree->GetEntry(l);
		for(column& col: columns){
			col.fill();
		}
		pc = ((100*l)/total);
		if(pc == k+10){
			k = pc;
			cout << pc << "\% complete\r" << flush;
		}
	}
	cout << "100" << "\% complete\r" << endl;
	soutTree->Write();
	sout->Write();
	sout->Close();
	in->Close();
	sw.Stop();
	cout << "-------------------------------------------------------" << endl;
	sw.Print();
	cout << "done." << endl;
}
