/* columnmaker: Part of the simpletools package
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
int main(int argc, char *argv[]) {
	if(argc != 6 ){
		cout << "columnmaker:   	creates a new column in an ntuple based on a formula" << endl;
		cout << "author:     	   Conor Fitzpatrick, 2008"<< endl;
		cout << "Syntax: " << argv[0] << " <input.root> <path/to/ntuple> <formula> <column name> <output.root>"<< endl;
		return EXIT_FAILURE;
	}

	TFile *in(0);
	TString inname = argv[1];   
	TString tpath = argv[2];   
	TString fname = argv[3];
	TString cname = argv[4];
	TFile *sout(0);
	TString soutname = argv[5];
	TStopwatch sw;

	cout << "--------COLUMNMAKER - Conor Fitzpatrick, 2008 ----------" << endl;
	cout << "adding column:		" << cname 	<< endl;
	cout << "using formula:		" << fname 	<< endl;
	cout <<	"to ntuple:		" << tpath 	<< endl;
	cout <<	"in file:		" << inname 	<< endl;
	cout << "output file:		" << soutname 	<< endl;

	cout << "-------------------------------------------------------" << endl;

	in = TFile::Open( inname );

	TTree* inTree = (TTree*)in->Get(tpath);
	UInt_t total = (UInt_t)inTree->GetEntries();

	//Ensure branch already named this way isn't copied: 
	inTree->SetBranchStatus(cname,0);

	TString slash = "/";
	TString name = tpath;
	tpath.Resize(std::max(tpath.First(slash),0)); 

	sout = new TFile(soutname,"RECREATE");
	if(tpath!=name){
	sout->mkdir(tpath);
	sout->cd(tpath);
	}else{
	sout->cd();
	}


	

	cout << "copying ntuple" << endl; sw.Start();	 
	TTree *soutTree = inTree->CloneTree(-1);
	Float_t val =0.0;
	TBranch *formbranch = soutTree->Branch(cname, &val, cname);
	TTreeFormula *formula = new TTreeFormula("formula",fname,soutTree);
	cout << "creating new column" << endl; sw.Start();	 
	int k=0;
	int pc =0;
	for(UInt_t l=0; l<total; l++){
		soutTree->GetEntry(l);
		val = formula->EvalInstance(0);
		formbranch->Fill();
		pc = ((100*l)/total);
		if(pc == k+10){
			k = pc;
			cout << pc << "\% complete\r" << flush;
		}
	}
	cout << "100" << "\% complete\r" << endl;
	soutTree->Write();	
	sout->Write();
	sw.Stop();
	cout << "-------------------------------------------------------" << endl;
	sw.Print();
	cout << "done." << endl;
}
