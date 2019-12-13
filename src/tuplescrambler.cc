/* tuplesampler: Part of the simpletools package
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
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TMVA/Tools.h>
#include <vector>
#include <deque>
#include <list>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>

	template <class T>
inline void INSERT_ELEMENTS (T& coll, int first, int last)
{
	for (int i=first; i<=last; ++i) {
		coll.insert(coll.end(),i);
	}
}
using namespace std;

int main(int argc, char *argv[]) {
	if(argc != 5 ){
		cout << "tuplescrambler:  randomly scrambles an ntuple into a new file"<< endl;
		cout << "author:        Conor Fitzpatrick, 2011"<< endl;
		cout << "Syntax: " << argv[0] << " <input.root> <path/to/ntuple> <seed> <output1.root>"<< endl;
		return EXIT_FAILURE;
	}

	TFile *in;
	TString inname = argv[1];
	TString tpath = argv[2];
	TString name = tpath;
	TMVA::RandomGenerator<TRandom3> rd(atoi(argv[3]));
	TFile *sout1(0);
	TString soutname1 = argv[4];

	TStopwatch sw;
	cout << "opening " << inname << endl;

	in = TFile::Open( inname );


	cout << "getting tree " << tpath << endl;
	TTree* inTree = (TTree*)in->Get(tpath);
	UInt_t total =  inTree->GetEntries();

	TString slash = "/";
	tpath.Resize(std::max(tpath.First(slash),0));
	sout1 = TFile::Open(soutname1,"RECREATE");
	if(name!=tpath){
	sout1->mkdir(tpath);
	sout1->cd(tpath);
	}else{
	sout1->cd();
	}
	TTree *soutTree1 = inTree->CloneTree(0);


	cout << "--------TUPLESCRAMBLER - Conor Fitzpatrick, 2008 ----------" << endl;
	cout << "sampling from:         " << tpath     << endl;
	cout << "in file:               " << inname      << endl;
	cout << "with seed:             " << argv[3]     << endl;
	cout << "to file:               " << soutname1   << endl;
	cout << "---------------------------------------------------------" << endl;
	//inTree->Print();


		sw.Start();
	vector<int> entries;
	INSERT_ELEMENTS(entries,0,total-1);
	std::shuffle(entries.begin(), entries.end(), rd);
	std::vector<int>::iterator it = entries.begin();
	UInt_t i = 0, pc =0, k=0;
	while( it != entries.end() ) {
		inTree->GetEntry(*it);
//		cout << i << " " << *it  << endl;
		soutTree1->Fill();
		++it;
		i++;
		pc = ((100*i)/total);
		if(pc == k+10){
			k = pc;
		cout << pc << "\% complete\r" << flush;
		}
	}


	cout << "100" << "\% complete\r" << endl;
	UInt_t out1 =  soutTree1->GetEntries();

	sout1->Write();
	sw.Stop();

	cout << "--------------------------------------------------------" << endl;
	cout << "Input contained   " << total << " events" << endl;
	cout << "Output 1 contains " << out1 << " events" << endl;
	cout << "--------------------------------------------------------" << endl;
	sw.Print();
//	sort (entries.begin(), entries.end());
//	it = entries.begin();
//	i = 0;
//	while( it != entries.end() ) {
//	cout << i-*it << " " << i << " " << *it  << endl;
//	i++;
//	++it;
//	}
	cout << "done." << endl;

}
