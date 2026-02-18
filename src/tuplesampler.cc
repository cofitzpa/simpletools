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
using std::cout;
using std::endl;
using std::flush;

int main(int argc, char *argv[]) {
	if(argc != 7 ){
		cout << "tuplesampler:  randomly samples an ntuple into two files of size specified by ratio"<< endl;
		cout << "author:        Conor Fitzpatrick, 2008"<< endl;
		cout << "Syntax: " << argv[0] << " <input.root> <path/to/ntuple> <ratio> <seed> <output1.root> <output2.root>"<< endl;
		cout << "If <ratio> is greater than 1, it is assumed to be the desired yield in output1." << endl;
		return EXIT_FAILURE;
	}

	TFile *in;
	TString inname = argv[1];
	TString tpath = argv[2];
	TString name = tpath;
	Double_t ratio = atof(argv[3]);
	TFile *sout1(0);
	TFile *sout2(0);
	TString soutname1 = argv[5];
	TString soutname2 = argv[6];
	TRandom3*  rndGen = new TRandom3(atoi(argv[4]));
	TNamed* seed = new TNamed("seed", argv[4]);
	TStopwatch sw;
	cout << "opening " << inname << endl;

	in = TFile::Open( inname );
	TString slash = "/";

	cout << "getting tree " << tpath << endl;
	TTree* inTree = (TTree*)in->Get(tpath);
	tpath.Resize(std::max(tpath.First(slash),0));
	UInt_t total =  inTree->GetEntries();
	if(ratio>1.0){
	ratio = ratio/total;
	}
	if(ratio>1.0){
	ratio=1.0;
	}
	sout1 = TFile::Open(soutname1,"RECREATE",0);
	seed->Write();
	if(name!=tpath){
	sout1->mkdir(tpath);
	sout1->cd(tpath);
	}else{
	sout1->cd();
	}
	TTree *soutTree1 = inTree->CloneTree(0);

	sout2 = TFile::Open(soutname2,"RECREATE");
	seed->Write();
	if(name!=tpath){
	sout2->mkdir(tpath);
	sout2->cd(tpath);
	}else{
		sout2->cd();
	}
	TTree *soutTree2 = inTree->CloneTree(0);

	cout << "--------TUPLESAMPLER - Conor Fitzpatrick, 2008 ----------" << endl;
	cout << "sampling from:         " << tpath     << endl;
	cout << "in file:               " << inname      << endl;
	cout << "with ratio:            " << ratio     << endl;
	cout << "to file:               " << soutname1   << endl;
	cout << "remainder in:          " << soutname2   << endl;
	cout << "---------------------------------------------------------" << endl;
	//inTree->Print();
	int k = 0;
	int pc = 0;
	sw.Start();
	for(UInt_t i = 0; i<total; i++){
		inTree->GetEntry(i);
		if(rndGen->Rndm() < ratio){
			soutTree1->Fill();
		}else{
			soutTree2->Fill();
		}
		pc = ((100*i)/total);
		if(pc == k+10){
			k = pc;
			cout << pc << "\% complete\r" << flush;
		}
	}
	cout << "100" << "\% complete\r" << endl;
	UInt_t out1 =  soutTree1->GetEntries();
	UInt_t out2 =  soutTree2->GetEntries();

	sout1->Write();
	sout2->Write();
	sw.Stop();

	cout << "--------------------------------------------------------" << endl;
	cout << "Input contained   " << total << " events" << endl;
	cout << "Output 1 contains " << out1 << " events" << endl;
	cout << "Output 2 contains " << out2 << " events" << endl;
	Double_t outRat = Double_t(out1)/Double_t(total);
	Double_t doutRat = sqrt(outRat*(1.0-outRat)/Double_t(total));
	cout << "Sampled ratio: " << outRat<<"+\\-"<< doutRat << endl;

	cout << "--------------------------------------------------------" << endl;
	sw.Print();
	cout << "done." << endl;

}
