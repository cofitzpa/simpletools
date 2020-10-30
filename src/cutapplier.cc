/* cutapplier: Part of the simpletools package
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
using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
	if(argc != 5 && argc != 6){
		cout << "cutapplier:          Applies a cut to an input ntuple"<< endl;
		cout << "author:              Conor Fitzpatrick, 2008"<< endl;
		cout << "Syntax:              " << argv[0] << " <input.root> [<input_friend.root>] <path/to/ntuple> <cut string> <accepted_output.root>"<< endl;
		cout << "<input_friend.root>: Optional, friend added to <input.root>. Not merged into <accepted_output.root>" << endl;
		return EXIT_FAILURE;
	}

	//PDG 2006 Unified conf. interval (95%)
	Double_t cl95m[11] = {0.00,0.05,0.36,0.82,1.37,1.84,2.21,2.58,2.94,4.36,4.75};
	Double_t cl95p[11] = {3.09,5.14,6.72,8.25,9.76,11.26,12.75,13.81,15.29,16.77,17.82};


	TString inname = "";
	TString ffriend = "";
	TString tpath = "";
	TString name = "";
	TString cname = "";
	TString soutname = "";

	if(argc == 5){
	  inname = argv[1];
	  ffriend = "";
	  tpath = argv[2];
	  name = tpath;
	  cname = argv[3];
	  soutname = argv[4];
	}
	else if(argc == 6){
	  inname = argv[1];
	  ffriend = argv[2];
	  tpath = argv[3];
	  name = tpath;
	  cname = argv[4];
	  soutname = argv[5];
	}
	Bool_t clerrs = false;
	Double_t accepted = 0, rejected = 0, total = 0;
	Double_t acceptedm = 0, rejectedm = 0, totalm = 0;
	Double_t acceptedp = 0, rejectedp = 0, totalp = 0;
	Double_t daccepted = 0, drejected = 0, dtotal = 0;
	Double_t eff = 0, rej = 0;
	Double_t effp = 0, rejp = 0;
	Double_t effm = 0, rejm = 0;
	Double_t deff = 0, drej = 0;
	TStopwatch sw;

	cout << "--------CUTAPPLIER - Conor Fitzpatrick, 2008 ----------" << endl;
	cout << "applying cut:		" << cname 	<< endl;
	cout <<	"to ntuple:		" << tpath 	<< endl;
	cout <<	"in file:		" << inname 	<< endl;
	cout <<	"friend (optional):	" << ffriend 	<< endl;
	cout << "output file:		" << soutname 	<< endl;
	cout << "-------------------------------------------------------" << endl;

	TFile* in = TFile::Open( inname );
	TString slash = "/";

	TTree* inTree;
	in->GetObject(tpath,inTree);
	if(!inTree){
		cout << "Error opening the specified ntuple- is the path you specified correct?" << endl;
		exit(1);
	}
	if(ffriend!=""){
	  inTree->AddFriend(tpath,ffriend);
	}
	tpath.Resize(std::max(tpath.First(slash),0));
	total = (Double_t)inTree->GetEntries();
	if(total>10){
		dtotal = sqrt((Double_t)total);
		totalm = total;
		totalp = total;
	}else{
		totalm = cl95m[(UInt_t)total];
		totalp = cl95p[(UInt_t)total];
		clerrs = true;
		cout << "WARNING: Total events in input ntuple < 10. Using 95\% confidence limits" << endl;
	}
	Int_t aftotal = inTree->GetEntries(cname);
	if(aftotal!=0){
		TFile* sout = TFile::Open(soutname,"RECREATE");
		if(name!=tpath){
			sout->mkdir(tpath);
			sout->cd(tpath);
		}else{
			sout->cd();
		}
		cout << "applying cut..." << endl; sw.Start();
		TTree *soutTree = inTree->CopyTree(cname);
		accepted = (Double_t)soutTree->GetEntries();
		if(argc == 6){
		  soutTree->GetListOfFriends()->Clear("nodelete");
		}
		sout->Write();
		sout->Close();
	}else{
	accepted = 0.0;
	}
	in->Close();

	if(accepted>10){
		daccepted = sqrt(accepted);
		acceptedm = accepted;
		acceptedp = accepted;
	}else{
		acceptedm =  cl95m[(UInt_t)accepted];
		acceptedp =  cl95p[(UInt_t)accepted];
		if(clerrs == false){
			clerrs = true;
			cout << "WARNING: Total accepted events < 10. Using 95\% confidence limits" << endl;
		}
	}

	rejected = total - accepted;

	if(rejected>10){
		drejected = sqrt(rejected);
		rejectedm = rejected;
		rejectedp = rejected;
	}else{
		rejectedm =  cl95m[(UInt_t)rejected];
		rejectedp =  cl95p[(UInt_t)rejected];
		if(clerrs == false){
			clerrs = true;
			cout << "WARNING: Total rejected events < 10. Using 95\% confidence limits" << endl;
		}
	}


	sw.Stop();


	eff = accepted/total;
	rej = rejected/total;
	deff = sqrt(eff*(1-eff)/total);
	drej = sqrt(rej*(1-rej)/total);
	if(clerrs == false){



		cout << "-------------------------------------------------------" << endl;

		cout << "total:    	" << total    << "	+/-	" << dtotal <<	endl;

		cout << "accepted:	" << accepted << "	+/-	" << daccepted << endl;

		cout << "rejected:	" << rejected << "	+/-	" << drejected << endl;

		cout << "efficiency:	" << eff      << "	+/-	" << deff << endl;

		cout << "reject. rate: 	" << rej      << "	+/-	" << drej << endl;
		cout << "-------------------------------------------------------" << endl;
		sw.Print();
		cout << "done" << endl;
		cout << "-------------------------------------------------------" << endl;
	}else{
		effp = acceptedp/totalp;
		effm = acceptedm/totalm;
		rejp = rejectedp/totalp;
		rejm = rejectedm/totalm;

		cout << "-------------------------------------------------------" << endl;

		if(total>10){
			cout << "total:    	" << total    << "	+/-	" << dtotal <<	endl;
		}else{
			cout << "total:	(" << total <<")	"<< totalm << "	< 95\%CL <	"<< totalp <<  endl;
		}
		if(accepted>10){
			cout << "accepted:	" << accepted << "	+/-	" << daccepted << endl;

		}else{
			cout << "accepted:	(" << accepted <<")	"<< acceptedm << "	< 95\%CL <	"<< acceptedp <<  endl;
		}
		if(rejected>10){
			cout << "rejected:	" << rejected << "	+/-	" << drejected << endl;
		}else{
			cout << "rejected:      (" << rejected <<")     "<< rejectedm << "      < 95\%CL <      "<< rejectedp <<  endl;
		}
		cout << "efficiency:	" << eff      << "	+/-	" << deff << endl;
		cout << "reject. rate: 	" << rej      << "	+/-	" << drej << endl;
		cout << "-------------------------------------------------------" << endl;
		sw.Print();
		cout << "done" << endl;
		cout << "-------------------------------------------------------" << endl;

	}
}
