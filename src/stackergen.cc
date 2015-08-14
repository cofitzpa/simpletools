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

using std::cout;
using std::endl;
cropdatastore * datastore;
cropvarensemble * varensemble;
int main(int argc, char *argv[]) {
	TCanvas *c = new TCanvas("null","null",0,0);
	if(argc != 6){	
		cout << "usage: stackergen <inputweightfile> <outputvarfile> <bins> <log/lin (x)> <log/lin (y)>" << endl;	
	}

	TString weightListName = argv[1];
	TString outname = argv[2];
	UInt_t bins = atoi(argv[3]);
	TString loglinx = argv[4];
	TString logliny = argv[5];
	bool usinglogsx = false;
	bool usinglogsy = false;
	if(loglinx == "log"){
		usinglogsx = true;
	}else{
		if(loglinx != "lin"){
			cout << "ERROR: argument 3 must be either \"log\" or \"lin\"" << endl;
			return EXIT_FAILURE;
		}
	}
	if(logliny == "log"){
		usinglogsy = true;
	}else{
		if(logliny != "lin"){
			cout << "ERROR: argument 4 must be either \"log\" or \"lin\"" << endl;
			return EXIT_FAILURE;
		}
	}
	cout << "INFO: Parsing weightfile..." << endl;
	datastore = new cropdatastore("test",weightListName);
	datastore->initStats();
	varensemble = new cropvarensemble(datastore, usinglogsx, usinglogsy, bins);
	varensemble->writeToFile(outname);

}
