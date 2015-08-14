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
if(argc <5 || argc >9){	
	mergevarsInfo();
	exit(EXIT_FAILURE);
	}
	TString varListName = argv[1];
	TString outname = argv[2];
	cout << "INFO: Parsing varfile..." << endl;
	varensemble = new cropvarensemble("input", argv[1]);
	if(argc == 5){varensemble->mergeVars(argv[3],argv[4]);}
	if(argc == 6){varensemble->mergeVars(argv[3],argv[4],argv[5]);}
	if(argc == 7){varensemble->mergeVars(argv[3],argv[4],argv[5],argv[6]);}
	if(argc == 8){varensemble->mergeVars(argv[3],argv[4],argv[5],argv[6],argv[7]);}
	if(argc == 9){varensemble->mergeVars(argv[3],argv[4],argv[5],argv[6],argv[7],argv[8]);}


	varensemble->writeToFile(outname);

}
