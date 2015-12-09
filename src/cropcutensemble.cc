#include "cropcutensemble.h"
using std::ifstream;
using std::ofstream;
using std::string;

cropcutensemble::cropcutensemble(){
	NCutVars = 0;
}
cropcutensemble::cropcutensemble(TString _name, TRandom3 *_rng){
	name = _name;
	rng = _rng;
	NCutVars = 0;
}

cropcutensemble::cropcutensemble(cropvarensemble *_varensemble, cropdatastore *_data){
	name = _varensemble->getName();
	NCutVars = 0;
	TString cut;
	for(UInt_t n = 0; n<_varensemble->Nvars; n++){
		if(!_varensemble->isUseless(n)){
			cut = _varensemble->getVar(n);

			if(cut.Contains("max(")){cut += "<";}else{
				if(cut.Contains("min(")){cut += ">";}else{
					if(cut.Contains("Max$(")){cut += "<";}else{
						if(cut.Contains("Min$(")){cut += ">";}else{
							TH1D *S = _varensemble->getSignalHisto(n,_data);
							TH1D *B = _varensemble->getBackgroundHisto(n,_data);
							if(S->GetMean()>B->GetMean()){
								cut += ">";
							}else{
								cut += "<";
							}
							S->Delete();
							B->Delete();
						}
					}
				}
			}
			addCut(cut, _varensemble->getVarMinVal(n),_varensemble->getVarMaxVal(n),_varensemble->getVarResolution(n));
		}

	}
}

void cropcutensemble::addCut(TString _cutval, Double_t _cutmin, Double_t _cutmax, Double_t _cutres){
	CutVars.push_back(_cutval);
	CutMinVals.push_back(_cutmin);
	CutMaxVals.push_back(_cutmax);
	CutResolutions.push_back(_cutres);
	cropcutspace space(_cutres,_cutmin,_cutmax);
	CutSpaces.push_back(space);
	CutOrder.push_back(NCutVars);
	NCutVars++;

}
void cropcutensemble::addCut(TString _cutval, Double_t _cutmin, Double_t _cutmax, UInt_t _cutsteps){
	CutVars.push_back(_cutval);
	CutMinVals.push_back(_cutmin);
	CutMaxVals.push_back(_cutmax);
	CutResolutions.push_back((_cutmax - _cutmin)/(Double_t)_cutsteps);
	cropcutspace space(_cutsteps,_cutmin,_cutmax);
	CutSpaces.push_back(space);
	CutOrder.push_back(NCutVars);
	NCutVars++;
}
void cropcutensemble::print() const{
	cout << "order	number	Sig. Eff	Bkg. Rej	FoM	cut"	<< endl;
	cout << "-----------------------------------------------------------"	<< endl;
	for(UInt_t i =0; i<NCutVars; i++){
		printOneLine(i);
	}

}
void cropcutensemble::printOneLine(UInt_t num) const{
	UInt_t n = CutOrder[num];
	TString optcut;
	getOptimalCut(num, &optcut);
	cout << num << "	" << n << "	" << prettyPrint(CutSpaces[n].getOptimalSeff())<<"+/-"<< prettyPrint(CutSpaces[n].getOptimaldSeff()) << "	" << prettyPrint(CutSpaces[n].getOptimalBrej()) << "+/-" << prettyPrint(CutSpaces[n].getOptimaldBrej()) << "	" << prettyPrint(CutSpaces[n].getMaxFoM()) << "+/-" << prettyPrint(CutSpaces[n].getdMaxFoM()) << "	" << optcut  << endl;
}

TString cropcutensemble::toLine(UInt_t n) const{
	TString out;
	UInt_t v = CutOrder[n];
	out = CutVars[v];
	out += cropcutensembleTab;
	out += CutMinVals[v];
	out += cropcutensembleTab;
	out += CutMaxVals[v];
	out += cropcutensembleTab;
	out += CutResolutions[v];
	return out;
}

void cropcutensemble::writeToFile(TString _filename)const{
	ofstream fileStream;
	fileStream.open(_filename);
	fileStream << "#cut\tmin\tmax\tres\n";
	for(UInt_t v = 0; v < NCutVars; v++){
		fileStream << toLine(v) << "\n";
	}
	fileStream.close();


}

bool cropcutensemble::isAGreaterThanCut(UInt_t n) const{
	return CutVars[n].EndsWith(">");
}

cropcutensemble::cropcutensemble(TString _name, TString _filename, TRandom3* _rng){
	rng = _rng;
	name = _name;
	NCutVars = 0;
	UInt_t linenum = 0;
	ifstream fileStream2;
	fileStream2.open(_filename);
	if (!fileStream2){
		cout << "FATAL: Error in opening cutfile" << _filename << endl;
		exit(EXIT_FAILURE);
	}
	TString line;
	while(!line.ReadLine(fileStream2).eof()){
		linenum++;
		if(!line.BeginsWith("#")){
			this->addCut(line, _filename,linenum);

		}
	}
	fileStream2.close();
}

void cropcutensemble::addCut(TString line, TString _filename, UInt_t linenum){
	TString cut;
	Double_t lower, upper, res;

	TObjArray* Strings = line.Tokenize(cropcutensembledelimiters);
	if(Strings->GetEntries() != 4 && Strings->GetEntries() != 5){
		cout << "FATAL: " << _filename << " line " << linenum << ": expected 4 entries, found " << Strings->GetEntries() << ". Check the line is correct and try again." << endl;
		exit(EXIT_FAILURE);
	}
	TIter iString(Strings);
	TObjString* os=(TObjString*)iString();
	try{
		cut = lexical_cast <string> (os->GetString());
	}catch(bad_lexical_cast &e){
		cout << "FATAL: "<< _filename << " line " << linenum << ": not a valid cut" << endl;
		exit(EXIT_FAILURE);
	}
	os=(TObjString*)iString();
	try{
		lower = lexical_cast <double> (os->GetString());
	}catch(bad_lexical_cast &e){
		cout << "FATAL: " << _filename << ": line " << linenum << ": not a valid lower bound" << endl;
		exit(EXIT_FAILURE);
	}
	os=(TObjString*)iString();
	try{
		upper = lexical_cast <double> (os->GetString());
	}catch(bad_lexical_cast &e){

		cout << "FATAL: " << _filename << ": line " << linenum << ": not a valid upper bound" << endl;
		exit(EXIT_FAILURE);
	}
	os=(TObjString*)iString();
	try{
		res = lexical_cast <double> (os->GetString());
	}catch(bad_lexical_cast &e){
		cout << "FATAL: " << _filename << ": line " << linenum << ": not a valid resolution" << endl;
		exit(EXIT_FAILURE);
	}
	addCut(cut, lower, upper, res);
}

void cropcutensemble::getOptimalEnsemble(UInt_t _skipme, TString *cut)const{
	*cut = "";
	TString add = "";
	for(UInt_t i=0; i<NCutVars; i++){
		if(i!=_skipme){
			getOptimalCut(i, &add);
			cut->Append(add);
			cut->Append(cropcutensembleAnd);
		}
	}
	cut->Remove(TString::kTrailing,*cropcutensembleAnd);
}

bool cropcutensemble::OrderAscendingSeff(){
	bool sane = false;
	Double_t firstseff = CutSpaces[0].getOptimalSeff();
	vector <bool> used(NCutVars,false);
	Double_t minseff;
	UInt_t minN = NCutVars+1;
	for(UInt_t i=0; i<NCutVars; i++){
		minseff =9999;
		for (UInt_t j = 0; j<NCutVars; j++){
			if(used[j]==false){
				Double_t Seff = CutSpaces[j].getOptimalSeff();
				if(Seff != firstseff){sane = true;}
				if(Seff<minseff){
					minseff = Seff;
					minN = j;
				}

			}
		}
		CutOrder[i] = minN;
		used[minN] = true;
	}
	return sane;
}
bool cropcutensemble::OrderDescendingBrej(){
	bool sane = false; 
	Double_t firstbrej = CutSpaces[0].getOptimalBrej();
	vector <bool> used(NCutVars,false);
	Double_t maxbrej;
	UInt_t maxN = NCutVars+1;
	for(UInt_t i=0; i<NCutVars; i++){
		maxbrej =-9999;
		for (UInt_t j = 0; j<NCutVars; j++){
			if(used[j]==false){
				Double_t Brej = CutSpaces[j].getOptimalBrej();
				if(Brej != firstbrej){sane = true;}
				if(Brej>maxbrej){
					maxbrej = Brej;
					maxN = j;
				}
			}
		}
		CutOrder[i] = maxN;
		used[maxN] = true;
	}
	return sane;
}

bool cropcutensemble::OrderAscendingSeffBeff(){
	bool sane = false;
	vector <bool> used(NCutVars,false);
	Double_t minseffbeff;
	Double_t firstseffbeff = CutSpaces[0].getOptimalSeff() * (1.0 - CutSpaces[0].getOptimalBrej());
	UInt_t minN = NCutVars+1;
	for(UInt_t i=0; i<NCutVars; i++){
		minseffbeff =9999;
		for (UInt_t j = 0; j<NCutVars; j++){
			if(used[j]==false){
				Double_t SeffBeff = CutSpaces[j].getOptimalSeff() * (1.0 - CutSpaces[j].getOptimalBrej());
				if(firstseffbeff != SeffBeff){sane = true;}
				if(SeffBeff<minseffbeff){
					minseffbeff = SeffBeff;
					minN = j;
				}
			}
		}		
		CutOrder[i] = minN;
		used[minN] = true;
	}
	return sane;
}

bool cropcutensemble::OrderRandom(){
	vector <bool> used(NCutVars,false);
	UInt_t rndN = UInt_t(rng->Rndm()*NCutVars);
	for(UInt_t i=0; i<NCutVars; i++){
		while(used[rndN]==true){
			rndN = UInt_t(rng->Rndm()*NCutVars);
		}
		CutOrder[i] = rndN;
		used[rndN] = true;
	}
	return true;
}
bool cropcutensemble::OrderOriginal(){
	for(UInt_t i=0; i<NCutVars; i++){
		CutOrder[i] = i;
	}
	return true;	
}
void cropcutensemble::writeAllPlots(TCanvas *c){
	for(UInt_t n=0; n<NCutVars; n++){
		plotCut(n,c);
		c->SetName(CutVars[CutOrder[n]]);
		c->Write();
	}
}

void cropcutensemble::buildorcut(Double_t *par, TString *orcut){
	*orcut = "(";
	for(UInt_t i = 0; i<NCutVars; i++){
		*orcut += "(";
		*orcut += CutVars[i];
		*orcut += par[i];
		*orcut += ")||";
	}
	orcut->Remove(TString::kTrailing,*cropcutensembleOr);
	*orcut+=")";
}
void cropcutensemble::buildorcut(std::vector <Double_t> par, TString *orcut){
	*orcut = "(";
	for(UInt_t i = 0; i<NCutVars; i++){
		*orcut += "(";
		*orcut += CutVars[i];
		*orcut += par[i];
		*orcut += ")||";
	}
	orcut->Remove(TString::kTrailing,*cropcutensembleOr);
	*orcut+=")";
}
