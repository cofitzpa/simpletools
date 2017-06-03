#include "cropvarensemble.h"

using std::ifstream;
using std::ofstream;
using std::max;
using std::min;
using std::string;
cropvarensemble::cropvarensemble(cropcutensemble *cuts, bool _log){
	name = "default";
	Nvars = 0;
	TString _var;
	TString _units = "\"\"";
	for(UInt_t n = 0; n<cuts->NCutVars; n++){
		cuts->getCutVar(n,&_var);
		if(cuts->isAGreaterThanCut(n)){
			_var.Remove(TString::kTrailing, *cropcutensembleGreaterThan);
		}else{
			_var.Remove(TString::kTrailing, *cropcutensembleLessThan);
		}
		addVar(_var,_units,cuts->getCutMinVal(n),cuts->getCutMaxVal(n),cuts->getCutSteps(n),_log,false);
	}

}

cropvarensemble::cropvarensemble(){
	name = "default";
	Nvars = 0;
}

cropvarensemble::cropvarensemble(TString _name){
	name = _name;
	Nvars = 0;
}
void cropvarensemble::addVar(TString _var, TString _units, Double_t _min, Double_t _max, Double_t _res, bool _log, bool _useless){
	vars.push_back(_var);
	units.push_back(_units);
	mins.push_back(_min);
	maxs.push_back(_max);
	ress.push_back(_res);
	bins.push_back((UInt_t)((_max-_min)/_res));
	varOrder.push_back(Nvars);
	logScale.push_back(_log);
	useless.push_back(_useless);
	Nvars++;

}
void cropvarensemble::addVar(TString _var, TString _units, Double_t _min, Double_t _max, UInt_t _bins, bool _log, bool _useless){
	vars.push_back(_var);
	units.push_back(_units);
	mins.push_back(_min);
	maxs.push_back(_max);
	bins.push_back(_bins);
	ress.push_back((_max - _min)/((Double_t)_bins));
	varOrder.push_back(Nvars);
	logScale.push_back(_log);
	useless.push_back(_useless);
	Nvars++;
}
void cropvarensemble::print() const{
	cout << "var	number	units	min	max	res	bins"	<< endl;
	cout << "-----------------------------------------------------------"	<< endl;
	for(UInt_t i =0; i<Nvars; i++){
		printOneLine(i);
	}

}
void cropvarensemble::printOneLine(UInt_t num) const{
	cout << toLine(num) << endl;

}

cropvarensemble::cropvarensemble(TString _name, TString _filename){
	name = _name;
	Nvars = 0;

	UInt_t linenum = 0;
	ifstream fileStream;
	fileStream.open(_filename);
	if (!fileStream){
		cout << "FATAL: Error in opening varfile" << _filename << endl;
		exit(EXIT_FAILURE);
	}
	TString line;
	while(!line.ReadLine(fileStream).eof()){
		linenum++;
		if(!line.BeginsWith("#")){
			this->addVar(line, _filename,linenum);

		}
	}
	fileStream.close();
}

void cropvarensemble::addVar(TString line, TString _filename, UInt_t linenum){
	TString _var, _units, logstring;
	Double_t _min, _max;
	UInt_t _bins;
	bool _log;

	TObjArray* Strings = line.Tokenize(cropvarensembledelimiters);
	if(Strings->GetEntries() != 6){
		cout << "FATAL: " << _filename << " line " << linenum << ": expected 6 entries, found " << Strings->GetEntries() << ". Check the line is correct and try again." << endl;
		exit(EXIT_FAILURE);
	}
	TIter iString(Strings);
	TObjString* os=(TObjString*)iString(); //var
	try{
		_var = lexical_cast <string> (os->GetString());
	}catch(bad_lexical_cast &e){
		cout << "FATAL: "<< _filename << " line " << linenum << ": not a valid var" << endl;
		exit(EXIT_FAILURE);
	}
	os=(TObjString*)iString(); //min
	try{
		_min = lexical_cast <double> (os->GetString());
	}catch(bad_lexical_cast &e){
		cout << "FATAL: " << _filename << ": line " << linenum << ": not a valid lower bound" << endl;
		exit(EXIT_FAILURE);
	}
	os=(TObjString*)iString(); //max
	try{
		_max = lexical_cast <double> (os->GetString());
	}catch(bad_lexical_cast &e){

		cout << "FATAL: " << _filename << ": line " << linenum << ": not a valid upper bound" << endl;
		exit(EXIT_FAILURE);
	}
	os=(TObjString*)iString(); //bins
	try{
		_bins = lexical_cast <unsigned int> (os->GetString());
	}catch(bad_lexical_cast &e){
		cout << "FATAL: " << _filename << ": line " << linenum << ": not a valid binning" << endl;
		exit(EXIT_FAILURE);
	}
	os=(TObjString*)iString(); //units
	try{
		_units = lexical_cast <string> (os->GetString());
	}catch(bad_lexical_cast &e){
		cout << "FATAL: " << _filename << ": line " << linenum << ": not a valid units" << endl;
		exit(EXIT_FAILURE);
	}

	os=(TObjString*)iString(); //log/lin
	try{
		logstring = lexical_cast <string> (os->GetString());
	}catch(bad_lexical_cast &e){
		cout << "FATAL: " << _filename << ": line " << linenum << ": log or lin not specified" << endl;
		exit(EXIT_FAILURE);
	}
	if(logstring == "log"){
		_log= true;
	}else{
		if(logstring == "lin"){
			_log = false;
		}else{
			cout << "FATAL: " << _filename << ": line " << linenum << ": log or lin not specified" << endl;
			exit(EXIT_FAILURE);
		}
	}
	addVar(_var,_units,_min,_max,_bins,_log, false);
}
cropvarensemble::cropvarensemble(cropdataset *_data, bool _logy, bool _logx, UInt_t _bins){
	name = "default";
	Nvars = 0;
	Double_t _min, _max;
	TString _var;
	TString _units = "\"\"";
	bool _useless;
	TTree * ntuple = _data->getTree();
	vector <TString> vars = VarsFromDataSet(_data);
	for(int i = 0; i < vars.size(); i++){
		_var = vars[i];
		if(_logx){
			_var += ")";
			_var.Prepend("log10(");
		}
		ntuple->Draw(_var+">>tmp",_data->getWeightVar());
		TH1F *tmp = (TH1F*)gDirectory->Get("tmp");
		_max = tmp->GetXaxis()->GetXmax();
		_min = tmp->GetXaxis()->GetXmin();
		if((_max == _min) || (tmp->GetRMS()==0.0)){
			_useless = true;
		}else{
			_useless = false;
		}
		tmp->Delete();
		addVar(_var, _units, _min, _max, _bins, _logy, _useless);
	}
}

vector <TString> cropvarensemble::VarsFromList(TString _filename){
	vector <TString> vars;
	ifstream fileStream;
	fileStream.open(_filename);
	if (!fileStream){
		cout << "FATAL: Error in opening varfile" << _filename << endl;
		exit(EXIT_FAILURE);
	}

	TString line;
	TString thisvar;

	while(!line.ReadLine(fileStream).eof()){
		if(!line.BeginsWith("#")){
			TObjArray* Strings = line.Tokenize(cropvarensembledelimiters);

			TIter iString(Strings);
			TObjString* os=(TObjString*)iString(); //var
			try{
				thisvar = lexical_cast <string> (os->GetString());
			}catch(bad_lexical_cast &e){
				cout << "FATAL: "<< _filename << " line " << line << ": not a valid var" << endl;
				exit(EXIT_FAILURE);
			}

			vars.push_back(thisvar);
		}
	}
	fileStream.close();
	return vars;
}

vector <TString> cropvarensemble::VarsFromDataSet(cropdataset *_data){
	TTree * ntuple = _data->getTree();
	vector <TString> vars;
	TObjArray *members = ntuple->GetListOfLeaves();
	for(int i = 0; i < members->GetEntries(); i++){
		vars.push_back(members->At(i)->GetName());
	}
	return vars;
}

cropvarensemble::cropvarensemble(cropdatastore *_datastore, TString _varfile, bool _logx, bool _logy, UInt_t _bins){
	name = "default";
	Nvars = 0;
	UInt_t nSignal = _datastore->getNSignalDatasets();
	UInt_t nBackground = _datastore->getNBackgroundDatasets();
	vector <TString> allvars;
	allvars = cropvarensemble::VarsFromList(_varfile);
	vector<TString>::const_iterator p, p_end;
	std::sort(allvars.begin(), allvars.end(), tstringlt);
	p_end = std::unique(allvars.begin(), allvars.end(),tstringcmp);
	TString _units = "\"\"";
	TString _var;
	for(p = allvars.begin(); p < p_end; p++){
		_var = *p;
		if(_logx){
			_var += ")";
			_var.Prepend("log10(");
		}
		Double_t _min=0, _max=0;
		Double_t tmpmax=0, tmpmin=0;
		bool drawn = false;
		bool _useless = false;
		Int_t success = 0;
		for(UInt_t s = 0; s<nSignal; s++){
			success = _datastore->getSignalDataset(s)->getTree()->Draw(_var+">>tmp",_datastore->getSignalDataset(s)->getWeightVar());
			if(success == -1){
				cout << "WARNING: finding range for " << _var  << " in dataset " << _datastore->getSignalDataset(s)->getName() << " failed" << endl;
				_useless = true;
			}else{
				TH1F *tmp = (TH1F*)gDirectory->Get("tmp");
				tmpmax = tmp->GetXaxis()->GetXmax();
				tmpmin = tmp->GetXaxis()->GetXmin();
				tmp->Delete();
				if(drawn){
					if(tmpmax > _max){_max = tmpmax;}
					if(tmpmin < _min){_min = tmpmin;}
				}else{
					_max = tmpmax;
					_min = tmpmin;
					drawn = true;
				}
			}
		}

		for(UInt_t b = 0; b<nBackground; b++){
			success = _datastore->getBackgroundDataset(b)->getTree()->Draw(_var+">>tmp",_datastore->getBackgroundDataset(b)->getWeightVar());
			if(success == -1){
				cout << "WARNING: finding range for " << _var  << " in dataset " << _datastore->getBackgroundDataset(b)->getName() << " failed" << endl;
				_useless = true;
			}else{
				TH1F *tmp = (TH1F*)gDirectory->Get("tmp");
				tmpmax = tmp->GetXaxis()->GetXmax();
				tmpmin = tmp->GetXaxis()->GetXmin();
				tmp->Delete();
				if(drawn){
					if(tmpmax > _max){_max = tmpmax;}
					if(tmpmin < _min){_min = tmpmin;}
				}else{
					_max = tmpmax;
					_min = tmpmin;
					drawn = true;
				}
			}


		}
		if(!_useless){
			if(_max == _min){
				_useless = true;

				cout << "WARNING: range of " << _var << " is zero for all datasets!" << endl;
			}else{
				TH1D *test = _datastore->getHisto(_var, _bins, _min, _max);
				if(test->GetRMS() == 0.0){
					_useless = true;

					cout << "WARNING: RMS of " << _var << " is zero!" << endl;
				}
			}
		}
		addVar(_var, _units, _min, _max, _bins, _logy, _useless);

	}
}

cropvarensemble::cropvarensemble(cropdatastore *_datastore, bool _logx, bool _logy, UInt_t _bins){
	name = "default";
	Nvars = 0;
	UInt_t nSignal = _datastore->getNSignalDatasets();
	UInt_t nBackground = _datastore->getNBackgroundDatasets();
	vector <TString> allvars, thesevars;
	TObjArray *members;
	if(nSignal>0){
		for(UInt_t s = 0; s<nSignal; s++){
			thesevars = VarsFromDataSet(_datastore->getSignalDataset(s));
			for(int i = 0; i < thesevars.size(); i++){
				allvars.push_back(thesevars[i]);
			}
		}

	}
	if(nBackground>0){
		for(UInt_t b = 0; b<nBackground; b++){
			thesevars = VarsFromDataSet(_datastore->getBackgroundDataset(b));
			for(int i = 0; i < thesevars.size(); i++){
				allvars.push_back(thesevars[i]);
			}

		}
	}
	vector<TString>::const_iterator p, p_end;
	std::sort(allvars.begin(), allvars.end(), tstringlt);
	p_end = std::unique(allvars.begin(), allvars.end(),tstringcmp);       // remove duplicates

	TString _units = "\"\"";
	TString _var;
	for(p = allvars.begin(); p < p_end; p++){
		_var = *p;
		if(_logx){
			_var += ")";
			_var.Prepend("log10(");

		}
		Double_t _min=0, _max=0;
		Double_t tmpmax=0, tmpmin=0;
		bool drawn = false;
		bool _useless = false;
		Int_t success = 0;
		for(UInt_t s = 0; s<nSignal; s++){
			success = _datastore->getSignalDataset(s)->getTree()->Draw(_var+">>tmp",_datastore->getSignalDataset(s)->getWeightVar());
			if(success == -1){
				cout << "WARNING: finding range for " << _var  << " in dataset " << _datastore->getSignalDataset(s)->getName() << " failed" << endl;
				_useless = true;
			}else{
				TH1F *tmp = (TH1F*)gDirectory->Get("tmp");
				tmpmax = tmp->GetXaxis()->GetXmax();
				tmpmin = tmp->GetXaxis()->GetXmin();
				tmp->Delete();
				if(drawn){
					if(tmpmax > _max){_max = tmpmax;}
					if(tmpmin < _min){_min = tmpmin;}
				}else{
					_max = tmpmax;
					_min = tmpmin;
					drawn = true;
				}
			}
		}

		for(UInt_t b = 0; b<nBackground; b++){
			success = _datastore->getBackgroundDataset(b)->getTree()->Draw(_var+">>tmp",_datastore->getBackgroundDataset(b)->getWeightVar());
			if(success == -1){
				cout << "WARNING: finding range for " << _var  << " in dataset " << _datastore->getBackgroundDataset(b)->getName() << " failed" << endl;
				_useless = true;
			}else{
				TH1F *tmp = (TH1F*)gDirectory->Get("tmp");
				tmpmax = tmp->GetXaxis()->GetXmax();
				tmpmin = tmp->GetXaxis()->GetXmin();
				tmp->Delete();
				if(drawn){
					if(tmpmax > _max){_max = tmpmax;}
					if(tmpmin < _min){_min = tmpmin;}
				}else{
					_max = tmpmax;
					_min = tmpmin;
					drawn = true;
				}
			}


		}
		if(!_useless){
			if(_max == _min){
				_useless = true;

				cout << "WARNING: range of " << _var << " is zero for all datasets!" << endl;
			}else{
				TH1D *test = _datastore->getHisto(_var, _bins, _min, _max);
				if(test->GetRMS() == 0.0){
					_useless = true;

					cout << "WARNING: RMS of " << _var << " is zero!" << endl;
				}
			}
		}
		addVar(_var, _units, _min, _max, _bins, _logy, _useless);

	}
}

void cropvarensemble::writeToFile(TString _filename) const{
	ofstream fileStream;
	fileStream.open(_filename);
	fileStream << "#var\tmin\tmax\tbins\tunits\tlog/lin\n";
	for(UInt_t v = 0; v < Nvars; v++){
		fileStream << toLine(varOrder[v]) << "\n";
	}
	fileStream.close();
}
TString cropvarensemble::toLine(UInt_t n) const{
	TString out;
	UInt_t v = n;
	if(useless[v]){
		out = cropvarensembleCstr;
		out += vars[v];
	}else{
		out = vars[v];
	}
	out += cropvarensembleTab;
	out += mins[v];
	out += cropvarensembleTab;
	out += maxs[v];
	out += cropvarensembleTab;
	out += bins[v];
	out += cropvarensembleTab;
	out += units[v];
	out += cropvarensembleTab;
	if(logScale[v]){
		out += "log";
	}else{
		out += "lin";
	}
	return out;
}

UInt_t cropvarensemble::findVar(TString * _var)const{
	vector<TString>::const_iterator is = std::find(vars.begin(), vars.end(), *_var);
	return is - vars.begin();
}


TH1D * cropvarensemble::getSignalHisto(UInt_t n, cropdatastore *_data){
	UInt_t v = varOrder[n];
	return _data->getSignalHisto(vars[v],bins[v],mins[v],maxs[v]);
}
TH1D * cropvarensemble::getBackgroundHisto(UInt_t n, cropdatastore *_data){
	UInt_t v = varOrder[n];
	return _data->getBackgroundHisto(vars[v],bins[v],mins[v],maxs[v]);

}

THStack * cropvarensemble::getStack(UInt_t n, cropdatastore *_data){
	UInt_t v = varOrder[n];
	return _data->getStack(vars[v],bins[v],mins[v],maxs[v]);
}

TH1D * cropvarensemble::getHisto(UInt_t n, cropdatastore *_data){
	UInt_t v = varOrder[n];
	return _data->getHisto(vars[v],bins[v],mins[v],maxs[v]);
}


void cropvarensemble::sortbySepPower(cropdatastore *_data){
	vector<point1D> seps;
	for(UInt_t v = 0; v < Nvars; v++){
		TH1D *S = _data->getSignalHisto(vars[v],bins[v],mins[v],maxs[v]);
		TH1D *B = _data->getBackgroundHisto(vars[v],bins[v],mins[v],maxs[v]);
		point1D point;
		point.x = v;
		point.y = GetSeparation(*S,*B);
		if(point.y == 0.0){
			cout << "WARNING: Separation power is zero for " << vars[v] << endl;
		}
		S->Delete();
		B->Delete();
		seps.push_back(point);
	}
	std::sort(seps.begin(), seps.end(), pointgt);
	cout << "var\tseparation" << endl;
	cout << "------------------------------------" << endl;
	Double_t min = seps[Nvars].y;
	Double_t max = seps[0].y;
	for(UInt_t v = 0; v < Nvars; v++){
		varOrder[v] = seps[v].x;
		cout << vars[seps[v].x] << "\t" << prettyPrint(seps[v].y) << endl;
		if(seps[v].y < 1.0e-10){
			useless[seps[v].x] = true;
		}
	}

}

void cropvarensemble::mergeVars(TString v1, TString v2, TString v3, TString v4, TString v5, TString v6){
	UInt_t _Nvars = Nvars;
	TString newVar1;
	TString newVar2;
	TString tmpString2, tmpString3, tmpString4, tmpString5, tmpString6;
	Int_t v2Pos =0;
	Int_t v3Pos =0;
	Int_t v4Pos =0;
	Int_t v5Pos =0;
	Int_t v6Pos =0;

	for(UInt_t n = 0; n< _Nvars; n++){
		if(vars[n].Contains(v1)){
			tmpString2 = vars[n];
			tmpString2.ReplaceAll(v1,v2);
			v2Pos = findVar(&tmpString2);
			if(v2Pos <0) continue;
			tmpString3 = vars[n];
			tmpString3.ReplaceAll(v1,v3);
			v3Pos = findVar(&tmpString3);
			if(v3Pos <0) continue;
			tmpString4 = vars[n];
			tmpString4.ReplaceAll(v1,v4);
			v4Pos = findVar(&tmpString4);
			if(v4Pos <0) continue;
			tmpString5 = vars[n];
			tmpString5.ReplaceAll(v1,v5);
			v5Pos = findVar(&tmpString5);
			if(v5Pos <0) continue;
			tmpString6 = vars[n];
			tmpString6.ReplaceAll(v1,v6);
			v6Pos = findVar(&tmpString6);
			if(v6Pos <0) continue;
			newVar1 = "max(max(max(max(max((" + vars[n] + "),(" + tmpString2 + ")),("+ tmpString3 + ")),("+ tmpString4 +")),("+ tmpString5 +")),("+ tmpString6 +"))";
			newVar2 = "min(min(min(min(min((" + vars[n] + "),(" + tmpString2 + ")),("+ tmpString3 + ")),("+ tmpString4 +")),("+ tmpString5 +")),("+ tmpString6 +"))";
			addVar(newVar1, units[n], max(max(max(max(max(mins[n],mins[v2Pos]),mins[v3Pos]),mins[v4Pos]),mins[v5Pos]),mins[v6Pos]), max(max(max(max(max(maxs[n],maxs[v2Pos]),maxs[v3Pos]),maxs[v4Pos]),maxs[v5Pos]),maxs[v6Pos]), bins[n], logScale[n], useless[n]);
			addVar(newVar2, units[n], min(min(min(min(min(mins[n],mins[v2Pos]),mins[v3Pos]),mins[v4Pos]),mins[v5Pos]),mins[v6Pos]), min(min(min(min(min(maxs[n],maxs[v2Pos]),maxs[v3Pos]),maxs[v4Pos]),maxs[v5Pos]),maxs[v6Pos]), bins[n], logScale[n], useless[n]);
			useless[n] = true;
		}
		if(vars[n].Contains(v2)||vars[n].Contains(v3)||vars[n].Contains(v4)||vars[n].Contains(v5)||vars[n].Contains(v6)){
			useless[n] = true;
		}
	}
}

void cropvarensemble::mergeVars(TString v1, TString v2, TString v3, TString v4, TString v5){
	UInt_t _Nvars = Nvars;
	TString newVar1;
	TString newVar2;
	TString tmpString2, tmpString3, tmpString4, tmpString5;
	Int_t v2Pos =0;
	Int_t v3Pos =0;
	Int_t v4Pos =0;
	Int_t v5Pos =0;
	for(UInt_t n = 0; n< _Nvars; n++){
		if(vars[n].Contains(v1)){
			tmpString2 = vars[n];
			tmpString2.ReplaceAll(v1,v2);
			v2Pos = findVar(&tmpString2);
			if(v2Pos <0) continue;
			tmpString3 = vars[n];
			tmpString3.ReplaceAll(v1,v3);
			v3Pos = findVar(&tmpString3);
			if(v3Pos <0) continue;
			tmpString4 = vars[n];
			tmpString4.ReplaceAll(v1,v4);
			v4Pos = findVar(&tmpString4);
			if(v4Pos <0) continue;
			tmpString5 = vars[n];
			tmpString5.ReplaceAll(v1,v5);
			v5Pos = findVar(&tmpString5);
			if(v5Pos <0) continue;
			newVar1 = "max(max(max(max((" + vars[n] + "),(" + tmpString2 + ")),("+ tmpString3 + ")),("+ tmpString4 +")),("+ tmpString5 +"))";
			newVar2 = "min(min(min(min((" + vars[n] + "),(" + tmpString2 + ")),("+ tmpString3 + ")),("+ tmpString4 +")),("+ tmpString5 +"))";
			addVar(newVar1, units[n], max(max(max(max(mins[n],mins[v2Pos]),mins[v3Pos]),mins[v4Pos]),mins[v5Pos]), max(max(max(max(maxs[n],maxs[v2Pos]),maxs[v3Pos]),maxs[v4Pos]),maxs[v5Pos]), bins[n], logScale[n], useless[n]);
			addVar(newVar2, units[n], min(min(min(min(mins[n],mins[v2Pos]),mins[v3Pos]),mins[v4Pos]),mins[v5Pos]), min(min(min(min(maxs[n],maxs[v2Pos]),maxs[v3Pos]),maxs[v4Pos]),maxs[v5Pos]), bins[n], logScale[n], useless[n]);
			useless[n] = true;
		}
		if(vars[n].Contains(v2)||vars[n].Contains(v3)||vars[n].Contains(v4)||vars[n].Contains(v5)){
			useless[n] = true;
		}
	}
}

void cropvarensemble::mergeVars(TString v1, TString v2, TString v3, TString v4){
	UInt_t _Nvars = Nvars;
	TString newVar1;
	TString newVar2;
	TString tmpString2, tmpString3, tmpString4;
	UInt_t v2Pos =0;
	UInt_t v3Pos =0;
	UInt_t v4Pos =0;
	for(UInt_t n = 0; n< _Nvars; n++){
		if(vars[n].Contains(v1)){
			tmpString2 = vars[n];
			tmpString2.ReplaceAll(v1,v2);
			v2Pos = findVar(&tmpString2);
			if(v2Pos == 0) continue;
			tmpString3 = vars[n];
			tmpString3.ReplaceAll(v1,v3);
			v3Pos = findVar(&tmpString3);
			if(v3Pos == 0) continue;
			tmpString4 = vars[n];
			tmpString4.ReplaceAll(v1,v4);
			v4Pos = findVar(&tmpString4);
			if(v4Pos == 0) continue;
			newVar1 = "max(max(max((" + vars[n] + "),(" + tmpString2 + ")),("+ tmpString3 + ")),("+ tmpString4 +"))";
			newVar2 = "min(min(min((" + vars[n] + "),(" + tmpString2 + ")),("+ tmpString3 + ")),("+ tmpString4 +"))";
			addVar(newVar1, units[n], max(max(max(mins[n],mins[v2Pos]),mins[v3Pos]),mins[v4Pos]), max(max(max(maxs[n],maxs[v2Pos]),maxs[v3Pos]),maxs[v4Pos]), bins[n], logScale[n], useless[n]);
			addVar(newVar2, units[n], min(min(min(mins[n],mins[v2Pos]),mins[v3Pos]),mins[v4Pos]), min(min(min(maxs[n],maxs[v2Pos]),maxs[v3Pos]),maxs[v4Pos]), bins[n], logScale[n], useless[n]);
			useless[n] = true;
		}
		if(vars[n].Contains(v2)||vars[n].Contains(v3)||vars[n].Contains(v4)){
			useless[n] = true;
		}
	}
}


void cropvarensemble::mergeVars(TString v1, TString v2, TString v3){

	UInt_t _Nvars = Nvars;
	TString newVar1;
	TString newVar2;
	TString tmpString2, tmpString3;
	UInt_t v2Pos =0;
	UInt_t v3Pos =0;
	for(UInt_t n = 0; n< _Nvars; n++){
		if(vars[n].Contains(v1)){
			tmpString2 = vars[n];
			tmpString2.ReplaceAll(v1,v2);
			v2Pos = findVar(&tmpString2);
			if(v2Pos == 0) continue;
			tmpString3 = vars[n];
			tmpString3.ReplaceAll(v1,v3);
			v3Pos = findVar(&tmpString3);
			if(v3Pos == 0) continue;
			newVar1 = "max(max((" + vars[n] + "),(" + tmpString2 + ")),("+ tmpString3 + "))";
			newVar2 = "min(min((" + vars[n] + "),(" + tmpString2 + ")),("+ tmpString3 + "))";
			addVar(newVar1, units[n], max(max(mins[n],mins[v2Pos]),mins[v3Pos]), max(max(maxs[n],maxs[v2Pos]),maxs[v3Pos]), bins[n], logScale[n], useless[n]);
			addVar(newVar2, units[n], min(min(mins[n],mins[v2Pos]),mins[v3Pos]), min(min(maxs[n],maxs[v2Pos]),maxs[v3Pos]), bins[n], logScale[n], useless[n]);
			useless[n] = true;
		}
		if(vars[n].Contains(v2)||vars[n].Contains(v3)){
			useless[n] = true;
		}
	}
}

void cropvarensemble::mergeVars(TString v1, TString v2){

	UInt_t _Nvars = Nvars;
	TString newVar1;
	TString newVar2;
	TString tmpString2;
	UInt_t v2Pos =0;
	for(UInt_t n = 0; n< _Nvars; n++){
		if(vars[n].Contains(v1)){
			tmpString2 = vars[n];
			tmpString2.ReplaceAll(v1,v2);
			v2Pos = findVar(&tmpString2);
			if(v2Pos == 0) continue;
			newVar1 = "max((" + vars[n] + "),(" + tmpString2 + "))";
			newVar2 = "min((" + vars[n] + "),(" + tmpString2 + "))";
			addVar(newVar1, units[n], max(mins[n],mins[v2Pos]), max(maxs[n],maxs[v2Pos]), bins[n], logScale[n], useless[n]);
			addVar(newVar2, units[n], min(mins[n],mins[v2Pos]), min(maxs[n],maxs[v2Pos]), bins[n], logScale[n], useless[n]);
			useless[n] = true;
		}
		if(vars[n].Contains(v2)){
			useless[n] = true;
		}
	}
}
