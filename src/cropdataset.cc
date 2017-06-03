#include "cropdataset.h"
using std::string;
cropdataset::cropdataset(){}

cropdataset::cropdataset(TString _name, TString _filePath, TString _ntuplePath, TString _preprocCutVar, TString _perEventWeightVar, UInt_t _color=0, UInt_t _fill =0, bool _signal=true){
	name = _name;
	filePath = _filePath;
	ntuplePath = _ntuplePath;
	signal = _signal;
	preprocCutVar = _preprocCutVar;
	fill = _fill;
	color = _color;
	if(_perEventWeightVar==""){
		perEventWeightVar="1";
	}else{
		perEventWeightVar = _perEventWeightVar;
	}
	originalTFile = TFile::Open(filePath,"READ");
	originalNtuple = (TTree*)originalTFile->Get(ntuplePath);
	totalEntries = (Double_t)originalNtuple->GetEntries();
	if(_perEventWeightVar=="A"){
		perEventWeightVar="1/";
		perEventWeightVar.Append(toString(totalEntries));
	}
	procTFile = originalTFile;
	procNtuple = originalNtuple;
	getWeightedEntries(new TString(""),&totalWeightedEntries, &d_totalWeightedEntries);
	procNtuple->Draw(">>elist_"+name,preprocCutVar,"entrylist");
	elist = (TEntryList*)gDirectory->Get("elist_"+name);
	elist->SetReapplyCut(kTRUE);
	procNtuple->SetEntryList(elist);
	if(preprocCutVar == "1" || preprocCutVar == "1.0"){preprocCut=false;}else{preprocCut=true;}
	procEntries = (Double_t)elist->GetN();
	getWeightedEntries(new TString(""),&procWeightedEntries, &d_procWeightedEntries);
	//print();
}

cropdataset::cropdataset(TString _name, TString _filePath, TString _ntuplePath, TString _preprocCutVar="", TString _perEventWeightVar = "", UInt_t _color=0, UInt_t _fill=0){
	name = _name;
	filePath = _filePath;
	ntuplePath = _ntuplePath;
	signal = true;
	preprocCutVar = _preprocCutVar;
	fill = _fill;
	color = _color;
	if(_perEventWeightVar==""){
		perEventWeightVar="1";
	}else{
		perEventWeightVar = _perEventWeightVar;
	}
	originalTFile = TFile::Open(filePath,"READ");
	originalNtuple = (TTree*)originalTFile->Get(ntuplePath);
	totalEntries = (Double_t)originalNtuple->GetEntries();
	if(_perEventWeightVar=="A"){
		perEventWeightVar="1/";
		perEventWeightVar.Append(toString(totalEntries));
	}
	procTFile = originalTFile;
	procTFile = originalTFile;
	procNtuple = originalNtuple;
	getWeightedEntries(new TString(""),&totalWeightedEntries, &d_totalWeightedEntries);
	procNtuple->Draw(">>elist_"+name,preprocCutVar,"entrylist");
	elist = (TEntryList*)gDirectory->Get("elist_"+name);
	elist->SetReapplyCut(kTRUE);
	procNtuple->SetEntryList(elist);
	preprocCut=true;
	procEntries = (Double_t)elist->GetN();
	getWeightedEntries(new TString(""),&procWeightedEntries, &d_procWeightedEntries);
	//print();
}

cropdataset::cropdataset(TString _name, TString _filePath, TString _ntuplePath, bool _signal, TString _preprocCutVar="", TString _perEventWeightVar = ""){
	fill = 1001;
	color = 1;
	name = _name;
	filePath = _filePath;
	ntuplePath = _ntuplePath;
	signal = _signal;
	preprocCutVar = _preprocCutVar;
	if(_perEventWeightVar==""){
		perEventWeightVar="1";
	}else{
		perEventWeightVar = _perEventWeightVar;
	}
	cout <<"trying to open: " << filePath << endl;
	originalTFile = TFile::Open(filePath,"READ");
	cout << "SUCCESS" << endl;
	originalNtuple = (TTree*)originalTFile->Get(ntuplePath);
	totalEntries = (Double_t)originalNtuple->GetEntries();
	if(_perEventWeightVar=="A"){
		perEventWeightVar="1/";
		perEventWeightVar.Append(toString(totalEntries));
	}
	procTFile = originalTFile;
	procTFile = originalTFile;
	procNtuple = originalNtuple;
	getWeightedEntries(new TString(""),&totalWeightedEntries, &d_totalWeightedEntries);
	procNtuple->Draw(">>elist_"+name,preprocCutVar,"entrylist");
	elist = (TEntryList*)gDirectory->Get("elist_"+name);
	elist->SetReapplyCut(kTRUE);
	procNtuple->SetEntryList(elist);
	preprocCut=true;
	procEntries = (Double_t)elist->GetN();
	getWeightedEntries(new TString(""),&procWeightedEntries, &d_procWeightedEntries);
	//print();
}





cropdataset::cropdataset(TString line, TString weightfilename, UInt_t linenum){
	signal = true;
	name = "default";
	filePath ="default";
	ntuplePath="default";
	preprocCutVar ="1";
	perEventWeightVar="1";
	special=1.0;
	color =0;
	fill =0;
	TObjArray* Strings = line.Tokenize(cropdatasetdelimiters);
//	cout << "Strings: " << Strings->GetEntries() << endl;
	if(Strings->GetEntries() >=8){
		TIter iString(Strings);
		TObjString* os=(TObjString*)iString(); //S/B:

		TString word = os->GetString();
		//S/B file ntuple weight cut legend fill color special
		if(word == cropdatasetBstr){
			signal = false;
		}else{
			if(word==cropdatasetSstr){
				signal = true;
			}else{
				cout << "FATAL: " << weightfilename << " line " << linenum << ": First character must be either S or B!" << endl;
				exit(EXIT_FAILURE);
			}
		}
		os=(TObjString*)iString(); //Filename
		try{
			filePath = lexical_cast <string> (os->GetString());
		}catch(bad_lexical_cast &e){
			cout << "FATAL: " << weightfilename << " line " << linenum <<  "Invalid filename" << endl;
			exit(EXIT_FAILURE);
		}
		os=(TObjString*)iString(); //Ntuplename
		try{
			ntuplePath = lexical_cast <string> (os->GetString());
		}catch(bad_lexical_cast &e){
			cout << "FATAL: " << weightfilename << " line " << linenum <<  "Invalid ntuplename" << endl;
			exit(EXIT_FAILURE);

		}
		os=(TObjString*)iString(); //Weight
		try{
			perEventWeightVar = lexical_cast <string> (os->GetString());

		}catch(bad_lexical_cast &e){
			cout << "FATAL: " << weightfilename << " line " << linenum <<  "Invalid weight" << endl;
			exit(EXIT_FAILURE);
		}
		os=(TObjString*)iString(); //Cut
		try{
			preprocCutVar = lexical_cast <string> (os->GetString());
		}catch(bad_lexical_cast &e){
			cout << "FATAL: " << weightfilename << " line " << linenum <<  "Invalid preproc. cut" << endl;
			exit(EXIT_FAILURE);
		}

		os=(TObjString*)iString(); //Legend
		try{
			name = lexical_cast <string> (os->GetString());

		}catch(bad_lexical_cast &e){
			cout << "FATAL: " << weightfilename << " line " << linenum <<  "Invalid legend" << endl;
			exit(EXIT_FAILURE);

		}

		os=(TObjString*)iString(); //fill
		try{
			fill = lexical_cast <unsigned int> (os->GetString());
		}catch(bad_lexical_cast &e){
			cout << "FATAL: " << weightfilename << " line " << linenum <<  "Invalid fill style" << endl;
			exit(EXIT_FAILURE);
		}

		os=(TObjString*)iString(); //color
		try{
			color = lexical_cast <unsigned int> (os->GetString());
		}catch(bad_lexical_cast &e){
			cout << "FATAL: " << weightfilename << " line " << linenum <<  "Invalid color" << endl;
			exit(EXIT_FAILURE);
		}
		if(Strings->GetEntries() ==9){
			os=(TObjString*)iString(); //special
			try{
				special = lexical_cast <double> (os->GetString());
			}catch(bad_lexical_cast &e){
				cout << "FATAL: " << weightfilename << " line " << linenum <<  "Invalid special" << endl;
				exit(EXIT_FAILURE);
			}
		}

	}else{

		if(Strings->GetEntries() != 5 && Strings->GetEntries() != 4){
			cout << "FATAL: " << weightfilename << " line " << linenum <<  " expected 8 entries, found " << Strings->GetEntries() << ". Check the line is correct and try again." << endl;
			exit(EXIT_FAILURE);
		}else{
			TIter iString(Strings);
			TObjString* os=(TObjString*)iString();
			TString word = os->GetString();
			if(word == cropdatasetBstr){
				signal = false;
				name = "Background";
				name += linenum;
			}else{

				if(word==cropdatasetSstr){
					signal = true;
					name = "Signal";
					name += linenum;
				}else{
					cout << "FATAL: " << weightfilename << " line " << linenum << ": First character must be either S or B!" << endl;
					exit(EXIT_FAILURE);
				}
			}

			os=(TObjString*)iString(); //Filename
			try{
				filePath = lexical_cast <string> (os->GetString());
			}catch(bad_lexical_cast &e){
				cout << "FATAL: " << weightfilename << " line " << linenum <<  "Invalid filename" << endl;
				exit(EXIT_FAILURE);
			}
			os=(TObjString*)iString(); //Ntuplename
			try{
				ntuplePath = lexical_cast <string> (os->GetString());
			}catch(bad_lexical_cast &e){
				cout << "FATAL: " << weightfilename << " line " << linenum <<  "Invalid ntuplename" << endl;
				exit(EXIT_FAILURE);

			}
			os=(TObjString*)iString(); //Weight
			try{
				perEventWeightVar = lexical_cast <string> (os->GetString());

			}catch(bad_lexical_cast &e){
				cout << "FATAL: " << weightfilename << " line " << linenum <<  "Invalid weight" << endl;
				exit(EXIT_FAILURE);
			}
			if((os=(TObjString*)iString())){
				try{
					preprocCutVar = lexical_cast <string> (os->GetString());
				}catch(bad_lexical_cast &e){
					cout << "FATAL: " << weightfilename << " line " << linenum <<  "Invalid preproc. cut" << endl;
					exit(EXIT_FAILURE);
				}
			}else{
				preprocCutVar = "1";
			}

		}
	}
	originalTFile = TFile::Open(filePath,"READ");
	originalNtuple = (TTree*)originalTFile->Get(ntuplePath);
	totalEntries = (Double_t)originalNtuple->GetEntries();
	if(perEventWeightVar=="A"){
		perEventWeightVar="1/";
		perEventWeightVar.Append(toString(totalEntries));
	}
	procTFile = originalTFile;
	procNtuple = originalNtuple;
	getWeightedEntries(new TString(""),&totalWeightedEntries, &d_totalWeightedEntries);
	procNtuple->Draw(">>elist_"+name,preprocCutVar,"entrylist");
	elist = (TEntryList*)gDirectory->Get("elist_"+name);
	elist->SetReapplyCut(kTRUE);
	procNtuple->SetEntryList(elist);
	if(preprocCutVar == "1" || preprocCutVar == "1.0"){preprocCut=false;}else{preprocCut=true;}
	procEntries = (Double_t)elist->GetN();
	getWeightedEntries(new TString(""),&procWeightedEntries, &d_procWeightedEntries);
	//print();
}

void cropdataset::getWeightedEntries(TString *cut, Double_t *Entries=0, Double_t *d_Entries=0) const{
	TH1D *tmpWeightHisto;
	tmpWeightHisto = new TH1D("tmpWeightHisto","tmpWeightHisto",1,-3.0,3.0);
	tmpWeightHisto->Sumw2();
	if(*cut != ""){
		procNtuple->Draw(*cut+">>tmpWeightHisto","("+*cut+")*("+perEventWeightVar+")");
	}else{
		procNtuple->Draw("1>>tmpWeightHisto","("+perEventWeightVar+")");
	}
	*Entries = tmpWeightHisto->Integral();
	*d_Entries = sqrt(tmpWeightHisto->GetSumw2()->At(1));
	tmpWeightHisto->Delete();
}

void cropdataset::print() const{


	cout << "INFO: ==========================================================================" << endl;
	cout << "INFO: "<< name << " is a dataset of type "; if(signal){cout << "Signal " << endl; }else{ cout << "Background " << endl;}
	cout << "INFO: --------------------------------------------------------------------------" << endl;
	cout << "INFO:	File Path: " << filePath << endl;
	cout << "INFO:	Ntuple Path: " << ntuplePath << endl;
	cout << "INFO:	Total Entries: " << totalEntries <<endl;
	cout << "INFO:	Using per-event weight of " << perEventWeightVar << endl;
	cout << "INFO:	Total Weighted Entries: " << totalWeightedEntries << "+/-" << d_totalWeightedEntries << endl;
	cout << "INFO:	Special : " << special << endl;
	if(preprocCut){
		cout << "INFO:	Using Preprocessing cut: " << preprocCutVar << endl;
		cout << "INFO:	Entries After Preprocessing cuts: " << procEntries << endl;
		cout << "INFO:	Weighted Entries After Preprocessing cuts: " << procWeightedEntries << "+/-" << d_procWeightedEntries << endl;
	}
	cout << "INFO: ==========================================================================" << endl;
}

void cropdataset::printOneLine() const{
	Double_t eff = 0.0;
	Double_t d_eff = 0.0;
	Double_t weff = 0.0;
	Double_t d_weff = 0.0;

	binomEff(procEntries,totalEntries,&eff, &d_eff);
	binomEff(procWeightedEntries,totalWeightedEntries,&weff, &d_weff);

	if(signal){cout << "S\t";}else{cout << "B\t";}

	cout << prettyPrint(totalEntries) << "\t" << prettyPrint(procEntries) << "\t" << prettyPrint(eff) << "+/-" << prettyPrint(d_eff) << "\t" << prettyPrint(totalWeightedEntries) << "+/-" << prettyPrint(d_totalWeightedEntries) << "\t" << prettyPrint(procWeightedEntries) << "+/-" << prettyPrint(d_procWeightedEntries) << "\t" << prettyPrint(weff) << "+/-" << prettyPrint(d_weff) << "\t" << special << "\t" << name << endl;
}

void cropdataset::printOneLine(TString *cut1, TString *cut2) const{
	Double_t incleff, dincleff, excleff, dexcleff, inclcands, dinclcands, exclcands, dexclcands;
	getWeightedEntries(cut2, &exclcands, &dexclcands);
	TString cut = *cut1;
	if(cut != ""){cut.Append("&&"+*cut2);}else{cut = *cut2;}
	getWeightedEntries(&cut, &inclcands, &dinclcands);
	getEfficiency(new TString(""), cut2, &excleff, &dexcleff);
	getEfficiency(cut1, cut2, &incleff, &dincleff);
	if(signal){cout << "S\t";}else{cout << "B\t";}
	cout <<prettyPrint(exclcands) << "+/-" << prettyPrint(dexclcands) << "\t" << prettyPrint(excleff) << "+/-" << prettyPrint(dexcleff) << "\t" <<  prettyPrint(inclcands) << "+/-" << prettyPrint(dinclcands) << "\t" << prettyPrint(incleff) << "+/-" << prettyPrint(dincleff) << "\t" << special << "\t" << name << endl;
}

TString cropdataset::getName() const{
	return name;
}

void cropdataset::setName(TString _name){
	name = _name;
}
void cropdataset::getEfficiency(TString *cut1, TString *cut2, Double_t *eff, Double_t *d_eff)const{
	Double_t S1, dS1;
	TString cut = *cut1;
	getWeightedEntries(&cut, &S1,&dS1);
	if(S1 != 0.0){
		if(cut != ""){cut.Append("&&"+*cut2);}else{cut = *cut2;}
		Double_t S2, dS2;
		getWeightedEntries(&cut, &S2,&dS2);
		binomEff(S2,S1,eff,d_eff);
	}else{
		*eff = 0.0;
		*d_eff = 0.0;
	}
}

TString cropdataset::toLine()const{
	//S/B file ntuple weight cut legend fill color
	TString out;
	if(signal){out = "S";}else{out = "B";}
	out += cropdatasetTab;
	out += filePath;
	out += cropdatasetTab;
	out += ntuplePath;
	out += cropdatasetTab;
	out += perEventWeightVar;
	out += cropdatasetTab;
	out += preprocCutVar;
	out += cropdatasetTab;
	out += name;
	out += cropdatasetTab;
	out += fill;
	out += cropdatasetTab;
	out += color;
	return out;
}

TH1D * cropdataset::getHisto(TString cut, TString varString, UInt_t bins, Double_t min, Double_t max)const{
	TH1D *Histo = new TH1D("temp",name,bins,min,max);
	TString wvarCut = "(";
	wvarCut += cut;
	wvarCut += ")";
	Histo->Sumw2();
	if(cut!=""){
		wvarCut += "*";
	}
	wvarCut += perEventWeightVar;
	//	cout << wvarCut << endl;
	procNtuple->Draw(varString+">>temp",wvarCut);
	Histo->SetLineColor(color);
	Histo->SetMarkerColor(color);
	Histo->SetFillStyle(fill);
	Histo->SetFillColor(color);
	Histo->SetTitle(name);
	Histo->SetName(varString+"_"+name);
	return Histo;
}
TH1D * cropdataset::getHisto(TString varString, UInt_t bins, Double_t min, Double_t max)const{
	TH1D *Histo = new TH1D("temp",name,bins,min,max);
	Histo->Sumw2();
	//      cout << wvarCut << endl;
	procNtuple->Draw(varString+">>temp",perEventWeightVar);
	Histo->SetLineColor(color);
	Histo->SetMarkerColor(color);
	Histo->SetFillStyle(fill);
	Histo->SetFillColor(color);
	Histo->SetTitle(name);
	Histo->SetName(varString+"_"+name);
	return Histo;
}


TH1D * cropdataset::getHisto(cropvarensemble *vars, UInt_t n)const{
	return getHisto(vars->getVar(n), vars->getVarBins(n), vars->getVarMinVal(n), vars->getVarMaxVal(n));
}

TH1D * cropdataset::getHisto(TString cut, cropvarensemble *vars, UInt_t n)const{
	return getHisto(cut, vars->getVar(n), vars->getVarBins(n), vars->getVarMinVal(n), vars->getVarMaxVal(n));
}

TH2D * cropdataset::getHisto(cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const{
	return getHisto(varsX->getVar(nX), varsX->getVarBins(nX), varsX->getVarMinVal(nX), varsX->getVarMaxVal(nX),varsY->getVar(nY), varsY->getVarBins(nY), varsY->getVarMinVal(nY), varsY->getVarMaxVal(nY));
}
TH2D * cropdataset::getHisto(TString cut, cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const{
	return getHisto(cut, varsX->getVar(nX), varsX->getVarBins(nX), varsX->getVarMinVal(nX), varsX->getVarMaxVal(nX),varsY->getVar(nY), varsY->getVarBins(nY), varsY->getVarMinVal(nY), varsY->getVarMaxVal(nY));
}
TH2D * cropdataset::getHisto(TString varX, UInt_t binsX, Double_t minX, Double_t maxX, TString varY, UInt_t binsY, Double_t minY, Double_t maxY)const{
	TH2D* Histo = new TH2D("temp",name,binsX,minX,maxX,binsY,minY,maxY);
	Histo->Sumw2();
	procNtuple->Draw(varY+":"+varX+">>temp",perEventWeightVar);
	Histo->SetTitle(name);
	Histo->SetName(varY+":"+varX+"_"+name);
	return Histo;
}
TH2D * cropdataset::getHisto(TString cut, TString varX, UInt_t binsX, Double_t minX, Double_t maxX, TString varY, UInt_t binsY, Double_t minY, Double_t maxY)const{
	TString wvarCut = "(";
	wvarCut += cut;
	wvarCut += ")";
	if(cut!=""){
		wvarCut += "*";
	}
	wvarCut += perEventWeightVar;
	TH2D* Histo = new TH2D("temp",name,binsX,minX,maxX,binsY,minY,maxY);
	Histo->Sumw2();
	procNtuple->Draw(varY+":"+varX+">>temp",wvarCut);
	Histo->SetTitle(name);
	Histo->SetName(varY+":"+varX+"_"+name);
	return Histo;
}


void cropdataset::setCut(TString *cut){
	procNtuple->Draw(">>tmpelist_"+name,*cut + "&&" + preprocCutVar,"entrylist");
	tmpelist = (TEntryList*)gDirectory->Get("tmpelist_"+name);
	procNtuple->SetEntryList(tmpelist);
	tmpelist->SetReapplyCut(kTRUE);
}

void cropdataset::resetCut(){
	procNtuple->SetEntryList(elist);
	tmpelist->Delete();
}
