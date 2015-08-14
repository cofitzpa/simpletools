#include "cropdatastore.h"
using std::ifstream;
using std::ofstream;


cropdatastore::cropdatastore(){
	name = "";
	NSignalDatasets=0, NBackgroundDatasets=0;
	totalSignalEntries=0.0, procSignalEntries=0.0;
	totalBackgroundEntries=0.0, procBackgroundEntries=0.0;
	totalSignalWeightedEntries=0.0, d_totalSignalWeightedEntries=0.0;
	totalBackgroundWeightedEntries=0.0, d_totalBackgroundWeightedEntries=0.0;
	procSignalWeightedEntries=0.0, d_procSignalWeightedEntries=0.0;
	procBackgroundWeightedEntries=0.0, d_procBackgroundWeightedEntries=0.0;

}
cropdatastore::cropdatastore(TString _name){
	name = _name;
	NSignalDatasets=0, NBackgroundDatasets=0;
	totalSignalEntries=0.0, procSignalEntries=0.0;
	totalBackgroundEntries=0.0, procBackgroundEntries=0.0;
	totalSignalWeightedEntries=0.0, d_totalSignalWeightedEntries=0.0;
	totalBackgroundWeightedEntries=0.0, d_totalBackgroundWeightedEntries=0.0;
	procSignalWeightedEntries=0.0, d_procSignalWeightedEntries=0.0;
	procBackgroundWeightedEntries=0.0, d_procBackgroundWeightedEntries=0.0;


}

cropdatastore::cropdatastore(TString _name, TString _filename){
	name = _name;
	NSignalDatasets=0, NBackgroundDatasets=0;
	totalSignalEntries=0.0, procSignalEntries=0.0;
	totalBackgroundEntries=0.0, procBackgroundEntries=0.0;
	totalSignalWeightedEntries=0.0, d_totalSignalWeightedEntries=0.0;
	totalBackgroundWeightedEntries=0.0, d_totalBackgroundWeightedEntries=0.0;
	procSignalWeightedEntries=0.0, d_procSignalWeightedEntries=0.0;
	procBackgroundWeightedEntries=0.0, d_procBackgroundWeightedEntries=0.0;

	UInt_t linenum = 0;
	ifstream fileStream2;
	fileStream2.open(_filename);
	if (!fileStream2){
		cout << "FATAL: Error in opening weightfile" << _filename << endl;
		exit(EXIT_FAILURE);
	}
	TString line;
	while(!line.ReadLine(fileStream2).eof()){
		linenum++;
		if(!line.BeginsWith("#")){
			this->addDataset(new cropdataset(line, _filename,linenum));

		}
	}
	fileStream2.close();

}
void cropdatastore::addDataset(cropdataset *dataset){
	if (dataset->isSignal()){
		SignalDatasets.push_back(dataset);
		NSignalDatasets++;
		totalSignalEntries += dataset->getTotalEntries();
		procSignalEntries += dataset->getProcEntries();
		totalSignalWeightedEntries += dataset->getTotalWeightedEntries();
		d_totalSignalWeightedEntries = sqrt(d_totalSignalWeightedEntries*d_totalSignalWeightedEntries + dataset->getTotalWeightedEntriesError()*dataset->getTotalWeightedEntriesError());

		procSignalWeightedEntries += dataset->getProcWeightedEntries();
		d_procSignalWeightedEntries = sqrt(d_procSignalWeightedEntries*d_procSignalWeightedEntries + dataset->getProcWeightedEntriesError()*dataset->getProcWeightedEntriesError());
	}else{
		BackgroundDatasets.push_back(dataset);
		NBackgroundDatasets++;
		totalBackgroundEntries += dataset->getTotalEntries();
		procBackgroundEntries += dataset->getProcEntries();
		totalBackgroundWeightedEntries += dataset->getTotalWeightedEntries();
		d_totalBackgroundWeightedEntries = sqrt(d_totalBackgroundWeightedEntries*d_totalBackgroundWeightedEntries + dataset->getTotalWeightedEntriesError()*dataset->getTotalWeightedEntriesError());
		procBackgroundWeightedEntries += dataset->getProcWeightedEntries();
		d_procBackgroundWeightedEntries = sqrt(d_procBackgroundWeightedEntries*d_procBackgroundWeightedEntries + dataset->getProcWeightedEntriesError()*dataset->getProcWeightedEntriesError());
	}

}

void cropdatastore::print() const{
	cout << "INFO:	========================================================================== " << endl;
	cout << "INFO: "<< name << " is a datastore containing " << NSignalDatasets << " signal and " << NBackgroundDatasets << " background datasets" << endl;
	cout << "INFO:	-------------------------------------------------------------------------- " << endl;
	cout << "INFO:	SIGNAL totals in this datastore:" << endl;
	cout << "INFO:  Total Entries: " << prettyPrint(totalSignalEntries) <<endl;
	cout << "INFO:	Weighted Entries: " << prettyPrint(totalSignalWeightedEntries) << "+/-" << prettyPrint(d_totalSignalWeightedEntries) << endl; 
	if(totalSignalWeightedEntries != procSignalWeightedEntries){
		cout << "INFO:  Entries after Preproc. Cuts: " << prettyPrint(procSignalEntries) <<endl;
		cout << "INFO:	Weighted Entries after Preproc. Cuts: " << prettyPrint(procSignalWeightedEntries) << "+/-" << prettyPrint(d_procSignalWeightedEntries) << endl; 		}
		cout << "INFO:	-------------------------------------------------------------------------- " << endl;
		cout << "INFO:  BACKGROUND totals in this datastore:" << endl;
		cout << "INFO:  Total Entries: " << prettyPrint(totalBackgroundEntries) <<endl;
		cout << "INFO:  Weighted Entries: " << prettyPrint(totalBackgroundWeightedEntries) << "+/-" << prettyPrint(d_totalBackgroundWeightedEntries) << endl; 
		if(totalBackgroundWeightedEntries!=procBackgroundWeightedEntries){
			cout << "INFO:  Entries after Preproc. Cuts: " << prettyPrint(procBackgroundEntries) <<endl;
			cout << "INFO:  Weighted Entries after Preproc. Cuts: " << prettyPrint(procBackgroundWeightedEntries) << "+/-" << prettyPrint(d_procBackgroundWeightedEntries) << endl; }
			cout << "INFO:	========================================================================== " << endl;
}


cropdataset* cropdatastore::getSignalDataset(UInt_t sig) const{
	return SignalDatasets[sig];
}
cropdataset* cropdatastore::getBackgroundDataset(UInt_t bkg) const{
	return BackgroundDatasets[bkg];
}

cropdataset* cropdatastore::getDataset(TString datasetName) const{
	for(UInt_t i =0; i<NSignalDatasets; i++){
		if(SignalDatasets[i]->getName()==datasetName){return SignalDatasets[i];}
	}
	for(UInt_t i =0; i<NBackgroundDatasets; i++){
		if(BackgroundDatasets[i]->getName()==datasetName){return BackgroundDatasets[i];}
	}
	cout << "FATAL:	Could not find a dataset with name " << datasetName << " in datastore "<< name <<"!!!" << endl;
	exit(EXIT_FAILURE);
}

void cropdatastore::getWeightedSignalEntries(TString *cut, Double_t *Entries, Double_t *d_Entries) const{
	Double_t tmpEntries;
	Double_t d_tmpEntries;
	*Entries = 0;
	*d_Entries = 0;
	for(UInt_t i =0; i<NSignalDatasets; i++){
		SignalDatasets[i]->getWeightedEntries(cut, &tmpEntries, &d_tmpEntries);
		*Entries += tmpEntries;
		*d_Entries = sqrt(*d_Entries**d_Entries + d_tmpEntries*d_tmpEntries);	
	}
}

void cropdatastore::getWeightedBackgroundEntries(TString *cut, Double_t *Entries, Double_t *d_Entries) const{
	Double_t tmpEntries;
	Double_t d_tmpEntries;
	*Entries = 0;
	*d_Entries = 0;
	for(UInt_t i =0; i<NBackgroundDatasets; i++){
		BackgroundDatasets[i]->getWeightedEntries(cut, &tmpEntries, &d_tmpEntries);
		*Entries += tmpEntries;
		*d_Entries = sqrt(*d_Entries**d_Entries + d_tmpEntries*d_tmpEntries);	
	}
}

void cropdatastore::initStats() const{
	cout << "Type	Total Entries	Total Proc. Entries	Proc. Cut Eff.	Total Weighted Entries	Total Proc. Weighted Entries	Proc Weight ratio	Name"<< endl;
	cout << "-------------------------------------------------------------------"<< endl;
	for(UInt_t i =0; i<NSignalDatasets; i++){
		SignalDatasets[i]->printOneLine();
	}
	printOneLine(true);
	for(UInt_t i =0; i<NBackgroundDatasets; i++){
		BackgroundDatasets[i]->printOneLine();
	}
	printOneLine(false);
}
void cropdatastore::finalStats(TString *cut1, TString *cut2)const{
	cout << *cut2 << ":" << endl;
	cout << "Type\tExcl. Cands.\tExcl. Eff.\tIncl. Cands.\tIncl. Eff.\tName" << endl;
	for(UInt_t i =0; i<NSignalDatasets; i++){
		SignalDatasets[i]->printOneLine(cut1, cut2);
	}
	printOneLine(true, cut1, cut2);
	for(UInt_t i =0; i<NBackgroundDatasets; i++){
		BackgroundDatasets[i]->printOneLine(cut1,cut2);
	}
	printOneLine(false, cut1, cut2);
}


void cropdatastore::printOneLine(bool _signal, TString *cut1, TString *cut2) const{
	Double_t incleff, dincleff, excleff, dexcleff, inclcands, dinclcands,exclcands,dexclcands;
	TString cut = *cut1;
       	if(cut != ""){cut.Append("&&"+*cut2);}else{cut = *cut2;}
	if(_signal){
		getWeightedSignalEntries(cut2, &exclcands, &dexclcands);
	        getWeightedSignalEntries(&cut, &inclcands, &dinclcands);
		getSignalEfficiency(cut1, cut2, &incleff, &dincleff);
		getSignalEfficiency(new TString(""), cut2, &excleff, &dexcleff);
		cout << "S\t" << prettyPrint(exclcands) << "+/-" << prettyPrint(dexclcands) << "\t" << prettyPrint(excleff) << "+/-" << prettyPrint(dexcleff) << "\t"<<  prettyPrint(inclcands) << "+/-" << prettyPrint(dinclcands) << "\t" <<  prettyPrint(incleff) << "+/-" << prettyPrint(dincleff) << "	TOTAL" << endl;
	}else{
		getWeightedBackgroundEntries(cut2, &exclcands, &dexclcands);
                getWeightedBackgroundEntries(&cut, &inclcands, &dinclcands);
		getBackgroundEfficiency(cut1, cut2, &incleff, &dincleff);
		getBackgroundEfficiency(new TString(""), cut2, &excleff, &dexcleff);
		cout << "B\t" << prettyPrint(exclcands) << "+/-" << prettyPrint(dexclcands) << "\t" << prettyPrint(excleff) << "+/-" << prettyPrint(dexcleff) << "\t"<<  prettyPrint(inclcands) << "+/-" << prettyPrint(dinclcands) << "\t" <<  prettyPrint(incleff) << "+/-" << prettyPrint(dincleff) << "	TOTAL" << endl;
	}
}

void cropdatastore::printOneLine(bool _signal) const{
	Double_t eff = 0.0;
	Double_t d_eff = 0.0;
	Double_t weff = 0.0;
	Double_t d_weff = 0.0;
	if(_signal){
		binomEff(procSignalEntries,totalSignalEntries,&eff, &d_eff);

		binomEff(procSignalWeightedEntries,totalSignalWeightedEntries,&weff, &d_weff);
		cout << "S\t" << prettyPrint(totalSignalEntries) << "\t" << prettyPrint(procSignalEntries) << "\t" << prettyPrint(eff) << "+/-" << prettyPrint(d_eff) << "\t" << prettyPrint(totalSignalWeightedEntries) << "+/-" << prettyPrint(d_totalSignalWeightedEntries) << "\t" << prettyPrint(procSignalWeightedEntries) << "+/-" << prettyPrint(d_procSignalWeightedEntries) << "\t" << prettyPrint(weff) << "+/-" << prettyPrint(d_weff) << "\tTOTAL" << endl;
	}else{
		binomEff(procBackgroundEntries,totalBackgroundEntries,&eff, &d_eff);

		binomEff(procBackgroundWeightedEntries,totalBackgroundWeightedEntries,&weff, &d_weff);
		cout << "B\t" << prettyPrint(totalBackgroundEntries) << "\t" << prettyPrint(procBackgroundEntries) << "\t" << prettyPrint(eff) << "+/-" << prettyPrint(d_eff) << "\t" << prettyPrint(totalBackgroundWeightedEntries) << "+/-" << prettyPrint(d_totalBackgroundWeightedEntries) << "\t" << prettyPrint(procBackgroundWeightedEntries) << "+/-" << prettyPrint(d_procBackgroundWeightedEntries) << "\t" << prettyPrint(weff) << "+/-" << prettyPrint(d_weff) << "\tTOTAL" << endl;

	}
}

void cropdatastore::getSignalEfficiency(TString *cut1, TString *cut2, Double_t *eff, Double_t *d_eff)const{
	Double_t S1, dS1;
	getWeightedSignalEntries(cut1, &S1,&dS1);
	if(S1 != 0.0){
		TString cut = *cut1;
		if(*cut1 != ""){cut.Append("&&"+*cut2);}else{cut = *cut2;}
		Double_t S2, dS2;
		getWeightedSignalEntries(&cut, &S2,&dS2);
		binomEff(S2,S1,eff,d_eff);
	}else{
		*eff = 0.0;
		*d_eff = 0.0;
	}
}

void cropdatastore::getBackgroundEfficiency(TString *cut1, TString *cut2, Double_t *eff, Double_t *d_eff)const{
	Double_t S1, dS1;
	getWeightedBackgroundEntries(cut1, &S1,&dS1);
	if(S1 != 0.0){
		TString cut = *cut1;
		if(*cut1 != ""){cut.Append("&&"+*cut2);}else{cut = *cut2;}
		Double_t S2, dS2;
		getWeightedBackgroundEntries(&cut, &S2,&dS2);
		binomEff(S2,S1,eff,d_eff);
	}else{
		*eff = 0.0;
		*d_eff = 0.0;
	}
}

TH2D * cropdatastore::getSignalHisto(cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const{
	return getSignalHisto(varsX->getVar(nX), varsX->getVarBins(nX), varsX->getVarMinVal(nX), varsX->getVarMaxVal(nX),varsY->getVar(nY), varsY->getVarBins(nY), varsY->getVarMinVal(nY), varsY->getVarMaxVal(nY));
}
TH2D * cropdatastore::getSignalHisto(TString cut, cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const{
	return getSignalHisto(cut, varsX->getVar(nX), varsX->getVarBins(nX), varsX->getVarMinVal(nX), varsX->getVarMaxVal(nX),varsY->getVar(nY), varsY->getVarBins(nY), varsY->getVarMinVal(nY), varsY->getVarMaxVal(nY));
}
TH2D * cropdatastore::getSignalHisto(TString varX, UInt_t binsX, Double_t minX, Double_t maxX,TString varY, UInt_t binsY, Double_t minY, Double_t maxY)const{
	TH2D* Histo = new TH2D("signal",name,binsX,minX,maxX,binsY,minY,maxY);
	Histo->Sumw2();
	for(UInt_t i =0; i<NSignalDatasets; i++){
		TH2D *tmpHisto = SignalDatasets[i]->getHisto(varX, binsX, minX, maxX, varY, binsY, minY, maxY);	                
		Histo->Add(tmpHisto);
		delete tmpHisto;
	}
	return Histo;

}
TH2D * cropdatastore::getSignalHisto(TString cut, TString varX, UInt_t binsX, Double_t minX, Double_t maxX,TString varY, UInt_t binsY, Double_t minY, Double_t maxY)const{
	TH2D* Histo = new TH2D("signal",name,binsX,minX,maxX,binsY,minY,maxY);
	Histo->Sumw2();
	for(UInt_t i =0; i<NSignalDatasets; i++){
		TH2D *tmpHisto = SignalDatasets[i]->getHisto(cut, varX, binsX, minX, maxX, varY, binsY, minY, maxY);
		Histo->Add(tmpHisto);
		delete tmpHisto;
	}
	return Histo;


}
TH2D * cropdatastore::getBackgroundHisto(cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const{
	return getBackgroundHisto(varsX->getVar(nX), varsX->getVarBins(nX), varsX->getVarMinVal(nX), varsX->getVarMaxVal(nX),varsY->getVar(nY), varsY->getVarBins(nY), varsY->getVarMinVal(nY), varsY->getVarMaxVal(nY));
}
TH2D * cropdatastore::getBackgroundHisto(TString cut, cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const{
	return getBackgroundHisto(cut, varsX->getVar(nX), varsX->getVarBins(nX), varsX->getVarMinVal(nX), varsX->getVarMaxVal(nX),varsY->getVar(nY), varsY->getVarBins(nY), varsY->getVarMinVal(nY), varsY->getVarMaxVal(nY));

}
TH2D * cropdatastore::getBackgroundHisto(TString varX, UInt_t binsX, Double_t minX, Double_t maxX,TString varY, UInt_t binsY, Double_t minY, Double_t maxY)const{
	TH2D* Histo = new TH2D("background",name,binsX,minX,maxX,binsY,minY,maxY);
	Histo->Sumw2();
	for(UInt_t i =0; i<NBackgroundDatasets; i++){
		TH2D *tmpHisto = BackgroundDatasets[i]->getHisto(varX, binsX, minX, maxX, varY, binsY, minY, maxY);
		Histo->Add(tmpHisto);
		delete tmpHisto;
	}
	return Histo;


}
TH2D * cropdatastore::getBackgroundHisto(TString cut, TString varX, UInt_t binsX, Double_t minX, Double_t maxX,TString varY, UInt_t binsY, Double_t minY, Double_t maxY)const{
	TH2D* Histo = new TH2D("background",name,binsX,minX,maxX,binsY,minY,maxY);
	Histo->Sumw2();
	for(UInt_t i =0; i<NBackgroundDatasets; i++){
		TH2D *tmpHisto = BackgroundDatasets[i]->getHisto(cut, varX, binsX, minX, maxX, varY, binsY, minY, maxY);
		Histo->Add(tmpHisto);
		delete tmpHisto;
	}
	return Histo;

}

TH2D * cropdatastore::getHisto(cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const{
return getHisto(varsX->getVar(nX), varsX->getVarBins(nX), varsX->getVarMinVal(nX), varsX->getVarMaxVal(nX),varsY->getVar(nY), varsY->getVarBins(nY), varsY->getVarMinVal(nY), varsY->getVarMaxVal(nY));
}

TH2D * cropdatastore::getHisto(TString cut, cropvarensemble *varsX, cropvarensemble *varsY, UInt_t nX, UInt_t nY)const{
return getHisto(cut, varsX->getVar(nX), varsX->getVarBins(nX), varsX->getVarMinVal(nX), varsX->getVarMaxVal(nX),varsY->getVar(nY), varsY->getVarBins(nY), varsY->getVarMinVal(nY), varsY->getVarMaxVal(nY));
}
TH2D * cropdatastore::getHisto(TString varX, UInt_t binsX, Double_t minX, Double_t maxX,TString varY, UInt_t binsY, Double_t minY, Double_t maxY)const{
	TH2D* Histo = new TH2D(varX+":"+varY+"_"+name,varX+":"+varY+"_"+name,binsX,minX,maxX,binsY,minY,maxY);  
	Histo->Sumw2();
	TH2D *signal = getSignalHisto(varX,binsX,minX,maxX,varY,binsY,minY,maxY);
	Histo->Add(signal);
	TH2D* background = getBackgroundHisto(varX,binsX,minX,maxX,varY,binsY,minY,maxY);
	Histo->Add(background);
	delete signal;
	delete background;
	return Histo;

}
TH2D * cropdatastore::getHisto(TString cut, TString varX, UInt_t binsX, Double_t minX, Double_t maxX,TString varY, UInt_t binsY, Double_t minY, Double_t maxY)const{
	TH2D* Histo = new TH2D(varX+":"+varY+"_"+name,varX+":"+varY+"_"+name,binsX,minX,maxX,binsY,minY,maxY);
	Histo->Sumw2();
	TH2D *signal = getSignalHisto(cut,varX,binsX,minX,maxX,varY,binsY,minY,maxY);
	Histo->Add(signal);
	TH2D* background = getBackgroundHisto(cut,varX,binsX,minX,maxX,varY,binsY,minY,maxY);
	Histo->Add(background);
	delete signal;
	delete background;
	return Histo;

}



TH1D* cropdatastore::getSignalHisto(TString cut, TString varString, UInt_t bins, Double_t min, Double_t max)const{
	TH1D * Histo = new TH1D("signal",name,bins,min,max);
	Histo->Sumw2();
	for(UInt_t i =0; i<NSignalDatasets; i++){
		TH1D * tmpHisto = SignalDatasets[i]->getHisto(cut, varString, bins, min, max);
		Histo->Add(tmpHisto);

		delete tmpHisto;
	}
	return Histo;
}

TH1D* cropdatastore::getSignalHisto(TString varString, UInt_t bins, Double_t min, Double_t max)const{
	TH1D * Histo = new TH1D("signal",name,bins,min,max);
	Histo->Sumw2();
	for(UInt_t i =0; i<NSignalDatasets; i++){
		TH1D * tmpHisto = SignalDatasets[i]->getHisto(varString, bins, min, max);
		Histo->Add(tmpHisto);
		delete tmpHisto;
	}
	return Histo;
}

TH1D* cropdatastore::getBackgroundHisto(TString cut, TString varString, UInt_t bins, Double_t min, Double_t max)const{
	TH1D * Histo = new TH1D("background",name,bins,min,max);
	Histo->Sumw2();
	for(UInt_t i =0; i<NBackgroundDatasets; i++){
		TH1D * tmpHisto = BackgroundDatasets[i]->getHisto(cut, varString, bins, min, max);
		Histo->Add(tmpHisto);
		tmpHisto->Delete();
	}
	return Histo;
}

TH1D* cropdatastore::getBackgroundHisto(TString varString, UInt_t bins, Double_t min, Double_t max)const{
	TH1D * Histo = new TH1D("background",name,bins,min,max);
	Histo->Sumw2();
	for(UInt_t i =0; i<NBackgroundDatasets; i++){
		TH1D * tmpHisto = BackgroundDatasets[i]->getHisto(varString, bins, min, max);
		Histo->Add(tmpHisto);
		tmpHisto->Delete();
	}
	return Histo;
}

TH1D* cropdatastore::getHisto(TString cut, TString varString, UInt_t bins, Double_t min, Double_t max)const{
	TH1D * Histo = new TH1D(varString + "_" + name,varString + "_" +name,bins,min,max);
	Histo->Sumw2();
	Histo->Add(getSignalHisto(cut,varString, bins, min, max));
	Histo->Add(getBackgroundHisto(cut,varString, bins, min, max));
	return Histo;
}

TH1D* cropdatastore::getHisto(TString varString, UInt_t bins, Double_t min, Double_t max)const{
	TH1D * Histo = new TH1D(varString + "_" +name,varString + "_" +name,bins,min,max);
	Histo->Sumw2();
	TH1D *signal = getSignalHisto(varString, bins, min, max);
	Histo->Add(signal);
	TH1D* background = getBackgroundHisto(varString, bins, min, max);
	Histo->Add(background);
	delete signal;
	delete background;
	return Histo;
}

TH1D* cropdatastore::getHisto(cropvarensemble *_data, UInt_t n)const{
	return getHisto(_data->getVar(n),_data->getVarBins(n),_data->getVarMinVal(n),_data->getVarMaxVal(n));
}

TH1D* cropdatastore::getHisto(TString cut, cropvarensemble *_data, UInt_t n)const{
	return getHisto(cut, _data->getVar(n),_data->getVarBins(n),_data->getVarMinVal(n),_data->getVarMaxVal(n));
}

TH1D* cropdatastore::getBackgroundHisto(cropvarensemble *_data, UInt_t n)const{
	return getBackgroundHisto(_data->getVar(n),_data->getVarBins(n),_data->getVarMinVal(n),_data->getVarMaxVal(n));
}
TH1D* cropdatastore::getBackgroundHisto(TString cut, cropvarensemble *_data, UInt_t n)const{
	return getBackgroundHisto(cut, _data->getVar(n),_data->getVarBins(n),_data->getVarMinVal(n),_data->getVarMaxVal(n));
}


TH1D* cropdatastore::getSignalHisto(cropvarensemble *_data, UInt_t n)const{
	return getSignalHisto(_data->getVar(n),_data->getVarBins(n),_data->getVarMinVal(n),_data->getVarMaxVal(n));
}
TH1D* cropdatastore::getSignalHisto(TString cut, cropvarensemble *_data, UInt_t n)const{
	return getSignalHisto(cut, _data->getVar(n),_data->getVarBins(n),_data->getVarMinVal(n),_data->getVarMaxVal(n));
} 

THStack* cropdatastore::getStack(cropvarensemble *_data, UInt_t n)const{
	return getStack(_data->getVar(n),_data->getVarBins(n),_data->getVarMinVal(n),_data->getVarMaxVal(n));
}
THStack* cropdatastore::getStack(TString cut, cropvarensemble *_data, UInt_t n)const{
	return getStack(cut, _data->getVar(n),_data->getVarBins(n),_data->getVarMinVal(n),_data->getVarMaxVal(n));
}

THStack * cropdatastore::getStack(TString varString, UInt_t bins, Double_t min, Double_t max)const{
	THStack *hs = new THStack(varString,varString);
	TH1D *tmpHisto;
	for(UInt_t i=0; i<NBackgroundDatasets; i++){
		tmpHisto = BackgroundDatasets[i]->getHisto(varString, bins, min, max);
		hs->Add(tmpHisto);
	}

	for(UInt_t i=0; i<NSignalDatasets; i++){
		tmpHisto = SignalDatasets[i]->getHisto(varString, bins, min, max);
		hs->Add(tmpHisto);
	}
	return hs;
}


THStack * cropdatastore::getStack(TString cut, TString varString, UInt_t bins, Double_t min, Double_t max)const{
	THStack *hs = new THStack(varString,varString);
	TH1D *tmpHisto;
	for(UInt_t i=0; i<NBackgroundDatasets; i++){
		tmpHisto = BackgroundDatasets[i]->getHisto(cut, varString, bins, min, max);
		hs->Add(tmpHisto);
	}
	for(UInt_t i=0; i<NSignalDatasets; i++){
		tmpHisto = SignalDatasets[i]->getHisto(cut, varString, bins, min, max);
		hs->Add(tmpHisto);
	}
	return hs;
}
void cropdatastore::writeToFile(TString _filename)const{
	ofstream fileStream;
	fileStream.open(_filename);
	fileStream << "#S/B\tfile\tntuple\tweight\tcut\tlegend\tfill\tcolor\n";
	for(UInt_t i=0; i<NSignalDatasets; i++){

		fileStream << SignalDatasets[i]->toLine() << "\n";
	}
	for(UInt_t i=0; i<NBackgroundDatasets; i++){

		fileStream << BackgroundDatasets[i]->toLine() << "\n";
	}
	fileStream.close();
}

void cropdatastore::setCut(TString *cut){
	for(UInt_t i =0; i<NSignalDatasets; i++){
		SignalDatasets[i]->setCut(cut);
	}
	for(UInt_t i =0; i<NBackgroundDatasets; i++){
		BackgroundDatasets[i]->setCut(cut);
	}
}

void cropdatastore::resetCut(){
	for(UInt_t i =0; i<NSignalDatasets; i++){
		SignalDatasets[i]->resetCut();
	}
	for(UInt_t i =0; i<NBackgroundDatasets; i++){
		BackgroundDatasets[i]->resetCut();
	}

}
