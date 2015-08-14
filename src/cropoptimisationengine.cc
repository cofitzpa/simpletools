#include "cropoptimisationengine.h"

cropoptimisationengine::cropoptimisationengine(cropdatastore *_datastore, cropcutensemble  *_cutensemble, TCanvas *_c, TRandom3 *_rng, void(*_fom)(cropdatastore*, TString*, Double_t*, Double_t*), TFunctor *_OrderMethod){
	
	datastore = _datastore;
	cutensemble = _cutensemble;
	fom = _fom;
	c = _c;
	rng = _rng;
	OrderMethod = _OrderMethod;
}

void cropoptimisationengine::initialise(){
	Double_t S = 0;
	Double_t d_S = 0;
	Double_t B = 0;
	Double_t d_B = 0;
	Double_t Seff = 0;
	Double_t d_Seff = 0;
	Double_t Brej = 0;
	Double_t d_Brej = 0;
	Double_t FoM = 0;
	Double_t d_FoM = 0;

	Double_t initS;
	Double_t initB;
	Double_t d_initS; 
	Double_t d_initB; 
	TString cut;
	
	datastore->getWeightedSignalEntries(new TString(""),&initS, &d_initS);
	datastore->getWeightedBackgroundEntries(new TString(""),&initB, &d_initB);
	cout << "order  number  Sig. Eff        Bkg. Rej        FoM     cut"    << endl;
	cout << "-----------------------------------------------------------"   << endl;
	for(UInt_t i = 0; i<cutensemble->NCutVars; i++){
		for(UInt_t j = 0; j<cutensemble->getCutSteps(i); j++){
			cutensemble->getCut(i,j,&cut);
			datastore->getWeightedSignalEntries(&cut,&S,&d_S);
			datastore->getWeightedBackgroundEntries(&cut,&B,&d_B);
			binomEff(S,initS,&Seff,&d_Seff);
			binomEff(B,initB,&Brej,&d_Brej);
			Brej = 1.0 - Brej;
			fom(datastore,&cut,&FoM,&d_FoM);
			cutensemble->addCutSpacePoint(i,j,&S,&d_S,&B,&d_B,&Seff,&d_Seff,&Brej,&d_Brej,&FoM,&d_FoM);
		}
		cutensemble->setOptimalStep(i);
		cutensemble->printOneLine(i);
		cutensemble->plotCut(i,c);
	}


	//if(!cutensemble->OrderAscendingSeffBeff()){
	if(!(OrderMethod->Call())){
	cutensemble->OrderRandom();
	}

	//cutensemble->print();
	cutensemble->writeAllPlots(c);
}

bool cropoptimisationengine::optimise(UInt_t maxsteps){
	UInt_t step = 0;
	Double_t S = 0;
	Double_t d_S = 0;
	Double_t B = 0;
	Double_t d_B = 0;
	Double_t Seff = 0;
	Double_t d_Seff = 0;
	Double_t Brej = 0;
	Double_t d_Brej = 0;
	Double_t FoM = 0;
	Double_t d_FoM = 0;
	Double_t initS;
	Double_t initB;
	Double_t d_initS; 
	Double_t d_initB; 

	TString cut;
	TString inclcut;
	TString stepcut;
	Double_t stepnum = 0.0;
	bool converged = false;
	bool thisconverged = false;
	while(!converged && step < maxsteps){
		bool allconverged = true;
		cout << "INFO: STEP NUMBER: " << step << " OF " << maxsteps << ":" << endl;
		cout << "convergency?	order  number  Sig. Eff        Bkg. Rej        FoM     cut"    << endl;
		cout << "---------------------------------------------------------------------------------------"   << endl;
		for(UInt_t i = 0; i<cutensemble->NCutVars; i++){
			cutensemble->getOptimalEnsemble(i,&inclcut);
			datastore->getWeightedSignalEntries(&inclcut,&initS, &d_initS);
			datastore->getWeightedBackgroundEntries(&inclcut,&initB, &d_initB);
			datastore->setCut(&inclcut);


//			inclcut.Append("&&");
			for(UInt_t j = 0; j<cutensemble->getCutSteps(i); j++){
				cutensemble->getCut(i,j,&stepcut);
				cut = stepcut;
//				cout << inclcut << "\t" << stepcut << endl;
				datastore->getWeightedSignalEntries(&cut,&S,&d_S);
				datastore->getWeightedBackgroundEntries(&cut,&B,&d_B);
				binomEff(S,initS,&Seff,&d_Seff);
				binomEff(B,initB,&Brej,&d_Brej);
				Brej = 1.0 - Brej;
				fom(datastore,&cut,&FoM,&d_FoM);
				cutensemble->addCutSpacePoint(i,j,&S,&d_S,&B,&d_B,&Seff,&d_Seff,&Brej,&d_Brej,&FoM,&d_FoM);
			}
			UInt_t lastOptimalStep = cutensemble->getOptimalStep(i);
			UInt_t thisOptimalStep = cutensemble->setOptimalStep(i);
			Step.push_back(stepnum);
			stepnum += 1.0;
			cutensemble->getOptimalCut(i,&cut);
			cutAtStep.push_back(cut);
			fomAtStep.push_back(cutensemble->getMaxFoM(i));
			d_fomAtStep.push_back(cutensemble->getdMaxFoM(i));
			if(lastOptimalStep == thisOptimalStep){thisconverged = true;}else{thisconverged=false;}
			cout << thisconverged << "\t\t";
			cutensemble->printOneLine(i);
			allconverged = thisconverged*allconverged;
			//cout << "lastOptimalStep: " << lastOptimalStep << "	thisOptimalStep: " << thisOptimalStep <<"	allconverged: " << allconverged << "	thisconverged: " << thisconverged << endl;
			cutensemble->plotCut(i,c);
			cutensemble->getCutVar(i,&cut);
			c->SetName(cut);
			c->Update();
			datastore->resetCut();
		}
		converged = allconverged;

		//if(!cutensemble->OrderAscendingSeffBeff()&&!converged){
		if(!(OrderMethod->Call())&&!converged){
		cout << "WARNING: Sanity check on ordering failed. Using random ordering for this step" << endl;
		cutensemble->OrderRandom();
		}
		step ++;
	}

	cutensemble->writeAllPlots(c);
	return converged;
}

void cropoptimisationengine::writeEvolutionPlot(){
	c->Clear();
	c->cd(0);
	c->SetName("FoM evolution");
	UInt_t size = fomAtStep.size();
	Double_t* fomArr = new Double_t [size];
	copy(fomAtStep.begin(),fomAtStep.end(),fomArr);
	Double_t* d_fomArr = new Double_t [size];
	copy(d_fomAtStep.begin(),d_fomAtStep.end(),d_fomArr);
	TString* cutArr = new TString [size];
	copy(cutAtStep.begin(),cutAtStep.end(),cutArr);
	Double_t* stepArr = new Double_t [size];
	copy(Step.begin(),Step.end(),stepArr);
	TGraphErrors * fomStepsGraph = new TGraphErrors(size,stepArr,fomArr,0,d_fomArr);
	fomStepsGraph->SetLineWidth(2);
	fomStepsGraph->Draw("AL1");
	fomStepsGraph->SetTitle("FoM evolution WRT number of reoptimisations");
	fomStepsGraph->GetYaxis()->SetTitle("FoM");
	for(UInt_t i = 0; i<size; i++){
		fomStepsGraph->GetXaxis()->SetBinLabel(fomStepsGraph->GetXaxis()->FindBin(i),cutArr[i]);

	}
	c->Update();
	c->Write();
}

void cropoptimisationengine::finalStats(){
	TString inclcut;
	TString exclcut;
	TString finalcut;
	cutensemble->OrderDescendingBrej();
	for(UInt_t i = 0; i<cutensemble->NCutVars; i++){
		cutensemble->getOptimalEnsemble(i,&inclcut);
		cout << endl;
		cutensemble->getOptimalCut(i,&exclcut);
		datastore->finalStats(&inclcut, &exclcut);
	}
	cutensemble->getOptimalEnsemble(cutensemble->NCutVars+1, &finalcut);
	cout << endl;
	cout << "Optimal cutstring: " << endl;
	datastore->finalStats(new TString(""),&finalcut);
}
