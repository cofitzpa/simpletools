#include "cropcutspace.h"
cropcutspace::cropcutspace(Double_t _res, Double_t _min, Double_t _max){
	steps = (UInt_t)((_max - _min)/_res) +1;
	res = _res;
	maxFoMStep = 0;
	Vals = new Double_t [steps];
	S = new Double_t [steps];
	B = new Double_t [steps];
	Seff = new Double_t [steps];
	Brej = new Double_t [steps];
	FoM = new Double_t [steps];
	d_Vals = new Double_t [steps];
	d_S = new Double_t [steps];
	d_B = new Double_t [steps];
	d_Seff = new Double_t [steps];
	d_Brej = new Double_t [steps];
	d_FoM = new Double_t [steps];
	setVals(_min);
}


cropcutspace::cropcutspace(UInt_t nSteps, Double_t _min, Double_t _max){
	steps = nSteps+1;
	res = (_max - _min)/((Double_t)nSteps);
	maxFoMStep = 0;
	Vals = new Double_t [steps];
	S = new Double_t [steps];
	B = new Double_t [steps];
	Seff = new Double_t [steps];
	Brej = new Double_t [steps];
	FoM = new Double_t [steps];
	d_Vals = new Double_t [steps];
	d_S = new Double_t [steps];
	d_B = new Double_t [steps];
	d_Seff = new Double_t [steps];
	d_Brej = new Double_t [steps];
	d_FoM = new Double_t [steps];
	setVals(_min);
}

cropcutspace::~cropcutspace(){}

void cropcutspace::addPoint(UInt_t * step, Double_t * _s, Double_t * _d_s, Double_t * _b, Double_t * _d_b, Double_t * _seff, Double_t * _d_seff,Double_t * _brej,Double_t * _d_brej,Double_t *_fom,Double_t * _d_fom){
	S[*step]=*_s;
	d_S[*step]=*_d_s;
	B[*step]=*_b;
	d_B[*step]=*_d_b;
	Seff[*step]=*_seff;
	d_Seff[*step]=*_d_seff;
	Brej[*step]=*_brej;
	d_Brej[*step]=*_d_brej;
	FoM[*step]=*_fom;
	d_FoM[*step]=*_d_fom;
}

UInt_t cropcutspace::setMaxFoMStep(bool gtCut){
	if(!gtCut){ 
		//LTCUT: x<number means start low and go high to make more events available so loop is positive
		//If the current step has a FoM equal to the previous maximum, make this the maxFoMstep as looser is safer in the ensemble
		Double_t max = FoM[0];
		maxFoMStep = 0;
		for(UInt_t i = 1; i<steps; i++){
			if(FoM[i] >= max){
				max = FoM[i]; 
				maxFoMStep = i;
			}
		}

	}else{
		//GTCUT: x>number means start high and go low to make more events available so loop is negative
		//If the current step has a FoM equal to the previous maximum, make this the maxFomStep to loosen
		Double_t max = FoM[steps-1];
		maxFoMStep = steps-1;
		for(Int_t j = steps-2; j>=0; j--){
			if(FoM[j] >= max){
				max = FoM[j]; 
				maxFoMStep = j;
			}
		}

	}
	return maxFoMStep; 
}

void cropcutspace::setVals(Double_t min){
	for(UInt_t n =0; n<steps; n++){
		Vals[n] = min + (((Double_t)n)*res);
		d_Vals[n] = res/2.0;
		S[n] = 0.0;
		d_S[n] = 0.0;
		B[n] = 0.0;
		d_B[n] = 0.0;
		Seff[n] = 0.0;
		d_Seff[n] = 0.0;
		Brej[n] = 0.0;
		d_Brej[n] = 0.0;
		FoM[n] = 0.0;
		d_FoM[n] = 0.0;
	}

}

void cropcutspace::plotFoM(TString name)const{
	TGraphErrors * FoMGraph = new TGraphErrors(steps,Vals,FoM,d_Vals,d_FoM);
	FoMGraph->SetLineWidth(2);
	FoMGraph->Draw("AL1");
	FoMGraph->SetTitle("Figure of Merit for cut "+name);
	FoMGraph->GetXaxis()->SetTitle(name);
	FoMGraph->GetYaxis()->SetTitle("FoM");

	TArrow* FoMarrow1 = new TArrow(Vals[maxFoMStep],FoMGraph->GetYaxis()->GetXmin(),Vals[maxFoMStep],FoM[maxFoMStep],0.01,"<-");
	FoMarrow1->SetLineWidth(2) ;
	FoMarrow1->SetLineColor(kBlue);
	FoMarrow1->Draw();

	TArrow* FoMarrow2 = new TArrow(FoMGraph->GetXaxis()->GetXmin(),FoM[maxFoMStep],Vals[maxFoMStep],FoM[maxFoMStep],0.01,"<-");
	FoMarrow2->SetLineWidth(2) ;
	FoMarrow2->SetLineColor(kBlue);
	FoMarrow2->Draw();
}

void cropcutspace::plotSeff(TString name)const{
	TGraphErrors * SeffGraph = new TGraphErrors(steps,Vals,Seff,d_Vals,d_Seff);
	SeffGraph->SetLineWidth(2);
	SeffGraph->Draw("AL1");
	SeffGraph->SetTitle("Signal Efficiency for cut "+name);
	SeffGraph->GetXaxis()->SetTitle(name);
	SeffGraph->GetYaxis()->SetTitle("#varepsilon_{S}");

	TArrow* Seffarrow1 = new TArrow(Vals[maxFoMStep],SeffGraph->GetYaxis()->GetXmin(),Vals[maxFoMStep],Seff[maxFoMStep],0.01,"<-");
	Seffarrow1->SetLineWidth(2) ; 
	Seffarrow1->SetLineColor(kBlue);
	Seffarrow1->Draw();

	TArrow* Seffarrow2 = new TArrow(SeffGraph->GetXaxis()->GetXmin(),Seff[maxFoMStep],Vals[maxFoMStep],Seff[maxFoMStep],0.01,"<-");
	Seffarrow2->SetLineWidth(2) ; 
	Seffarrow2->SetLineColor(kBlue);
	Seffarrow2->Draw();
}


void cropcutspace::plotBrej(TString name)const{
	TGraphErrors * BrejGraph = new TGraphErrors(steps,Vals,Brej,d_Vals,d_Brej);
	BrejGraph->SetLineWidth(2);
	BrejGraph->Draw("AL1");
	BrejGraph->SetTitle("Background rejection for cut "+name);
	BrejGraph->GetXaxis()->SetTitle(name);
	BrejGraph->GetYaxis()->SetTitle("1-#varepsilon_{B}");

	TArrow* Brejarrow1 = new TArrow(Vals[maxFoMStep],BrejGraph->GetYaxis()->GetXmin(),Vals[maxFoMStep],Brej[maxFoMStep],0.01,"<-");
	Brejarrow1->SetLineWidth(2) ; 
	Brejarrow1->SetLineColor(kBlue);
	Brejarrow1->Draw();

	TArrow* Brejarrow2 = new TArrow(BrejGraph->GetXaxis()->GetXmin(),Brej[maxFoMStep],Vals[maxFoMStep],Brej[maxFoMStep],0.01,"<-");
	Brejarrow2->SetLineWidth(2) ; 
	Brejarrow2->SetLineColor(kBlue);
	Brejarrow2->Draw();
}



void cropcutspace::plotRoC(TString name)const{
	TGraphErrors * RoCGraph = new TGraphErrors(steps,Seff,Brej,d_Seff,d_Brej);
	RoCGraph->SetLineWidth(2);
	RoCGraph->Draw("AL1");
	RoCGraph->SetTitle("Signal efficiency vs. Background rejection for cut "+name);
	RoCGraph->GetXaxis()->SetTitle("#varepsilon_{S}");
	RoCGraph->GetYaxis()->SetTitle("1-#varepsilon_{B}");

	TArrow* RoCarrow1 = new TArrow(Seff[maxFoMStep],RoCGraph->GetYaxis()->GetXmin(),Seff[maxFoMStep],Brej[maxFoMStep],0.01,"<-");
	RoCarrow1->SetLineWidth(2) ; 
	RoCarrow1->SetLineColor(kBlue);
	RoCarrow1->Draw();

	TArrow* RoCarrow2 = new TArrow(RoCGraph->GetXaxis()->GetXmin(),Brej[maxFoMStep],Seff[maxFoMStep],Brej[maxFoMStep],0.01,"<-");
	RoCarrow2->SetLineWidth(2) ; 
	RoCarrow2->SetLineColor(kBlue);
	RoCarrow2->Draw();
}

void cropcutspace::plotAll(TString name, TCanvas* c)const{
	c->Clear();
	c->Divide(2,2);
	c->cd(1);
	plotFoM(name);
	c->Update();
	c->cd(2);
	plotRoC(name);
	c->Update();
	c->cd(3);
	plotSeff(name);
	c->Update();
	c->cd(4);
	plotBrej(name);
	c->Update();
}
