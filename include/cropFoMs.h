#ifndef CROP_foms_HH
#define CROP_foms_HH
#include "cropdatastore.h"
#include "cropmiscfunctions.h"

//Punzi FoM at the request of Steve-O and Young John. Significance set to 3.0 sigma, change and recompile if you want different. 
inline void cropFoM_Punzi(cropdatastore* data =0, TString *cut =0, Double_t *FoM=0, Double_t *d_FoM=0){
	Double_t S;
	Double_t d_S;
	data->getWeightedSignalEntries(cut, &S, &d_S);
	if(S==0.0){
		*FoM=0;
		*d_FoM=0;
	}else{
		Double_t B;
		Double_t d_B;
		data->getWeightedBackgroundEntries(cut, &B, &d_B);
		Double_t sqrtb =0.0;
		Double_t d_sqrtb =0.0;
		if(B>0){
		quadPow(B, d_B, 0.5, &sqrtb, &d_sqrtb);
		}
		Double_t punzidenom = 3.0/2.0; //Change 3.0 to desired signif... 
		punzidenom += sqrtb;
		quadEff(S,d_S,punzidenom,d_sqrtb,FoM,d_FoM);
//	cout << *FoM << "	" << *d_FoM << endl;
	}
}

inline void cropFoM_SoRSB(cropdatastore* data =0, TString *cut =0, Double_t *SoRSB=0, Double_t *d_SoRSB=0){
	Double_t S;
	Double_t d_S;
	data->getWeightedSignalEntries(cut, &S, &d_S);
	if(S==0.0){
		*SoRSB=0;
		*d_SoRSB=0;
	}else{
		Double_t B;
		Double_t d_B;
		data->getWeightedBackgroundEntries(cut, &B, &d_B);
		Double_t splusb;
		Double_t d_splusb;
		quadSum(S,d_S,B,d_B,&splusb,&d_splusb);
		Double_t sqrtsplusb;
		Double_t d_sqrtsplusb;
		quadPow(splusb, d_splusb, 0.5, &sqrtsplusb, &d_sqrtsplusb);
		quadEff(S,d_S,sqrtsplusb,d_sqrtsplusb,SoRSB,d_SoRSB);
	}
}

inline void cropFoM_betas(cropdatastore* data=0, TString *cut=0, Double_t *FoM=0, Double_t *d_FoM=0){
//TString fomhstring = "B_s_TAU";
//TString fomfstring = "((1.0-2.0*B_s_TAGOMEGA)**2)*exp(-(20.0*B_s_TAUERR)**2)";
TString fomhstring = "time";
//TString fomfstring = "(((1.0-2.0*B_s0_TAGOMEGA)**2)*exp(-((20000.0*B_s0_TAUERR)**2)))";

//TString fomfstring = "(((1.0-2.0*mistag)**2)*exp(-(20.0*(time*(B_s0_LOKI_DTF_CTAUS/B_s0_LOKI_DTF_CTAU)))**2.0))";
TString fomfstring = "((1.0-2.0*mistag)**2)";
UInt_t fomhbins = 6;
Double_t fomhmin = -3.0;
Double_t fomhmax = 12.0;
TH1D* sigh;
TH1D* bkgh;
	//      TCanvas *canv = new TCanvas("canv","canv",800,600);
	//      //      canv->Divide(2,3);
	//      //      canv->cd(1);
	sigh = data->getSignalHisto(*cut,fomhstring,fomhbins,fomhmin,fomhmax);
	sigh->SetName("sigh");
	//      sigh->Print("all");
	//      //      Double_t sig = sigh->Integral();
	//      //      Double_t dsig = sqrt(sigh->GetSumOfWeights());
	//      //      sigh->Draw();
	//      //      canv->Update();
	//      //      canv->cd(2);
	//
	bkgh = data->getBackgroundHisto(*cut,fomhstring,fomhbins,fomhmin,fomhmax);
	bkgh->SetName("bkgh");
	//      bkgh->Print("all");
	//      //      Double_t bkg = bkgh->Integral();
	//      //      Double_t dbkg = sqrt(bkgh->GetSumOfWeights());
	//      //      bkgh->Draw();
	//      //      canv->Update();
	//      //      canv->cd(3);
	//
	TString errForm = "(";
	errForm+= *cut;
	errForm += ")*";
	errForm += fomfstring;
	//      cout << errForm << endl;
	TH1D * localSignificance = data->getSignalHisto(errForm,fomhstring,fomhbins,fomhmin,fomhmax);
	localSignificance->SetName("localSignificance");
	//      localSignificance->Print("all");
	//      //      Double_t localSig = localSignificance->Integral();
	//      //      Double_t dlocalSig = sqrt(localSignificance->GetSumOfWeights());
	//      //      localSignificance->Draw();
	//      //      canv->Update();
	//      //      canv->cd(4);
	//
	localSignificance->Multiply(localSignificance,sigh); //S^2
	//      localSignificance->Print("all");
	//      //      Double_t sigsq = localSignificance->Integral();
	//      //      Double_t dsigsq =  sqrt(localSignificance->GetSumOfWeights());
	//      //      localSignificance->Draw();
	//      //      canv->Update();
	//      //      canv->cd(5);
	sigh->Add(bkgh); //S+B
	//      sigh->Print("all");
	//      //      Double_t splusb = sigh->Integral();
	//      //      Double_t dsplusb = sqrt(sigh->GetSumOfWeights());
	//      //      sigh->Draw();
	//      //      canv->Update();
	//      //      canv->cd(6);
	localSignificance->Divide(sigh); //S^{2}/S+B
	//      localSignificance->Print("all");
	//      //      Double_t fom = localSignificance->Integral();
	//      //      Double_t dfom = sqrt(localSignificance->GetSumOfWeights());
	//      //      localSignificance->Draw();
	//      //      canv->Update();
	//
	//
	//      //      *FoM = localSignificance->Integral();
	*FoM = localSignificance->IntegralAndError(0,(Int_t)fomhbins,*d_FoM);
	//      *d_FoM = sqrt(localSignificance->GetSumOfWeights());
	sigh->Delete();
	bkgh->Delete();
	localSignificance->Delete();
	//      cout << "sigh =" << sig << "+/-" << dsig  << " bkgh =" << bkg << "+/-" << dbkg << " s+b =" << splusb << "+/-" << dsplusb << " localsig=" << localSig << "+/-" << dlocalSig << " s*localsig=" << sigsq << "+/-" << dsigsq << " fom="  << fom << "+/-" << dfom << endl;
}


#endif
