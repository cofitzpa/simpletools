#ifndef cropoptimisationengine_HH
#define cropoptimisationengine_HH


#include "cropdatastore.h"
#include "cropcutensemble.h"
#include "cropFoMs.h"
#include "TRandom3.h"
#include "TFunctor.h"

using std::cout;
using std::endl;
class cropoptimisationengine {
	public:
		cropoptimisationengine(cropdatastore *,cropcutensemble *,TCanvas *, TRandom3 *, void (*)(cropdatastore *, TString *, Double_t *, Double_t *), TFunctor *);
		void initialise();
		bool optimise(UInt_t);
		void writeEvolutionPlot();
		void finalStats();

	private:
		cropdatastore * datastore;
		cropcutensemble * cutensemble;
		void (*fom)(cropdatastore *, TString *, Double_t *, Double_t *);
		TCanvas *c;
		TRandom3 *rng;
		vector<TString> cutAtStep;
		vector<Double_t> fomAtStep;
		vector<Double_t> d_fomAtStep;
		vector<Double_t> Step;
		TFunctor *OrderMethod;
};
#endif
