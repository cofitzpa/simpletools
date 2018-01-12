#ifndef TFunctor_HH
#define TFunctor_HH

class TFunctor{
	public:

		// two possible functions to call member function. virtual cause derived
		// classes will use a pointer to an object and a pointer to a member function
		// to make the function call
		virtual ~TFunctor(){}
		virtual bool operator()()=0;  // call using operator
		virtual bool Call()=0;        // call using function



};

// derived template class
template <class TClass> class TSpecificFunctor : public TFunctor
{
	private:
		bool (TClass::*fpt)();   // pointer to member function
		TClass* pt2Object;                  // pointer to object

	public:
		// constructor - takes pointer to an object and pointer to a member and stores
		// them in two private variables
		~TSpecificFunctor(){}
		TSpecificFunctor(TClass* _pt2Object, bool(TClass::*_fpt)())
		{ pt2Object = _pt2Object;  fpt=_fpt; };

		// override operator "()"
		virtual bool operator()()
		{ return (*pt2Object.*fpt)();};              // execute member function
		// override function "Call"
		virtual bool Call()
		{ return (*pt2Object.*fpt)();};             // execute member function
};
#endif
