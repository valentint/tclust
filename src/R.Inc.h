
#ifdef _MSC_VER
	#define EXPORT extern "C" __declspec(dllexport)
#else
	#define EXPORT extern "C"
#endif //#ifdef _MSC_VER

