#ifndef _H_MANIFOLDSAPI_H_
#define _H_MANIFOLDSAPI_H_

#ifdef WIN32
    #define MANIFOLDS_DLLIMPORT __declspec(dllimport)
    #define MANIFOLDS_DLLEXPORT __declspec(dllexport)
#else
    #define MANIFOLDS_DLLIMPORT
    #define MANIFOLDS_DLLEXPORT
#endif

#ifdef MANIFOLDS_EXPORT
    #define MANIFOLDS_API MANIFOLDS_DLLEXPORT
#else
    #define MANIFOLDS_API MANIFOLDS_DLLIMPORT
#endif

#endif
