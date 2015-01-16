#ifndef _H_PGSAPI_H_
#define _H_PGSPAPI_H_

#ifdef WIN32
    #define PGS_DLLIMPORT __declspec(dllimport)
    #define PGS_DLLEXPORT __declspec(dllexport)
#else
    #define PGS_DLLIMPORT
    #define PGS_DLLEXPORT
#endif

#ifdef PGS_EXPORT
    #define PGS_API PGS_DLLEXPORT
#else
    #define PGS_API PGS_DLLIMPORT
#endif
    
#endif