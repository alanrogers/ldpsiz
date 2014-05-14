/**
 * @file typedefs.h
 * @author Alan R. Rogers
 * @brief Typedefs of all classes
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_TYPEDEFS_H
#define LDPSIZ_TYPEDEFS_H

#define INIFILE "ldpsiz.ini"
#define FILENAMESIZE 200
#define VERBOSE 1

typedef struct Assignment Assignment;
typedef struct Boot Boot;
typedef struct BootConf BootConf;
typedef struct Chain Chain;
typedef struct EpochLink EpochLink;
typedef struct FileIndex FileIndex;
typedef struct HillData HillData;
typedef struct Ini Ini;
typedef struct JobQueue JobQueue;
typedef struct LocRef LocRef;
typedef struct MinimizerResult MinimizerResult;
typedef struct Model Model;
typedef struct ModelList ModelList;
typedef struct ObjLink ObjLink;
typedef struct ODE ODE;
typedef struct PopHist PopHist;
typedef struct SNP SNP;
typedef struct SNPLoc SNPLoc;
typedef struct SNPstore SNPstore;
typedef struct SeqDatPar SeqDatPar;
typedef struct SimplexInfo SimplexInfo;
typedef struct SimAnn SimAnn;
typedef struct Tabulation Tabulation;
typedef struct ThreadBounds ThreadBounds;
typedef struct Threadpool Threadpool;
typedef struct ThreadQueue ThreadQueue;
typedef struct Tokenizer Tokenizer;
typedef struct Window Window;

struct dydt_params {
    double      c, u, twoN, duration;
};

#endif
