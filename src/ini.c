/**
 * @file ini.c
 * @author Alan R. Rogers
 * @brief Functions for objects of class Ini, which reads parameters
 * from an initialization file  
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "misc.h"
#include "pophist.h"
#include "tokenizer.h"
#include "assign.h"
#include "ini.h"

/*
 *     Parameters         Used in
 *
 *     basepairs          ms fitld predld
 *     PopHist            fitld predld
 *     blocksize          obsld
 *     bootfilename       obsld fitld
 *     bootreps           obsld
 *     confidence         fitld
 *     doEquilibria       predld
 *     hiC               predld
 *     loC               predld
 *     methods            fitld predld
 *     nbins              obsld predld
 *     nthreads           obsld fitld
 *     samplingInterval  obsld
 *     twoNsmp            fitld predld
 *     u                  fitld predld
 *     verbose            obsld
 *     windowsize_cm      obsld
 */

/** The information read from an initialization file. */
struct Ini {
    EpochLink  *epochList;      /**< population history */
    Assignment *a;              /**< linked list of assignments */
};

/**
 * Open and read an initialization file, putting the information
 * therein into a newly-allocated object of type Ini.
 *
 * @param[in] ifname Name of input file
 *
 * @returns Newly-allocate Ini object containing info from input file.
 */
Ini        *Ini_new(const char *ifname) {
    FILE       *ifp = fopen(ifname, "r");
    int         inPopHist = 0;

    if(ifp == NULL)
        return NULL;

    Ini        *ini = malloc(sizeof(Ini));

    checkmem(ini, __FILE__, __LINE__);
    memset(ini, 0, sizeof(Ini));
    ini->a = NULL;
    ini->epochList = NULL;

    Tokenizer  *tkz = Tokenizer_new(100);
    char        buff[1000];
    int         lineno = 0, ntokens;

    while(fgets(buff, sizeof(buff), ifp) != NULL) {

        ++lineno;

        if(!strchr(buff, '\n') && !feof(ifp))
            eprintf("ERR@%s:%d: Buffer overflow. buff=\"%s\"\n",
                    __FILE__, __LINE__, buff);

        /* skip blank lines and comments */
        stripComment(buff);
        if(strempty(buff))
            continue;

        if(inPopHist) {
            Tokenizer_split(tkz, buff, " \t");  /* tokenize */
            ntokens = Tokenizer_strip(tkz, " \t\n");    /* strip extraneous */
            if(ntokens != 2)
                break;

            double      t = strtod(Tokenizer_token(tkz, 0), NULL);
            double      twoN = strtod(Tokenizer_token(tkz, 1), NULL);

            ini->epochList = EpochLink_new(ini->epochList, t, twoN);
            if(!isfinite(t))
                break;
        } else if(strchr(buff, '=')) {
            Tokenizer_split(tkz, buff, "=");    /* tokenize */
            ntokens = Tokenizer_strip(tkz, " \t\n");    /* strip extraneous */
            if(ntokens != 2)
                eprintf("ERR@:%s:%d:"
                        "Broken assignment @ line %u"
                        " of initialization file",
                        __FILE__, __LINE__, lineno);

            ini->a = Assignment_new(ini->a,
                                    Tokenizer_token(tkz, 0),
                                    Tokenizer_token(tkz, 1));

        } else {
            Tokenizer_split(tkz, buff, " \t");  /* tokenize */
            ntokens = Tokenizer_strip(tkz, " \t\n");    /* strip
                                                         * extraneous */
            if(ntokens == 0)
                continue;
            if(ntokens != 1)
                eprintf("ERR@:%s:%d:"
                        "Broken command @ line %u"
                        " of initialization file."
                        " inPopHist=%d; ntokens=%d\n",
                        __FILE__, __LINE__, lineno, inPopHist, ntokens);
            if(!strcmp("PopHist", Tokenizer_token(tkz, 0)))
                inPopHist = 1;
            else
                ini->a = Assignment_new(ini->a, Tokenizer_token(tkz, 0), "1");
        }
    }

    Tokenizer_free(tkz);
    fclose(ifp);
    return ini;
}

/** Free object of type Ini */
void Ini_free(Ini * ini) {

    if(ini == NULL)
        return;
    EpochLink_free(ini->epochList);
    Assignment_free(ini->a);
    free(ini);
}

/** Print object of type Ini to file ofp. */
void Ini_print(Ini * ini, FILE * ofp) {
    assert(ini);
    Assignment_print(ini->a, ofp);
    fprintf(ofp, "PopHist:\n");
    EpochLink_print(0, ini->epochList, ofp);
}

/**
 * Return the value of twoN in the 0th epoch.
 */
double Ini_twoN0(const Ini * ini) {
    assert(ini);
    assert(ini->epochList);
    return EpochLink_twoN(ini->epochList);
}

/**
 * Initialize *ptr to value corresponding to given key.
 *
 * @param[in] ini A structure containing the assigned values.
 *
 * @param[in] mandatory If mandatory!=0, the function aborts execution
 * if "key" is missing from Ini structure.
 *
 * @param[out] ptr Points to a double variable, which will be
 * initialized. If variable was set in initialization file, then this
 * value will be written into *ptr. Otherwise, *ptr is unaltered.
 *
 * @returns 1 if variable was initialized, 0 otherwise.
 */
int Ini_setDbl(const Ini * ini, const char *key, double *ptr, int mandatory) {
    assert(ini);
    return Assignment_setDbl(ini->a, key, ptr, mandatory);
}

/**
 * Initialize *ptr do value corresponding to given key.
 *
 * @param[in] ini A structure containing the assigned values.
 *
 * @param[in] mandatory If mandatory!=0, the function aborts execution
 * if "key" is missing from Ini structure.
 *
 * @param[out] ptr Points to a variable of type "long int", which will
 * be initialized. If variable was set in initialization file, then
 * this value will be written into *ptr. Otherwise, *ptr is unaltered.
 *
 * @returns 1 if variable was initialized, 0 otherwise.
 */
int Ini_setLong(const Ini * ini, const char *key, long *ptr, int mandatory) {
    assert(ini);
    return Assignment_setLong(ini->a, key, ptr, mandatory);
}

/**
 * Initialize *ptr do value corresponding to given key.
 *
 * @param[in] ini A structure containing the assigned values.
 *
 * @param[in] mandatory If mandatory!=0, the function aborts execution
 * if "key" is missing from Ini structure.
 *
 * @param[out] ptr Points to a variable of type "int", which will
 * be initialized. If variable was set in initialization file, then
 * this value will be written into *ptr. Otherwise, *ptr is unaltered.
 *
 * @returns 1 if variable was initialized, 0 otherwise.
 */
int Ini_setInt(const Ini * ini, const char *key, int *ptr, int mandatory) {
    assert(ini);
    return Assignment_setInt(ini->a, key, ptr, mandatory);
}

/**
 * Initialize *ptr do value corresponding to given key.
 *
 * @param[in] ini A structure containing the assigned values.
 *
 * @param[in] mandatory If mandatory!=0, the function aborts execution
 * if "key" is missing from Ini structure.
 *
 * @param[out] ptr Points to a variable of type "unsigned int", which will
 * be initialized. If variable was set in initialization file, then
 * this value will be written into *ptr. Otherwise, *ptr is unaltered.
 *
 * @returns 1 if variable was initialized, 0 otherwise.
 */
unsigned Ini_setUnsignedInt(const Ini * ini, const char *key, unsigned *ptr,
                            int mandatory) {
    assert(ini);
    return Assignment_setUnsignedInt(ini->a, key, ptr, mandatory);
}

/**
 * Initialize *ptr do value corresponding to given key.
 *
 * @param[in] ini A structure containing the assigned values.
 *
 * @param[in] mandatory If mandatory!=0, the function aborts execution
 * if "key" is missing from Ini structure.
 *
 * @param[out] ptr Points to a character array containing "size"
 * bytes. If the variable "key" is assigned within "ini", its value
 * will be copied into "ptr". Otherwise, *ptr is unaltered.
 *
 * @returns 1 if variable was initialized, 0 otherwise.
 */
int Ini_setString(const Ini * ini, const char *key, char *ptr, int size,
                  int mandatory) {
    assert(ini);
    return Assignment_setString(ini->a, key, ptr, size, mandatory);
}

/**
 * Initialize linked list of Epochs from value (if any) specified in
 * initialization file.
 *
 * @param[in] ini A structure containing the values specified in the
 * initialization file.
 *
 * @param[in] mandatory If mandatory!=0, the function aborts execution
 * if "key" is missing from Ini structure.
 *
 * @param[out] ptr Points to variable to be initialized. If variable
 * was set in initialization file, then this value will be written
 * into *ptr. Otherwise, *ptr is unaltered. If ptr is initialized, it
 * will point to a newly-allocated linked list of EpochLink
 * structures, which must be freed by the calling routine.
 *
 * @returns 1 if variable was initialized, 0 otherwise.
 */
int Ini_setEpochLink(const Ini * ini, EpochLink ** ptr, int mandatory) {
    assert(ini);
    if(ini->epochList) {
        *ptr = EpochLink_dup(ini->epochList);
        return 1;
    }
    if(mandatory)
        eprintf("ERR@%s:%d: \"%s\" not set in initialization file",
                __FILE__, __LINE__, "PopHist");
    return 0;
}

/* Local Variables: */
/* mode: c */
/* End: */
