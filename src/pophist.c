/**
 * @file pophist.c
 * @author Alan R. Rogers
 * @brief Functions for objects of type PopHist, which represent the history
 *        of population size.
 *
 * The entire history of population size, represented as a sequence of Epochs.
 * Population parameters may change at Epoch boundaries but are
 * constant within epochs.
 *
 * In this version, all parameters are in a single array, and twoN
 * parameters are coded as 1/twoN.
 *
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "misc.h"
#include "pophist.h"
#include <assert.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** The number of parameters in each Epoch */
#define PAR_PER_EPOCH 2

/** Names of the parameters within each lpoch */
const static char *epochParamName[] = { "twoN", "t" };

/**
 * PopHist represents the history of a single population.
 * It is a structure allocated on a single block of memory, using the
 * "struct hack" of C programming. This makes it possible to avoid
 * pointers within PopHist, so that these objects can be copied using
 * memcpy. (Otherwise, you end up with pointers in one object pointing
 * to memory in another.)
 */
struct PopHist {
    unsigned    nepoch;    /**< number of epochs */
    size_t      size;      /**< size of memory block allocated */
    double      p[PAR_PER_EPOCH];   /**< array of parameter values */
};

/**
 * Population history represented as a linked list of epochs.
 */
struct EpochLink {
    struct EpochLink *next;     /* pointer to next (earlier) epoch */
    double      t, twoN;
};

static inline int isPopSize(int ndx);

/// Return true if index refers to a population size parameter.
static inline int isPopSize(int ndx) {
    if(ndx % PAR_PER_EPOCH)
        return false;
    return true;
}

void PopHist_sanityCheck(PopHist * ph, const char *file, int lineno) {
#ifndef NDEBUG
    REQUIRE(ph != NULL, file, lineno);
    REQUIRE(ph->nepoch > 0, file, lineno);

    void       *end = (void *) &ph->p[2 * ph->nepoch];
    int         i;

    REQUIRE(ph->nepoch >= 1, file, lineno);
    REQUIRE(ph->size == ((size_t) end) - (size_t) ph, file, lineno);
    REQUIRE(ph->size == sizeof(PopHist)
            + PAR_PER_EPOCH * (ph->nepoch - 1) * sizeof(double), file,
            lineno);
    for(i = 0; i < ph->nepoch; ++i) {
        REQUIRE(0 <= PopHist_twoN(ph, i), file, lineno);
        REQUIRE(0 <= PopHist_duration(ph, i), file, lineno);
    }
#endif
}

/**
 * Print a linked list.
 *
 * @param ndx A state variable, which should be zero in the top-level
 * call.
 * @param link Pointer to the head of the linked list.
 * @param ofp Pointer to output file.
 */
void EpochLink_print(int ndx, const EpochLink * link, FILE * ofp) {
    if(link == NULL)
        return;
    printf("%4d:t=%lg 2N=%lg\n", ndx, link->t, link->twoN);
    EpochLink_print(++ndx, link->next, ofp);
    return;
}

/**
 * Insert a new link into the list.
 *
 * Add item to linked list, allocating as necessary.
 * @returns pointer to head of list.
 */
EpochLink  *EpochLink_new(EpochLink * head, double t, double twoN) {
    EpochLink  *tail, *curr;

    curr = (EpochLink *) malloc(sizeof(EpochLink));
    checkmem(curr, __FILE__, __LINE__);

    curr->next = NULL;
    curr->t = t;
    curr->twoN = twoN;

    if(head == NULL)
        return curr;
    for(tail = head; tail->next != NULL; tail = tail->next) ;
    tail->next = curr;
    return head;
}

/** Return a newly allocated copy of a linked list of epochs. */
EpochLink  *EpochLink_dup(EpochLink * head) {
    if(head == NULL)
        return NULL;

    EpochLink  *new = malloc(sizeof(EpochLink));

    checkmem(new, __FILE__, __LINE__);

    new->t = head->t;
    new->twoN = head->twoN;
    new->next = EpochLink_dup(head->next);

    return new;
}

/**
 * free linked list
 */
void EpochLink_free(EpochLink * link) {
    if(link == NULL)
        return;
    EpochLink_free(link->next);
    free(link);
    return;
}

/**
 * Count the number of links in a linked list of EpochLink objects.
 *
 * @param head pointer to beginning of linked list.
 * @returns Number of links in the chain beginning with "head".
 */
int EpochLink_nlinks(const EpochLink * head) {
    if(head == NULL)
        return 0;

    return 1 + EpochLink_nlinks(head->next);
}

/**
 * Return a pointer to the next link in a chain of EpochLinks.
 */
const EpochLink *EpochLink_next(const EpochLink * link) {
    myassert(link);
    return link->next;
}

/** Return time value from given link. */
double EpochLink_duration(const EpochLink * link) {
    myassert(link);
    return link->t;
}

/** Return twoN value from given link. */
double EpochLink_twoN(const EpochLink * link) {
    myassert(link);
    return link->twoN;
}

#ifndef NDEBUG
void EpochLink_test(void) {
    EpochLink  *el = NULL;

    el = EpochLink_new(el, 100.0, 100000.0);
    el = EpochLink_new(el, 10.0, 100.0);
    el = EpochLink_new(el, 50.0, 10000.0);  /* 50 should become Inf */
    assert(el);
    assert(el->next);
    assert(el->next->next);
    assert(3 == EpochLink_nlinks(el));
    assert(2 == EpochLink_nlinks(EpochLink_next(el)));
    assert(100.0 == EpochLink_duration(el));
    assert(100000.0 == EpochLink_twoN(el));

    EpochLink  *el2 = EpochLink_dup(el);

    assert(el->t == el2->t);
    assert(el->twoN == el2->twoN);
    assert(el->next->t == el2->next->t);
    assert(el->next->twoN == el2->next->twoN);
    assert(el->next->next->t == el2->next->next->t);
    assert(el->next->next->twoN == el2->next->next->twoN);

    EpochLink_free(el);
    EpochLink_free(el2);
}
#endif

/**
 * Set up a PopHist object, which has previously been allocated
 * as a contiguous memory block.
 */
void PopHist_init(PopHist * ph, unsigned nepoch, size_t size) {
    ph->nepoch = nepoch;
    ph->size = size;
    memset(ph->p, 0, nepoch * PAR_PER_EPOCH * sizeof(ph->p[0]));
    ph->p[1 + (ph->nepoch - 1) * PAR_PER_EPOCH] = strtod("Inf", 0);
    PopHist_sanityCheck(ph, __FILE__, __LINE__);
}

/**
 * Allocate an empty PopHist.
 *
 * The parameters of each epoch are initialized with zeroes, except for
 * the duration of the final (earliest) epoch, which is infinite.
 *
 * @param[in] nepoch Number of epochs in the new PopHist.
 *
 * @returns a pointer to a newly allocated PopHist object with default values.
 */
PopHist    *PopHist_newEmpty(unsigned nepoch) {
    size_t      size = PopHist_calcSize(nepoch);
    PopHist    *ph = malloc(size);

    checkmem(ph, __FILE__, __LINE__);
    PopHist_init(ph, nepoch, size);
    return ph;
}

/** Return the size of the memory block on which `ph` is allocated. */
size_t PopHist_size(const PopHist * ph) {
    return ph->size;
}

/**
 * Return the number of bytes occupied by a PopHist with "nepoch"
 * epochs.
 */
size_t PopHist_calcSize(unsigned nepoch) {
    if(nepoch == 0) {
        dostacktrace(__FILE__, __LINE__, stderr);
        eprintf("ERR in%s@%s:%d: PopHist must have at least one epoch\n",
                __func__, __FILE__, __LINE__);
    }
    return sizeof(PopHist) + PAR_PER_EPOCH * (nepoch - 1) * sizeof(double);
}

/**
 * Return duration of Epoch i.
 */
double PopHist_duration(const PopHist * ph, int i) {
    return ph->p[1 + i * PAR_PER_EPOCH];
}

/**
 * Return time from leaves to the most recent end of Epoch i. Thus,
 * PopHist_age(ph, 0) is 0 and PopHist_age(ph, PopHist_nepoch(ph))
 * gives the sum of all epoch durations except the final infinite one.
 */
double PopHist_age(const PopHist *ph, int i) {
	double t=0.0;
	int j;
	for(j=0; j<i; ++j)
		t += PopHist_duration(ph, j);
	return t;
}

/** Return the haploid population size, 2N, during epoch i. */
double PopHist_twoN(const PopHist * ph, int i) {
    return 1.0 / ph->p[i * PAR_PER_EPOCH];
}

/**
 * Find epoch that covers a given age.
 */
int PopHist_findEpoch(const PopHist *ph, double age) {
	double thusfar=0.0;
	int epoch=0;
    double duration = PopHist_duration(ph, epoch);

    while( thusfar + duration < age) {
        thusfar += duration;
        ++epoch;
        duration = PopHist_duration(ph, epoch);
    }
    
	return epoch;
}

/** Return 1/2N for the relevant epoch */
double PopHist_twoNinv(const PopHist * ph, int i) {
    return ph->p[i * PAR_PER_EPOCH];
}

/**
 * Set duration of the i'th epoch of PopHist structure ph.
 *
 * @param[out] ph PopHist object to be modified.
 * @param[in] i index of epoch to be modified. Must
 *            be less than ph->nepoch-1.
 * @param[in] duration new value of duration of i'th epoch.
 */
void PopHist_setDuration(PopHist * ph, int i, double duration) {
    if(i + 1 == ph->nepoch)
        eprintf("ERR@%s:%s:%d: Can't set duration of initial epoch\n",
                __func__, __FILE__, __LINE__);
    ph->p[1 + i * PAR_PER_EPOCH] = duration;
}

/**
 * Set population size during the i'th epoch of PopHist structure ph.
 *
 * @param[out] ph PopHist object to be modified.
 * @param[in] i index of epoch to be modified. Must
 *            be less than ph->nepoch-1.
 * @param[in] t new value of duration of i'th epoch.
 */
void PopHist_setTwoN(PopHist * ph, int i, double twoN) {
    ph->p[i * PAR_PER_EPOCH] = 1.0 / twoN;
}

/**
 * Set 1/2N during the i'th epoch of PopHist structure ph.
 *
 * @param[out] ph PopHist object to be modified.
 * @param[in] i index of epoch to be modified. Must
 *            be less than ph->nepoch-1.
 * @param[in] t new value of duration of i'th epoch.
 */
void PopHist_setTwoNinv(PopHist * ph, int i, double twoNinv) {
    ph->p[i * PAR_PER_EPOCH] = twoNinv;
}

/**
 * Allocate a new PopHist and initialize from a chain of EpochLink
 * objects.
 *
 * @param[in] head Beginning of chain of EpochLink objects.
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
PopHist    *PopHist_fromEpochLink(const EpochLink * head) {
    int         i, nepoch;
    const EpochLink *link;
    PopHist    *ph;

    nepoch = EpochLink_nlinks(head);

    if(nepoch == 0)
        return NULL;

    ph = PopHist_newEmpty(nepoch);

    for(i = 0, link = head; link != NULL; ++i, link = EpochLink_next(link)) {
        myassert(i < nepoch);

        if(link->next != NULL)
            PopHist_setDuration(ph, i, EpochLink_duration(link));
        PopHist_setTwoN(ph, i, EpochLink_twoN(link));
    }

    myassert(PopHist_duration(ph, ph->nepoch - 1) == HUGE_VAL);
    PopHist_sanityCheck(ph, __FILE__, __LINE__);

    return ph;
}

#pragma GCC diagnostic pop

/** return the number of epochs */
unsigned PopHist_nepoch(const PopHist * ph) {
    return ph->nepoch;
}

/**
 * Print PopHist, preceding each line with a comment string.
 *
 * @param ph PopHist to print.
 * @param comstr Comment string to prepend to each line of output.
 * @param outfile Output file.
 */
void PopHist_print_comment(const PopHist * ph, const char *comstr,
                           FILE * outfile) {
    int         i;

    myassert(ph != NULL);
    myassert(PopHist_duration(ph, ph->nepoch - 1) == HUGE_VAL);

    fprintf(outfile, "%s%5s %8s %10s\n", comstr, "epoch", "duration", "1/2N");
    for(i = 0; i < ph->nepoch; ++i) {
        fprintf(outfile, "%s%5d %8.3lg %10.3lg\n",
                comstr, i, PopHist_duration(ph, i), PopHist_twoNinv(ph, i));
    }
    return;
}

/**
 * Print a PopHist object.
 *
 * @param ph PopHist to print.
 * @param outfile Output file.
 */
void PopHist_print(const PopHist * ph, FILE * outfile) {
    PopHist_print_comment(ph, "", outfile);
    return;
}

/** Destroy a PopHist object. */
void PopHist_free(PopHist * ph) {
    free(ph);
}

/**
 * This is just like PopHist_nParams, except that you don't need to
 * allocate a PopHist object before calling it. All you need is the
 * number of epochs.
 */
int PopHist_calc_nParams(int nepoch) {
    return (PAR_PER_EPOCH * nepoch) - 1;
}

/**
 * Return the number of parameters in a PopHist.
 *
 * There are two adjustable parameters per epoch: N and t, except that
 * t is fixed at infinity in the earliest epoch.  The number of
 * adjustable parameters is therefore 2*nepochs - 1.
 */
unsigned PopHist_nParams(const PopHist * ph) {
    return PopHist_calc_nParams(ph->nepoch);
}

/**
 * Return the name of a given PopHist parameter.
 *
 * @param[in] ph A pointer to a PopHist. It is used only to determine
 * the number of epochs.
 * @param[out] buff A character buffer into which the parameter's name will
 * be written. If bufflen is too short, the name will be truncated.
 * @param[in] bufflen The length of the buffer.
 * @param[in] ndx The index of the parameter, a positive integer less than
 * the number of PopHist parameters, as given by PopHist_nParams.
 * @returns 0 on success, 1 on failure.
 */
int PopHist_paramName(const PopHist * ph, char *buff, int bufflen, int ndx) {
    div_t       quotrem;
    int         rval;

    myassert(ndx >= 0);
    myassert(ndx < PopHist_nParams(ph));

    // quotrem.quot is the index of the desired epoch.
    // quotrem.rem = the index of the parameter within the epoch.
    quotrem = div(ndx, PAR_PER_EPOCH);

    rval = snprintf(buff, (unsigned) bufflen,
                    "%s%d", epochParamName[quotrem.rem], quotrem.quot);
    if(rval + 1 == bufflen) {
        fprintf(stderr, "WARNING@:%s:%d: snprintf filled buffer (%d bytes)"
                " buff=%s\n", __FILE__, __LINE__, rval, buff);
        return 1;
    }
    return 0;
}

/**
 * Return the value of a given PopHist parameter.
 *
 * @param[in] ph A pointer to a PopHist.
 * @param[in] ndx The index of the parameter, a positive integer less than
 * the number of PopHist parameters, as given by PopHist_nParams.
 * @returns Value of specified parameter within ph.
 */
double PopHist_paramValue(const PopHist * ph, int ndx) {

    myassert(ndx >= 0);
    myassert(ndx < PopHist_nParams(ph));

    if(isPopSize(ndx))
        return 1.0/ph->p[ndx];
    return ph->p[ndx];
}

/**
 * Put PopHist parameters into a gsl_vector.
 *
 * @param[in] PopHist object.
 * @param[out] v gsl_vector into which paramters will be written.
 */
void PopHist_to_vector(gsl_vector * v, const PopHist * ph) {
    int         i;

    for(i = 0; i < PopHist_nParams(ph); ++i)
        gsl_vector_set(v, i, ph->p[i]);

    return;
}

/**
 * Put PopHist parameters into an ordinary C array.
 *
 * @param[in] PopHist object.
 * @param[out] v  An array into which paramters will be written.
 * @param[in] dim The dimension of array v
 */
void PopHist_to_C_array(int dim, double v[dim], const PopHist * ph) {

    assert(dim == PopHist_nParams(ph));
    memcpy(v, ph->p, dim * sizeof v[0]);

    return;
}

/**
 * Assign the same value to all positions in array corresponding to
 * twoNinv values in the parameter vector of a PopHist object.
 */
void PopHist_setAllTwoNinv(double *x, int dim, double value) {
    int         i;

    for(i = 0; i < dim; i += 2)
        x[i] = value;
}

/**
 * Assign the same value to all positions in array corresponding to
 * duration values in the parameter vector of a PopHist object.
 */
void PopHist_setAllDuration(double *x, int dim, double value) {
    int         i;

    for(i = 1; i < dim; i += 2)
        x[i] = value;
}

/**
 * Set scale vector used by simplex algorithm. Each entry
 * receives a value that represents the expected magnitude of the
 * corresponding parameter in PopHist.
 *
 * @param[in] PopHist object.
 * @param[out] scale gsl_vector into which paramters will be written.
 */
void PopHist_setSimplexScale(gsl_vector * scale, const PopHist * ph) {
    double      twoN, duration;
    int         nepoch = ph->nepoch;
    unsigned    i;

    /*
     * This version sets scale of each argument based on the size of
     * the corresponding value in PopHist argument, unless that value
     * is very small. If it is small, scale is set using the fixed
     * constants defined here.
     */
    static const double minTwoNscale = 100.0;
    static const double minDurationScale = 10.0;
    static const double inflation = 0.9;

    for(i = 0; i < nepoch - 1; ++i) {
        twoN = inflation * gsl_vector_get(scale, 2 * i);
        twoN = fmax(twoN, minTwoNscale);
        gsl_vector_set(scale, 2 * i, twoN);

        duration = inflation * gsl_vector_get(scale, 2 * i + 1);
        duration = fmax(duration, minDurationScale);
        gsl_vector_set(scale, 2 * i + 1, duration);
    }
    myassert(i == nepoch - 1);
    gsl_vector_set(scale, 2 * i, twoN);

    return;
}

/**
 * Adjust the parameters in PopHist structure ph to reflect the values
 * in gsl vector p.
 */
void vector_to_PopHist(PopHist * ph, const gsl_vector * v) {
    unsigned    i;
    int         nepoch = ph->nepoch;
    double      t, twoNinv;

    myassert(v->size == PopHist_nParams(ph));
    for(i = 0; i < nepoch - 1; ++i) {
        twoNinv = gsl_vector_get(v, 2 * i);
        t = gsl_vector_get(v, 2 * i + 1);
        PopHist_setDuration(ph, i, t);
        PopHist_setTwoNinv(ph, i, twoNinv);
    }
    PopHist_setTwoNinv(ph, nepoch - 1, gsl_vector_get(v, 2 * i));

    return;
}

/**
 * Adjust the parameters in PopHist structure ph to reflect the values
 * in C array p, of dimension "dim".
 */
void C_array_to_PopHist(PopHist * ph, int dim, const double *v) {

    if(dim != PopHist_nParams(ph)) {
        printf("%s:%s:%d: dim=%d nParams=%d\n",
               __FILE__, __func__, __LINE__, dim, PopHist_nParams(ph));
    }
    assert(dim == PopHist_nParams(ph));
    memcpy(ph->p, v, dim * sizeof v[0]);

    return;
}

/** Return a pointer to a newly allocated copy of PopHist structure ph */
PopHist    *PopHist_dup(const PopHist * ph) {
    if(ph == NULL)
        return NULL;

    PopHist    *new = malloc(ph->size);

    checkmem(new, __FILE__, __LINE__);

    memcpy(new, ph, ph->size);
    return new;
}

/**
 * Copy PopHist src into dest.
 *
 * The two must have equal values of nepoch. Otherwise, the function
 * aborts.
 * @param src The source PopHist.
 * @param dest The destination PopHist.
 */
void PopHist_copy(PopHist * dest, const PopHist * src) {

    if(dest->size != src->size)
        eprintf("ERR in %s@%s:%d: src and dest have unequal sizes\n",
                __func__, __FILE__, __LINE__);

    memcpy(dest, src, src->size);

    PopHist_sanityCheck(dest, __FILE__, __LINE__);
}

/**
 * Fill ph with random values based on ph0.
 *
 * ph0 is first copied into ph. Then ph is perturbed as described in
 * the documentation to PopHist_perturbInPlace.
 *
 * @param[out] ph PopHist into which new values will be written.
 * @param[in] ph0 PopHist whose values are perturbed to obtain new
 * values. The values within ph0 are unchanged.
 * @param[in] dt controls the magnitude of perturbations in the
 * duration variable.
 * @param[in] dNinv controls the magnitude of perturbations in population
 * size.
 * @param[in] rng random number generator.
 * @returns ph
 */
PopHist    *PopHist_perturb(PopHist * ph, const PopHist * ph0, double dt,
                            double dNinv, gsl_rng * rng) {
    PopHist_copy(ph, ph0);
    return PopHist_perturbInPlace(ph, dt, dNinv, rng);
}

/**
 * Perturb parameters in PopHist structure ph.
 *
 * Each parameter is perturbed away from its initial value. The
 * distribution of the perturbations depends on the macros
 * PERTURB_GAUSSIAN and PERTURB_TDIST. If the first of these is
 * defined at compile time, then perturbations are
 * Gaussian. Otherwise, if the second macro is defined, perturbations
 * are drawn from a t distribution. Otherwise they are uniform.
 *
 * 1/2N values are reflected back and forth so that the perturbed
 * value lies within [loTwoNinv, hiTwoNinv].
 *
 * @param[in,out] ph The PopHist to be perturbed.
 * @param[in] dt controls the magnitude of perturbations in the
 * duration variable.
 * @param[in] dNinv controls the magnitude of perturbations in 1/2N.
 * @param[in] rng random number generator.
 * @returns ph
 */
PopHist    *PopHist_perturbInPlace(PopHist * ph, double dt, double dNinv,
                                   gsl_rng * rng) {
    double      x, curr;
    static const double loTwoNinv = 1.0e-8;
    static const double hiTwoNinv = 1.0;
    int         i;

    myassert(ph);

    /* Perturb 1/2N values */
    for(i = 0; i < ph->nepoch; ++i) {
        curr = PopHist_twoNinv(ph, i);
#ifdef PERTURB_GAUSSIAN
        x = gsl_ran_gaussian(rng, dNinv) + curr;
#elif defined(PERTURB_TDIST)
        x = dNinv * gsl_ran_tdist(rng, 1) + curr;
#else
        x = gsl_ran_flat(rng, curr - dNinv, curr + dNinv);
#endif
        x = reflect(x, loTwoNinv, hiTwoNinv);
        PopHist_setTwoNinv(ph, i, x);
    }

    /* Perturb t values */
    for(i = 0; i < ph->nepoch - 1; ++i) {
        curr = PopHist_duration(ph, i);
#ifdef PERTURB_GAUSSIAN
        x = gsl_ran_gaussian(rng, dt) + curr;
#elif defined(PERTURB_TDIST)
        x = dt * gsl_ran_tdist(rng, 1) + curr;
#else
        x = gsl_ran_flat(rng, curr - dt, curr + dt);
#endif
        PopHist_setDuration(ph, i, fabs(x));
    }

    PopHist_sanityCheck(ph, __FILE__, __LINE__);
    return ph;
}

/**
 * Measure the distance between two PopHist objects.
 */
double PopHist_distance(PopHist * ph1, PopHist * ph2) {
    int         i;
    int         nepoch = PopHist_nepoch(ph1);
    double      dist = 0.0;

    myassert(PopHist_nepoch(ph2) == nepoch);

    for(i = 0; i < nepoch; ++i) {
        double      a = PopHist_twoN(ph1, i);
        double      b = PopHist_twoN(ph2, i);

        dist += fabs(a - b);
        if(i + 1 < nepoch) {
            a = PopHist_duration(ph1, i);
            b = PopHist_duration(ph2, i);
            dist += fabs(a - b);
        }
    }
    return dist;
}

