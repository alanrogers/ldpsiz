/**
 * @file misc.c
 * @author Alan R. Rogers
 * @brief Miscellaneous functions.
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#ifdef _WIN32
#include <windows.h>
#elif defined(MACOS)
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif
#include <execinfo.h>
#include "misc.h"

// Divide num by denom and round the result to the nearest integer.
long LInt_div_round(long num, long denom) {
	assert(denom != 0);
	ldiv_t quotrem = ldiv(num, denom);
	if(2L * quotrem.rem > denom)
		return 1 + quotrem.quot;
	return quotrem.quot;
}

/**
 * Return 1 if string arg matches either shortOpt or longOpt. ShortOpt
 * and longOpt are ignored if their values are NULL. If arg doesn't
 * match either shortOpt or longOpt, then return 0.
 */
int isopt(const char *shortOpt, const char *longOpt, const char *arg) {
    size_t arglen = strlen(arg);
    
    if(shortOpt
       && arglen == strlen(shortOpt)
       && strncmp(arg, shortOpt, arglen) == 0)
        return 1;
    if(longOpt
       && arglen == strlen(longOpt)
       && strncmp(arg, longOpt, arglen) == 0)
        return 1;
    return 0;
}

/*
 * Describe an option. For use in "usage" functions.
 */
void tellopt(const char *opt, const char *description) {
    fprintf(stderr, "   %s\n      %s\n", opt, description);
    return;
}

/** Convert NULL-terminated string to lower case */
char       *strlowercase(char *s) {
    char       *p = s;

    for(p = s; *p != '\0'; ++p)
        *p = tolower(*p);
    return s;
}

/// Older C compilers lack strcasecmp, which compares strings ignoring
/// case. For portability, I provide this here under a different name.
#define ccmp(a,b) ((a) == (b) ? 0 : ((a) > (b) ? 1 : -1))
int mystrcasecmp(const char *s1, const char *s2) {
    char c1, c2;
    for ( ; ; ) {
       if (*s1 == '\0' || *s2 == '\0')
            return ccmp(*s1,*s2);
        c1= (isascii(*s1) && isupper(*s1)) ? tolower(*s1) : *s1;
        c2= (isascii(*s2) && isupper(*s2)) ? tolower(*s2) : *s2;
        if (c1 != c2)
            return ccmp(c1,c2);
        s1++;
        s2++;
    }
}
#undef ccmp

void checkmem( /*@null@ */ void *obj, const char *file, int line) {
    if(obj == NULL)
        die("allocation error", file, line);
    return;
}

void assertFiniteArray(const double *x, size_t dim, const char *file,
                       int line) {
    if(!isfinite(x[cblas_idamax(dim, x, 1)])) {
        fprintf(stderr, "ERR@%s:%d: Array is not finite: [", file, line);
        while(dim > 1) {
            fprintf(stderr, "% g,", *x++);
            --dim;
        }
        if(dim)
            fprintf(stderr, " %g", *x);
        fprintf(stderr, "]\n");
        dostacktrace(file, line, stderr);
        exit(EXIT_FAILURE);
    }
    return;
}

void printsqrmat(const char *msg, unsigned dim, double m[][dim]) {
    int         i, j;

    if(msg != NULL)
        printf("%s:\n", msg);
    for(i = 0; i < dim; ++i) {
        for(j = 0; j < dim; ++j)
            printf(" %12.3g", m[i][j]);
        putchar('\n');
    }
}

void printgslmat(const char *msg, size_t dim, gsl_matrix * m) {
    size_t      i, j;

    if(msg != NULL)
        printf("%s:\n", msg);
    for(i = 0; i < dim; ++i) {
        for(j = 0; j < dim; ++j)
            printf(" %8.4f", gsl_matrix_get(m, i, j));
        putchar('\n');
    }
}

int matIsFinite(unsigned dim, double m[][dim]) {
    unsigned    i, j;

    for(i = 0; i < dim; ++i)
        for(j = 0; j < dim; ++j)
            if(!isfinite(m[i][j]))
                return 0;
    return 1;
}

/*
 * Calculate the relative absolute difference between two vectors.
 */
double getreldiff(int dim, double x[], double y[], int verbose) {
    int         i;
    double      absdiff, abssum, relerr;

    abssum = absdiff = 0.0;
    for(i = 0; i < dim; ++i) {
        abssum += fabs(y[i]);
        absdiff += fabs(y[i] - x[i]);
        if(verbose)
            printf("%12.4g %12.4g\n", x[i], y[i]);
    }
    relerr = absdiff / abssum;

    if(verbose)
        printf("getreldiff: absdiff=%g abssum=%g relerr=%g\n",
               absdiff, abssum, relerr);

    return relerr;
}

/**
 * Center string "text" in a field of width "width". The centered string
 * is written into the character string "buff", whose size is
 * "buffsize".
 */
char       *strcenter(const char *text, unsigned width,
                      char *buff, unsigned buffsize) {
    int         i, j, lpad = 0, rpad = 0, txtwid;

    txtwid = strlen(text);
    if(txtwid >= buffsize) {
        snprintf(buff, buffsize, "%s", text);
        return buff;
    }
    if(width > txtwid)
        lpad = (width - txtwid) / 2;
    rpad = width - txtwid - lpad;
    for(i = 0; i < lpad; ++i)
        buff[i] = '-';
    for(j = 0; j < txtwid; ++j)
        buff[i + j] = text[j];
    for(j = 0; j < rpad; ++j)
        buff[i + txtwid + j] = '-';
    buff[lpad + txtwid + rpad] = '\0';
    return buff;
}

/*
 * An almost platform-independent function that returns the number of
 * CPU cores on the current machine.
 * Source: http://stackoverflow.com/questions/150355/
 * programmatically-find-the-number-of-cores-on-a-machine. If I
 * understand the webpage correctly, this code was written by Dirk-Jan
 * Kroon.
 */
int getNumCores(void) {
#ifdef WIN32
    SYSTEM_INFO sysinfo;

    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#elif defined(MACOS)
    int         nm[2];
    size_t      len = 4;
    uint32_t    count;

    nm[0] = CTL_HW;
    nm[1] = HW_AVAILCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);

    if(count < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &count, &len, NULL, 0);
        if(count < 1) {
            count = 1;
        }
    }
    return count;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

/*
 * Vector v must be sorted in ascending order before this function is
 * called.  Function returns index of first element in sorted array
 * that is >= val.  The function assumes without checking that the
 * input is sorted. If val > vec[len-1], the function returns len.
 */
long long_first_geq(long val, long *v, long len) {
    register long lo, mid, hi;

    myassert(len > 0);
    lo = 0;
    hi = len - 1;
    if(val > v[hi])
        return len;
    while(lo < hi) {
        mid = lo + (hi - lo) / 2;
        if(mid == lo)
            break;
        if(v[mid] < val)
            lo = mid;
        else
            hi = mid;
    }
    if(v[lo] >= val)
        hi = lo;

    myassert(hi >= 0);
    myassert(hi < len);
    myassert(v[hi] >= val);
    myassert(hi == 0 || v[hi - 1] < val);

    return hi;
}

/*
 * Vector v must be sorted in ascending order before this function is
 * called.  Function returns index of last element in sorted array
 * that is <= val.  The function assumes without checking that the
 * input is sorted. If val < vec[0], the function returns -1.
 */
long long_last_leq(long val, long *v, long len) {
    register long lo, mid, hi;

    myassert(len > 0);
    lo = 0;
    hi = len - 1;
    if(val < v[0])
        return -1;
    while(lo < hi) {
        mid = hi - (hi - lo) / 2;
        if(mid == hi)
            break;
        if(v[mid] > val)
            hi = mid;
        else
            lo = mid;
    }
    if(v[hi] <= val)
        lo = hi;

    myassert(lo >= 0);
    myassert(lo < len);
    myassert(v[lo] <= val);
    myassert(lo == len - 1 || v[lo + 1] > val);

    return lo;
}

/* print message and exit */
void die(const char *msg, const char *file, int line) {
    fflush(stdout);
    fprintf(stderr, "ERR@%s:%d: %s \n", file, line, msg);
    exit(EXIT_FAILURE);
}

/* eprintf: print error message and exit */
void eprintf(const char *fmt, ...) {
    va_list     args;

    fflush(stdout);

    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);

    if(fmt[0] != '\0' && fmt[strlen(fmt) - 1] == ':')
        fprintf(stderr, " %s", strerror(errno));
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}

/* duplicate memory block */
void       *memdup(const void *p, size_t n) {
    void       *q;

    myassert(p != NULL);
    myassert(n > 0);

    q = malloc(n);
    checkmem(q, __FILE__, __LINE__);
    memcpy(q, p, n);
    return q;
}

/*
 * In string str, count the number of contiguous chunks of characters
 * belonging to set.
 */
int strCountSetChunks(const char *str, const char *sep) {
    int         nchunks = 0, i;

    while(*str != '\0') {
        i = strcspn(str, sep);  /* skip chars not in set */
        str += i;
        i = strspn(str, sep);   /* skip chars in set */
        if(i > 0) {
            ++nchunks;
            str += i;
        }
    }
    return nchunks;
}

/**
 * Return true if the first non-white char in string s is '#'; false
 * otherwise.
 */
int strcomment(const char *s) {
    const char *p = s;

    while(isspace(*p))
        ++p;
    if(*p == '#')
        return true;
    return false;
}

/** strip comment ('#' to eol) from a string */
char       *stripComment(char *s) {
    char       *p = strchr(s, '#');

    if(p && (*p == '#'))
        *p = '\0';

    return s;
}

/** Return true if string contains only whitespace; false otherwise. */
int strempty(const char *s) {
    const char *p = s;

    while(isspace(*p))
        ++p;
    if(*p == '\0')
        return true;
    return false;
}

/**
 * Test equality of x and y. In this test NaN == NaN, Inf == Inf, and
 * -Inf == -Inf.
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
int dblEquals(double x, double y) {
    int         type_x, type_y;

    type_x = fpclassify(x);
    type_y = fpclassify(y);

    if(type_y != type_x)
        return 0;               /* x != y because types differ */

    assert(type_x == type_y);

    if(type_y == FP_NAN)
        return 1;

    return (x == y);
}

#pragma GCC diagnostic pop

#define CALLSTACK_SIZE 128
void dostacktrace(const char *file, int line, FILE * ofp) {
    void       *callstack[CALLSTACK_SIZE];
    int         nsymbols = backtrace(callstack, CALLSTACK_SIZE);

    fprintf(ofp, "backtrace returned %d\n", nsymbols);
    fprintf(ofp, "dostacktrace called from %s:%d:\n", file, line);
    backtrace_symbols_fd(callstack, nsymbols, fileno(ofp));
}

/**
 * Fold x back and forth across the boundaries "lo" and "hi" to obtain a value
 * y such that lo <= y <= hi.
 */
double reflect(double x, double lo, double hi) {
    assert(hi > lo);
    x -= lo;
    hi -= lo;

    double      z = fabs(fmod(x, 2.0 * hi));

    /*    printf("initially z=%lg\n", z); */

    if(z > hi)
        z = 2.0 * hi - z;

    z += lo;

    assert(z >= lo && z <= hi + lo);

    return z;
}

/**
 * Replace suffix (part after last period) of string "str" with new
 * suffix. If original string has no suffix, then append one to
 * string. Abort if string isn't long enough.
 */
void replaceSuffix(char *str, size_t str_size, const char *suffix, size_t suffix_len) {
    char *s = strrchr(str, '.');
    if(s == NULL)
        s = strrchr(str, '\0');
    str_size -= s-str;
    if(str_size < suffix_len+1)
        eprintf("%s:%s:%d: string is too small for suffix\n",
                __func__, __FILE__, __LINE__);
    snprintf(s, str_size,"%s", suffix);
    
}

void unitTstResult(const char *facility, const char *result) {
    printf("%-26s %s\n", facility, result);
}

/**
 * Compare two long ints.
 *
 * Function interprets its arguments as pointers to long ints.
 *
 * @param void_x,void_y Pointers to the two integers, cast as pointers
 * to voids.
 * @returns -1, 0, or 1 depending on whether the first arg is <,
 * ==, or > the second.
 */
int compareLongs(const void *void_x, const void *void_y) {
    const long *x = (const long *) void_x;
    const long *y = (const long *) void_y;

    return (*x < *y) ? -1 : (*x > *y) ? 1 : 0;
}

/**
 * Compare two doubles.
 *
 * Function interprets its arguments as pointers to doubles.
 *
 * @param void_x,void_y Pointers to the two doubles, cast as pointers
 * to voids.
 * @returns -1, 0, or 1 depending on whether the first arg is <,
 * ==, or > the second.
 */
int compareDoubles(const void *void_x, const void *void_y) {
    const double *x = (const double *) void_x;
    const double *y = (const double *) void_y;

    return (*x < *y) ? -1 : (*x > *y) ? 1 : 0;
}

int pr_gsl_vector(FILE *fp, const char *fmt, const gsl_vector * v) {
    int rval=0;
    size_t i, n = v->size;
    putc('[', fp);
    if(n>0) {
        rval = fprintf(fp, fmt, gsl_vector_get(v, 0));
        if(rval <= 0)
            return rval;
        
    }
    for(i=1; i < n; ++i) {
        fputs(", ", fp);
        rval = fprintf(fp,fmt,  gsl_vector_get(v, i));
        if(rval <= 0)
            return rval;
    }
    putc(']', fp);
    return rval;
}

/// Hash a character string. Used for generating unique file names.
unsigned hash(const char *s) {
    unsigned hashval;

    for(hashval=0; *s != '\0'; ++s)
        hashval += *s + 31 * hashval;

    return hashval % 1001;
}
