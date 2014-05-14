/**
 * @file xassign.c
 * @author Alan R. Rogers
 * @brief Test assign.c.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "assign.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {

    int         ok = 1;
    int         verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xassign [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xassign [-v]\n");
    }

    double      dbl1 = 1.0, dbl2 = -999.123, dbl3 = 3.456789e10;
    int         int1 = 1, int2 = -999, int3 = 345678910;
    long        longval;
    const char *str1 = "hello";
    const char *str2 = "123";
    const char *str3 =
        "this is a long string this is a long string this is a long string "
        "this is a long string this is a long string this is a long string "
        "this is a long string this is a long string this is a long string "
        "this is a long string this is a long string this is a long string "
        "this is a long string this is a long string this is a long string "
        "this is a long string this is a long string this is a long string "
        "this is a long string this is a long string this is a long string "
        "this is a long string this is a long string this is a long string "
        "this is a long string this is a long string this is a long string";

    double      x = 0;
    int         rval, i, cmp;
    char        buff[1000];

    memset(buff, 0, sizeof(buff));

    Assignment *a = NULL;

    if(verbose) {
        printf("Before assignments, list is:\n");
        Assignment_print(a, stdout);
    }

    snprintf(buff, sizeof(buff), "%lf", dbl1);
    a = Assignment_new(a, "dbl1", buff);

    snprintf(buff, sizeof(buff), "%lf", dbl2);
    a = Assignment_new(a, "dbl2", buff);

    snprintf(buff, sizeof(buff), "%lf", dbl3);
    a = Assignment_new(a, "dbl3", buff);

    snprintf(buff, sizeof(buff), "%d", int1);
    a = Assignment_new(a, "int1", buff);

    snprintf(buff, sizeof(buff), "%d", int2);
    a = Assignment_new(a, "int2", buff);

    snprintf(buff, sizeof(buff), "%d", int3);
    a = Assignment_new(a, "int3", buff);

    a = Assignment_new(a, "str1", str1);
    a = Assignment_new(a, "str2", str2);
    a = Assignment_new(a, "str3", str3);

    if(verbose) {
        printf("After assignments, list is:\n");
        Assignment_print(a, stdout);
    }

    rval = Assignment_setDbl(a, "dbl1", &x, MANDATORY);
    assert(rval == 1);
    if(x != dbl1) {
        ok = 0;
        printf("Bad assignment@%s:%d: key=%s true=%lf got=%lf\n",
               __FILE__, __LINE__, "dbl1", dbl1, x);
    }

    rval = Assignment_setDbl(a, "dbl2", &x, MANDATORY);
    assert(rval == 1);
    if(x != dbl2) {
        ok = 0;
        printf("Bad assignment@%s:%d: key=%s true=%lf got=%lf\n",
               __FILE__, __LINE__, "dbl2", dbl2, x);
    }

    rval = Assignment_setDbl(a, "dbl3", &x, MANDATORY);
    assert(rval == 1);
    if(x != dbl3) {
        ok = 0;
        printf("Bad assignment@%s:%d: key=%s true=%lf got=%lf\n",
               __FILE__, __LINE__, "dbl3", dbl3, x);
    }

    rval = Assignment_setInt(a, "int1", &i, MANDATORY);
    assert(rval == 1);
    if(i != int1) {
        ok = 0;
        printf("Bad assignment@%s:%d: key=%s true=%d got=%d\n",
               __FILE__, __LINE__, "int1", int1, i);
    }

    rval = Assignment_setInt(a, "int2", &i, MANDATORY);
    assert(rval == 1);
    if(i != int2) {
        ok = 0;
        printf("Bad assignment@%s:%d: key=%s true=%d got=%d\n",
               __FILE__, __LINE__, "int2", int2, i);
    }

    rval = Assignment_setInt(a, "int3", &i, MANDATORY);
    assert(rval == 1);
    if(i != int3) {
        ok = 0;
        printf("Bad assignment@%s:%d: key=%s true=%d got=%d\n",
               __FILE__, __LINE__, "int3", int3, i);
    }

    rval = Assignment_setLong(a, "int3", &longval, MANDATORY);
    assert(rval == 1);
    if(longval != int3) {
        ok = 0;
        printf("Bad assignment@%s:%d: key=%s true=%d got=%ld\n",
               __FILE__, __LINE__, "int3", int3, longval);
    }

    rval = Assignment_setString(a, "str1", buff, sizeof(buff), MANDATORY);
    assert(rval == 1);
    cmp = strcmp(buff, str1);
    if(cmp != 0) {
        ok = 0;
        printf("Bad assignment@%s:%d: key=%s true=\"%s\" got=\"%s\"\n",
               __FILE__, __LINE__, "str1", str1, buff);
    }

    rval = Assignment_setString(a, "str2", buff, sizeof(buff), MANDATORY);
    assert(rval == 1);
    cmp = strcmp(buff, str2);
    if(cmp != 0) {
        ok = 0;
        printf("Bad assignment@%s:%d: key=%s true=\"%s\" got=\"%s\"\n",
               __FILE__, __LINE__, "str2", str2, buff);
    }

    rval = Assignment_setString(a, "str3", buff, sizeof(buff), MANDATORY);
    assert(rval == 1);
    cmp = strcmp(buff, str3);
    if(cmp != 0) {
        ok = 0;
        printf("Bad assignment@%s:%d: key=%s true=\"%s\" got=\"%s\"\n",
               __FILE__, __LINE__, "str3", str3, buff);
    }

    x = 123.4;
    rval = Assignment_setDbl(a, "notthere", &x, !MANDATORY);
    assert(x = 123.4);
    assert(rval == 0);

    i = 1234;
    rval = Assignment_setInt(a, "notthere", &i, !MANDATORY);
    assert(i = 1234);
    assert(rval == 0);

    longval = 1234;
    rval = Assignment_setLong(a, "notthere", &longval, !MANDATORY);
    assert(longval = 1234);
    assert(rval == 0);

    snprintf(buff, sizeof(buff), "1234");
    rval =
        Assignment_setString(a, "notthere", buff, sizeof(buff), !MANDATORY);
    assert(strcmp(buff, "1234") == 0);
    assert(rval == 0);

    Assignment_free(a);

    if(ok) {
        unitTstResult("Assignment", "OK");
    } else
        unitTstResult("Assignment", "FAIL");

    return (ok ? 0 : EXIT_FAILURE);
}

/* Local Variables: */
/* mode: c */
/* End: */
