/**
 * @file xmisc.c
 * @author Alan R. Rogers
 * @brief Test misc.c
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "misc.h"
#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include <float.h>

int main(int argc, char **argv) {

    int         verbose = 0;

    long        v[] = { 0, 0, 1, 1, 1, 2, 2 };
    long        len = 7;
    const char *str1 = "abc1cd23efgh4";
    const char *str2 = "999999999abc1cd23efgh4";
    const char *str3 = "abc1cd23efgh";
    const char *str4 = "abccdefgh";
    const char *set = "0123456789";

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xmisc [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xmisc [-v]\n");
    }

    assert(dblEquals(0.0 / 0.0, 0.0 / 0.0));
    assert(!dblEquals(0.0 / 0.0, DBL_MAX));
    assert(dblEquals(1.0 / 0.0, 1.0 / 0.0));
    assert(!dblEquals(-1.0 / 0.0, 1.0 / 0.0));
    assert(dblEquals(1.23, 1.23));
    unitTstResult("dblEquals", "OK");

    assert(long_last_leq(-1, v, len) == -1);
    assert(long_last_leq(0, v, len) == 1);
    assert(long_last_leq(1, v, len) == 4);
    assert(long_last_leq(2, v, len) == 6);
    assert(long_last_leq(3, v, len) == 6);
    assert(long_last_leq(-1, v, 1) == -1);
    assert(long_last_leq(0, v, 1) == 0);
    assert(long_last_leq(1, v, 1) == 0);
    unitTstResult("long_last_leq", "OK");

    assert(long_first_geq(-1, v, len) == 0);
    assert(long_first_geq(0, v, len) == 0);
    assert(long_first_geq(1, v, len) == 2);
    assert(long_first_geq(2, v, len) == 5);
    assert(long_first_geq(3, v, len) == 7);
    assert(long_first_geq(-1, v, 1) == 0);
    assert(long_first_geq(0, v, 1) == 0);
    assert(long_first_geq(1, v, 1) == 1);
    unitTstResult("long_first_geq", "OK");

    assert(strCountSetChunks(str1, set) == 3);
    assert(strCountSetChunks(str2, set) == 4);
    assert(strCountSetChunks(str3, set) == 2);
    assert(strCountSetChunks(str4, set) == 0);
    assert(strCountSetChunks("", set) == 0);
    assert(strCountSetChunks(set, set) == 1);
    unitTstResult("strCountSetChunks", "OK");

    double      x[] = { 1.0, 3.0, 4.3 };
    assertFiniteArray(x, 3, __FILE__, __LINE__);
    unitTstResult("assertFiniteArray", "OK");

#if 0
    x[1] = 1.0 / 0.0;
    printf("next line should die\n");
    assertFiniteArray(x, 3, __FILE__, __LINE__);
#endif

    assert(strcomment("   ab cde") == 0);
    assert(strcomment("   #ab cde") == 1);
    unitTstResult("strcomment", "OK");

    char        buff[30], *s;

    snprintf(buff, sizeof(buff), " asdfaf #comment");

    s = stripComment(buff);
    assert(s == buff);
    assert(strcmp(buff, " asdfaf ") == 0);
    assert(strlen(buff) == strlen(" asdfaf "));
    unitTstResult("stripComment", "OK");

    assert(reflect(-0.5, 1.0, 2.0) == 1.5);
    assert(reflect(0.0, 1.0, 2.0) == 2.0);
    assert(reflect(1.5, 1.0, 2.0) == 1.5);
    assert(reflect(1.0, 1.0, 2.0) == 1.0);
    assert(reflect(2.0, 1.0, 2.0) == 2.0);
    assert(reflect(2.25, 1.0, 2.0) == 1.75);
    assert(reflect(3.25, 1.0, 2.0) == 1.25);
    assert(reflect(4.75, 1.0, 2.0) == 1.25);
    unitTstResult("reflect", "OK");

    assert(encode01('0') == 0);
    assert(encode01('1') == 1);
    assert(encode01('2') == 255);
    assert(encode01('h') == 255);
    unitTstResult("encode01", "OK");

    unsigned char gtype[100];
    unsigned    i, nGtype;
    const char *gtypeString = "1001h1100";

    nGtype = encodeDiploid(gtype, sizeof(gtype), gtypeString);
    if(verbose) {
        printf("diploid input: %s\n", gtypeString);
        printf("nGtype=%u\ngtype:", nGtype);
        for(i = 0; i < nGtype; ++i)
            printf(" %u", (unsigned) gtype[i]);
        putchar('\n');
    }
    assert(nGtype == 5);
    assert(gtype[0] == 2);
    assert(gtype[1] == 1);
    assert(gtype[2] == UNPHASED_HETEROZYGOTE);
    assert(gtype[3] == 3);
    assert(gtype[4] == 0);

    unitTstResult("encodeDiploid", "OK");

    gtypeString = "10011";
    nGtype = encodeHaploid(gtype, sizeof(gtype), gtypeString);
    if(verbose) {
        printf("haploid input: %s\n", gtypeString);
        printf("nGtype=%u\ngtype:", nGtype);
        for(i = 0; i < nGtype; ++i)
            printf(" %u", (unsigned) gtype[i]);
        putchar('\n');
    }
    assert(nGtype == 5);
    assert(gtype[0] == 1);
    assert(gtype[1] == 0);
    assert(gtype[2] == 0);
    assert(gtype[3] == 1);
    assert(gtype[4] == 1);
    unitTstResult("encodeHaploid", "OK");

    char suffix[FILENAMESIZE] = {'\0'};
    char s1[FILENAMESIZE] = {'\0'};
    char s2[FILENAMESIZE] = {'\0'};

    snprintf(suffix, sizeof(suffix),".bar");
    snprintf(s1, sizeof(suffix),"pear.apple");
    snprintf(s2, sizeof(suffix),"bananna");

    replaceSuffix(s1, sizeof(s1), suffix, strlen(suffix));
    if(verbose)
        printf("\"%s\" should equal \"pear.bar\"\n", s1);
    assert(strcmp(s1, "pear.bar") == 0);

    replaceSuffix(s2, sizeof(s2), suffix, strlen(suffix));
    if(verbose)
        printf("\"%s\" should equal \"bananna.bar\"\n", s2);
    assert(strcmp(s2, "bananna.bar") == 0);

    unitTstResult("replaceSuffix", "OK");

    assert(0 == mystrcasecmp("AbCd", "aBcD"));
    assert(0 != mystrcasecmp("AbCde", "aBcD"));
    unitTstResult("mystrcasecmp", "OK");

	assert(3 == LInt_div_round(9L, 3L));
	assert(3 == LInt_div_round(10L, 3L));
	assert(4 == LInt_div_round(11L, 3L));
	assert(4 == LInt_div_round(12L, 3L));
    unitTstResult("LInt_div_round", "OK");

    if(verbose) {
        printf("hash(%s) = %04x\n", "a", hash("a"));
        printf("hash(%s) = %04x\n", "abc", hash("abc"));
        printf("hash(%s) = %04x\n", "abcxxxxxffff", hash("abcxxxxxffff"));
    }
    unitTstResult("hash", "OK");

    return 0;
}
