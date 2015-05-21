/**
 * @file xsnp.c
 * @author Alan R. Rogers
 * @brief Test snp.c.
 * @copyright Copyright (c) 2015, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "snp.h"
#include <stdio.h>
#include <string.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {
    int         verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xsnp [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xsnp [-v]\n");
    }

    SNP_test(verbose);

    return 0;
}
