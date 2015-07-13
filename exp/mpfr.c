#include <stdio.h>
#include <stdarg.h>
#include <mpfr.h>
#include <math.h>

int main(void) {
    // MPFR_RNDN : round to nearest
    // MPFR_RNDZ : round toward zero
    // mpfr_t : floating point type
    // mpfr_prec_t : floating point precision type
    // mpfr_free_cache : call before terminating thread
    // mpfr can be compiled as thread safe (TLS for thread-local storage)

    mpfr_t x, y, z;
    double xd, yd, zd;
    //mpfr_prec_t prec;
    mpfr_rnd_t rnd = MPFR_RNDN; // round to nearest
#if 1
    mpfr_inits2(256, x, y, z, (mpfr_ptr) 0);

    mpfr_set_d(x, 2.0, rnd);
    mpfr_set_d(y, 3.0, rnd);
    xd = mpfr_get_d(x, rnd);
    yd = mpfr_get_d(y, rnd);

    mpfr_add(z, x, y, rnd);
    zd = mpfr_get_d(z, rnd);
    printf("Add: %lf+%lf = %lf\n", xd, yd, zd);

    mpfr_sub(z, x, y, rnd);
    zd = mpfr_get_d(z, rnd);
    printf("Subtr: %lf-%lf = %lf\n", xd, yd, zd);

    mpfr_mul(z, x, y, rnd);
    zd = mpfr_get_d(z, rnd);
    printf("Mult: %lf*%lf = %lf\n", xd, yd, zd);

    mpfr_div(z, x, y, rnd);
    zd = mpfr_get_d(z, rnd);
    printf("Div: %lf/%lf = %lf\n", xd, yd, zd);

    mpfr_pow(z, x, y, rnd);
    zd = mpfr_get_d(z, rnd);
    printf("Pow: %lf^%lf = %lf\n", xd, yd, zd);

    mpfr_sqr(z, y, rnd);
    zd = mpfr_get_d(z, rnd);
    printf("Sqr: %lf^2 = %lf\n", yd, zd);

    mpfr_sqrt(z, y, rnd);
    zd = mpfr_get_d(z, rnd);
    printf("Sqrt: sqrt(%lf) = %lf\n", yd, zd);

    mpfr_neg(z, y, rnd);
    zd = mpfr_get_d(z, rnd);
    printf("Neg: -%lf = %lf\n", yd, zd);

    mpfr_abs(z, y, rnd);
    zd = mpfr_get_d(z, rnd);
    printf("Abs: abs(%lf) = %lf\n", yd, zd);

    mpfr_exp(z, y, rnd);
    zd = mpfr_get_d(z, rnd);
    printf("Exp: exp(%lf) = %lf\n", yd, zd);

    mpfr_log(z, y, rnd);
    zd = mpfr_get_d(z, rnd);
    printf("Log: log(%lf) = %lf\n", yd, zd);

    // Write z in full precision
    size_t charsWritten = mpfr_out_str(stdout, 10, 0, z, rnd);
    putchar('\n');
    printf("Wrote %zu chars\n", charsWritten);

    // Write 25 significant digits
    charsWritten = mpfr_out_str(stdout, 10, 25, z, rnd);
    putchar('\n');
    printf("Wrote %zu chars\n", charsWritten);

#endif
    mpfr_clears(x, y, (mpfr_ptr) 0);
    mpfr_free_cache();
    return 0;
}
