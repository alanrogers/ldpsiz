#include <stdio.h>
#include <string.h>
int main(void) {
    const long nSamples = 8;
    const long dim = nSamples - 1;
    long i, j, k, ii, jj;
    double rvec[dim][dim];
    double cvec[dim][dim];
    double A[dim][dim];

    memset(rvec, 0, sizeof rvec);
    memset(cvec, 0, sizeof cvec);
    memset(A, 0, sizeof cvec);

    for(j=2; j <= nSamples; ++j) {
        jj = j-2;
        double lambda = -j*(j-1);
        cvec[jj][jj] = rvec[jj][jj] = 1;
        for(i=j-1; i > 1; --i) {
            ii = i-2;
            cvec[ii][jj] = cvec[ii+1][jj]*(i*(i+1.0)/(i*(i-1.0) + lambda));
        }
        for(i=j+1; i <= nSamples; ++i) {
            ii = i-2;
            rvec[jj][ii] = rvec[jj][ii-1]*(i*(i-1.0)/(i*(i-1.0) + lambda));
        }
    }

    printf("Row eigenvectors:\n");
    for(ii=0; ii<dim; ++ii) {
        for(jj=0; jj<dim; ++jj)
            printf(" %8.4lf", rvec[ii][jj]);
        putchar('\n');
    }

    printf("Column eigenvectors:\n");
    for(ii=0; ii<dim; ++ii) {
        for(jj=0; jj<dim; ++jj)
            printf(" %8.4lf", cvec[ii][jj]);
        putchar('\n');
    }

    for(i=0; i<dim; ++i)
        for(j=0; j<dim; ++j)
            for(k=0; k<dim; ++k)
                A[i][j] += cvec[i][k] * rvec[k][j];
        
    printf("Product should be identity matrix:\n");
    for(ii=0; ii<dim; ++ii) {
        for(jj=0; jj<dim; ++jj)
            printf(" %8.4lf", A[ii][jj]);
        putchar('\n');
    }
    return 0;
}
