#include <stdio.h>
#include <float.h>
int main(void) {
    printf("%a\n", 1.0/7.0);
    printf("%.*le\n", DBL_DIG + 3, 1.0/7.0);
    printf("%*s\n", DBL_DIG + 9, "asf");
    return 0;
}
