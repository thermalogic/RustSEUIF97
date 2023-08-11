/*
   The C example to call the shared library
     gcc demo.c -o demo  -L../target/release  -Wl,-rpath=../target/release  -lseuif97 -lm
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define OH 4
#define OS 5

extern double pt(double p,double t,short o_id);

int main(void)
{

    double p = 16.0;
    double t = 530.0;
    double h, s;

    h = pt(p, t, OH);
    s = pt(p, t, OS);
    printf("p,t %f,%f h= %f s= %f\n", p, t, h, s);
    return EXIT_SUCCESS;
}
