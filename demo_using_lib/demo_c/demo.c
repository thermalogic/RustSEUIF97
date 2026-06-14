
/*
   The C example to call the shared library: cdel
    cargo build -r --features cdecl
    make  
*/

#include <stdlib.h>
#include <stdio.h>


#define OH 4

extern double pt(double p,double t,short o_id);
extern double pt2s(double p,double t);

int main(void)
{

    double p = 16.0;
    double t = 530.0;
    double h, s;
    // universal functions with o_id parameter
    h = pt(p, t, OH);
    // direct property functions
    s = pt2s(p, t);
    printf("p,t %f,%f h= %f s= %f\n", p, t, h, s);
    return EXIT_SUCCESS;
}
