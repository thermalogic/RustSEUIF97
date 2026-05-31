
/*
   The C example to call the shared library
   Windows with MinGW-W64：
      gcc -o demo.exe demo.c -L../target/release/  -lseuif97
   Linux
      gcc demo.c -o demo  -L../target/release  -Wl,-rpath=../target/release  -lseuif97 -
     
    make  

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define OH 4

extern double pt(double p,double t,short o_id);
extern double pt2s(double p,double t);

int main(void)
{

    double p = 16.0;
    double t = 530.0;
    double h, s;

    h = pt(p, t, OH);
    s = pt2s(p, t);
    printf("p,t %f,%f h= %f s= %f\n", p, t, h, s);
    return EXIT_SUCCESS;
}
