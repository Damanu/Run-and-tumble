//usr/bin/gcc
//#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "./../Lib/random.h"
void main() {
	printf("learning c again\n");
	double x;
	int i ;
	double y;
	long int seed;
	seed = 1234567;
	y=0;
	for(i = 0; i < 200; i++ ) {
		x = ran3(&seed);
		printf("random double [0.1]: %f \n",x);
		y=y+x;
	}
	printf("y: %f\n",y);
	printf("<y>: %f\n", y/200.0);
}

