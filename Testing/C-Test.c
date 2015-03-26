//usr/bin/gcc

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
void main() {
	printf("learning c again\n");
	int x,i ;
	int y;
	y=0;
	for(i = 0; i < 200; i++ ) {
		srand( time(NULL)+i);
		x = rand() %100 ;
		printf("random integer: %d \n",x);
		y=y+x;
	}
	printf("y: %i",y);
	printf("<y>: %f\n", y/200.0);
}

