//usr/bin/gcc

//Creating a 1D lattice with entries of 0,1 and -1
// 1 lables Particles that move right
// -1 lables Particles that move left

//--------------Libaries-----------------------
#include <stdio.h>
#include <stdlib.h>
#include "./../Lib/random.h"
#include <math.h>
#include <string.h>
//--------------Prototypes-----------------------
int * init_lat(int N,float phi);
int rand_index(int arraylength);
//----------------Main Program--------------------
void main() {
	printf("--------------------\n");
	printf("1D-Testing started\n");
	float alph, phi;			//alph: propability for tumbling event; phi: particel concentration
	int M,N;				//M: total number of particles; N: total number of sites (or length of lattice array)
	char word;
	printf("Number of sites (N): ");
	scanf("\n%d", &N);			//get number of sites
	printf("Particle concentration (phi): ");
	scanf("\n%f", &phi);			//get concentration
	printf("Probability for tumbling (alpha): ");
	scanf("\n%f", &alph);			//get tumbling probability
	int lattice[N];				//create lattice with N sites
	int i=0;			
	for(i=0;i<N;i++) {
		printf("%d\n",lattice[i]);
	}
	init_lat(N,phi);
}
//-----------------Functions-----------------------

//--init_lat--
//initialises the Lattice with cells on it
//input:int N (number of sites),float phi (particle concentration)
//output: int array lattice (lattice with particles at timestep 0)
int * init_lat(int N,float phi) {
	long int seed = 123456789;
	int lattice[N];				
	float M_ =(float)(N)*phi; 				//M (number of Particles) --> if N*phi >= n.5 (with n natrual number) there is an error. This error is negligible for big N
	int M = roundf(M_);
	printf("%d*%f= %d",N,phi,M);
	double interval=1/N;			//separate the space 0-1 into N pieces with length interval
	int i = 0;
	double rndnum;
	int len = strlen(lattice);
	printf("length of array:  %d",len);			
	for(i=0;i<M;i++) { 			//loop to find random indizes
		rand_index(lattice)
		 //here is a problem, i wont get the right index cause i cut something off... need to fix this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}
	printf("chosen indizes: %d",lattice);
	return lattice
}

//--rand_index--
//returns a random index of given array
int rand_index(int arraylenght) {
	long int seed = 123456789;
	double interval=1/arraylength;		//separate the space 0-1 into N pieces with length interval
	int i = 0;
	double rndnum;			
	rndnum = ran3(&seed);			//pick a random number out of the spacing 0-1
	index = roundf(rndnum/interval);	//evaluate the index of the lattice, this index is a randomly chosen one		
}

