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
#include <time.h>
//#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>
//--------------Prototypes-----------------------
int * init_lat(int N,int M,float phi);
int rand_index(double arraylength);
int length(int * array);
int rnddirection(); 
int * timestep(int * lattice,int N, int M, double alph);
char * geturand();
//----------------Main Program--------------------
int main() 
{
	geturand();
	printf("--------------------\n");
	printf("1D-Testing started\n");
	float alph, phi;			//alph: propability for tumbling event; phi: particle concentration
	int M,N,tottime;			//M: total number of particles; N: total number of sites (or length of lattice array)
	char word;
	printf("Number of sites (N): ");
	scanf("\n%d", &N);			//get number of sites
	printf("Particle concentration (phi): ");
	scanf("\n%f", &phi);			//get concentration
	printf("Probability for tumbling (alpha): ");
	scanf("\n%f", &alph);			//get tumbling probability
	printf("Total Time for evolution (T): ");
	scanf("\n%d",&tottime);
	float M_ =(float)(N)*phi; 				//M (number of Particles) --> if N*phi >= n.5 (with n natrual number) there is an error. This error is negligible for big N
	M = roundf(M_);
	int *lattice;				//declare lattice 
	int i=0,ii=0;			
	lattice = init_lat(N,M,phi);
	for(ii=0;ii<tottime;ii++)
	{
		for(i=0;i<N;i++) 
		{
			printf("%d",*(lattice+i));
		}
		printf("\n");
		lattice=timestep(lattice,N,M,alph);
	}
}
//-----------------Functions-----------------------

//--init_lat--
//initialises the Lattice with cells on it
//input:int N (number of sites),float phi (particle concentration)
//output: int array lattice (lattice with particles at timestep 0)
int * init_lat(int N,int M,float phi) 
{
	long int seed = 123456789;
	static int lattice[3000];				//allocating 3000*sizeof(integer)bits space for the lattice array --> should be allocated dynamically, but didnt work till now	
//	printf("%d*%f= %d",N,phi,M);
	double interval=1/N;			//separate the space 0-1 into N pieces with length interval
	int i = 0;
	double rndnum;
	int len = N;
	int ind;
	printf("length of array:  %d\n",len);			
	for(i=0;i<M;i++)  			//loop to find random indizes
	{
		ind=rand_index(len);
//		printf("chosen indizes: %d\n",ind);
		if(lattice[ind]==0) 
		{
			lattice[ind]=rnddirection();
		}
		else
		{
			i--;
		}
	}
	return lattice;
}

//--rand_index--
//returns a random index between 0 and arraylength
//input: double arraylength
//output: int random index
int rand_index(double arraylength) 
{
	long int seed = 123456789;
	double interval = 1/arraylength;	//separate the space 0-1 into N pieces with length interval
//	printf("interval: %f\n",interval);
	int i = 0;
	double rndnum;			
	rndnum = ran3(&seed);			//pick a random number out of the spacing 0-1
	int index = roundf(rndnum/interval);	//evaluate the index of the lattice, this index is a randomly chosen one		
	return index;
}
//--direction--
//randomly chooses a -1 or 1
//input none
//output: -1 or 1 as int
int rnddirection() 
{		
	time_t t= time(0);
	double rndnum;
	long int seed=t;
	printf("seed = %ld\n",seed);
	rndnum = ran3(&seed);
	if (rndnum < 0.5)
	{
//		printf("-1");
		return -1;
	}
	else 
	{
//		printf("1");
		return 1;
	}
}
//--timestep--
//input: int pointer to first element of lattice (allocated in init)
//output: int * lattice (evolution of lattice after one timestep)
int * timestep(int * lattice,int N,int M, double alph)
{
	int i = 0;
	int ind;
//	printf("M= %d\n",M);
	for(i = 0; i < M; i++)
	{
		ind=rand_index(N);
//		printf("timestep .. lattice[ind]= %d\n",lattice[ind]);
		if(lattice[ind] == 0) //if there is no particle try it again
		{
			i--;
		}
		else
		{
			if(ind == N-1) //if upper periodic boundary 
			{
				if(lattice[ind] > 0) 		//find out the direction
				{
					if(lattice[0] == 0) 	//find out if "the way if free"
					{
						lattice[ind]=0; //move
						lattice[0]=1;
					}
				}
				else
				{
					if(lattice[ind-1] == 0)	//if the way is free
					{
						lattice[ind]=0;		//move
						lattice[ind-1]=-1;
					}	
				}
			}
			else
			{
				if(ind == 0)	//if lower periodic boundary
				{
					if(lattice[ind] > 0) 			//find out the direction
					{
						if(lattice[ind+1] == 0) 	//find out if "the way if free"
						{
							lattice[ind]=0;		//move
							lattice[ind+1]=1;
						}
					}
					else
					{
						if(lattice[N-1] == 0) 		//if the way if free
						{
							lattice[ind]=0;		//move
							lattice[N-1]=-1;
						}	
					
					}
				}
				else
				{
					if(lattice[ind] > 0) 			//find out the direction
					{
						if(lattice[ind+1] == 0) 	//find out if "the way if free"
						{
							lattice[ind]=0;		//move
							lattice[ind+1]=1;
						}
					}
					else	
					{
						if(lattice[ind-1] == 0)		//if the way is free
						{
							lattice[ind]=0;		//move
							lattice[ind-1]=-1;
						}
					}	
				}
			}
		}
	}
	return lattice;
}
/*
char * geturand()
{
	int randomData = open("/dev/random", O_RDONLY);		//opens a random and returns int reffering to random
	static char myRandomData[50];
	size_t randomDataLen = 0;
	while (randomDataLen < sizeof myRandomData)
	{
	    ssize_t result = read(myRandomData, myRandomData + randomDataLen, (sizeof myRandomData) - randomDataLen);
	    if (result < 0)
	    {
		// error, unable to read /dev/random 
	    }
	    randomDataLen += result;
	}
	close(randomData);
//	printf("",myrandomdata);
	return myRandomData;
}
*/
char * geturand()
{
	char ch,file_name[25];
	FILE *fp;
	int len = 10;
	//char ranch[len];
	//  printf("Enter the name of file you wish to see\n");
	//file_name="/dev/urand";

	fp = fopen("/dev/urandom","r"); // read mode

	if( fp == NULL )
	{
		perror("Error while opening the file.\n");
		exit(EXIT_FAILURE);
	}

	//printf("The contents of %s file are :\n", file_name);
	int i;
	for(i=0;i<len;i++)
	{
		ch = fgetc(fp);
		printf("%d",ch);
	}
	fclose(fp);
	return 0;
}
