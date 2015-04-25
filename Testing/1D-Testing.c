//usr/bin/gcc

//Creating a 1D lattice with entries of 0,1 and -1
// 1 lables Particles that move right
// -1 lables Particles that move left


//dev/urandom ... rand nummer gen linux

//--------------Libaries-----------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#include <fcntl.h>
//#include <stdio.h>
//#include <time.h>
#include <stdlib.h>
#include <stdlib.h>
#include "./../Lib/random.h"
//#include <unistd.h>
//--------------Prototypes-----------------------
int * init_lat(int N,int M,float phi);
int rand_index(double arraylength);
int length(int * array);
int rnddirection(); 
int * timestep(int * lattice,int N, int M, double alph);
//char * geturand();
int tumble(int dir,double alph);
int cluster_counter(int * lattice,int N);
//--------------Main Program---------------------
int main() 
{
	srand(time(0));
	//double r = rand()/(double)RAND_MAX;
	long t=time(0);
	double r=ran3(&t);
	printf("%lf",r);
	printf("\n--------------------\n");
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
	int i,ii=0;			
	lattice = init_lat(N,M,phi);		//initalize lattice
	for(ii=0;ii<tottime;ii++)
	{
		for(i=0;i<N;i++) 
		{
			//printf("%d",*(lattice+i));
			if(lattice[i]==0)
			{
				printf(" ");
			}
			else
			{		
				printf("|");
			}
		}
		printf("\n");
		lattice=timestep(lattice,N,M,alph);
	}
	int cl_count = cluster_counter(lattice,N);	
	double Lc = (double)M/(double)cl_count;		//Lc is the average cluster size .... mistake: does not consider particles that are not in a cluster (all particles are in cluster ... wrong, but gets better for bigger numbers of particles)
	printf("Cluster count: %d \n",cl_count);
	printf("Lc: %lf \n",Lc);
}
//-----------------Functions-----------------------

//--init_lat--
//initialises the Lattice with cells on it
//input:int N (number of sites),float phi (particle concentration)
//output: int array lattice (lattice with particles at timestep 0)
int * init_lat(int N,int M,float phi) 
{	
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
	double interval = 1/arraylength;	//separate the space 0-1 into N pieces with length interval
//	printf("interval: %f\n",interval);
	int i = 0;
	double rndnum;			
	long seed = time(NULL);
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
	long seed = time(NULL);
	//printf("seed = %ld\n",seed);
	double rndnum = ran3(&seed);
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
	long seed = time(NULL);
	int i = 0;
	int ind;
//	printf("M= %d\n",M);

//--this algorithm has to be optimized, there must be a shorter way--
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
			lattice[ind] = tumble(lattice[ind],alph);
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
/*
char * geturand()
{
	int len = 12;
	static char str[20];
	char file_name[25];
	FILE *fp;
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
		str[i] = fgetc(fp);
	//	printf("%c",str[i]);
	}
	fclose(fp);
	return str;	
}
*/

//--tumble()--
//input: dir (direction of the particle), alph (probability to tumble)
//output: direction of the particle (with probability alph changed)
int tumble(int dir,double alph)
{	
	int dir_ = dir;
	long seed = time(NULL);
	double rndnum = ran3(&seed);
/*
	//this is more complicated than needed here, but with two or three dimensions it is easier like that
	if(rndnum <= alph) 	//if random number is smaller than alph
	{
		while(dir == dir_)	// as long as the direction didnt change
		{
			dir_ = rnddirection(); 		//randomly choose new direction
		}
		return dir_;
	}
	else
	{
		return dir; //return old direction if random number is above probability
	}
*/
//not shure if the first way is right, so here an other one where tumbling can also not change the direction --> prob 1/2 to change it

	if(rndnum <= alph)
	{
		return rnddirection();
	}
	else
	{
		return dir;
	}
}


//--cluster_counter--
//input: the lattice
//output: number of clusters on lattice (a clusters (in 1D) consists particles that for shure won't change their position at the next time step if they don't tumble, the smallest cluster are two particles with direction towards each other)
int cluster_counter(int * lattice,int N)
{	
	int i, ii, iii, count=0;
//!!!!!!!!!!!!!-one big mistake I make here ist, that periodic boudary conditions are not considered-!!!!!!!!!!!!!!!!!!!!!!!!!
//this alters the count only by a small number so I dont care right now... I will change this though
	for(i=0;i<N;i++)
	{
		if(lattice[i] != 0)	//if there is a particle
		{
			ii = 1;
			while(lattice[i+ii] != 0) 	//look if there are particles to the right of it
			{
				ii++;
			}
			if(lattice[i+ii-1] == -lattice[i]) //if the last particle right of that particle has the opposite direction
			{
				i=i+ii;
				count++;
//				printf(" %d \n", count);	
				//!!!!!!it is possible to include measurement of cluster size here!!!!!!!
			}
			else	//if it has the same direction (wich is possible if the particle to its right tumbled and that particle didn't move at that timestep) search the next particle (with opposite direction to the first) left to this one.
			{
				iii=1;
				while((lattice[i+ii-iii] != -lattice[i]) && ((ii-iii) > 0)) //search the next left particle that has opposite direction to the first one
				{
					iii++;
				}
				if(ii-iii != 0) //if it is found add one to the counter
				{
					i=i+ii;
					count++;
//					printf(" %d \n", count);	
				}	
				else	//otherwise this is no cluster
				{
					i=i+ii;
				}			
			}
		}
	}
//	printf("%d \n",count);
	return count;
}

//--mean_dist--
//input: the lattice, NumOfTSteps (number of times the timestep function is executed, NumOfSweeps (number of times the whole lattice runs for NumOFTSteps times.  
//output: the mean distance for NumOfSweeps
//Important: this function only makes sense if there is only one Particle!
int mean_dist(int NumOfTSteps, int NumOfSweeps,int N, int M, float phi,double alph)
{	
	int * lattice;
	lattice = init_lat(N,M,phi);
	int i;
	for(i=0; i<NumOfTSteps; i++)
	{
		lattice = timestep(lattice,M,N,alph);
	}
}









