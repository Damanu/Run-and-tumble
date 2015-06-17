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
#include <assert.h>
#include "./../Lib/random.h"
#include "./../Lib/hk.h"
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
int mean_dist(int NumOfTSteps, int NumOfSweeps,int N, int M, float phi,double alph);
int find(int * lattice,int N);
struct particle * init_lat_2(int N,int M,float phi); 
struct particle * timestep_2(struct particle * lattice,int N,int M, double alph);
int * transform(struct particle * lattice,int M, int N);
double mean_dist_2(struct particle * lattice,int M,int N,double alph,int T);
struct particle * timestep_alone(struct particle * lattice, int N, int M, double alph);
int ** transform_2d(struct particle * lattice, int M,int m, int n);
//----------------Structures----------------------
struct particle
{
int ind;		//index on lattice
int dir;		//direction
int wallcount;		//counter how often it passed the end of lattice (right end is + left end is -)
int type;		//(optional) a type of particle (maybe not needed)
};

//-------------Global Variables--------------
static long int LATTICESIZE=10000000;
//--------------Main Program---------------------
int main() 
{
//-----------------Parameter input---------------------------
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
int i, ii;
struct particle * lattice;
int * r_lattice;
double m_d;
int ** matrix;
//----------------------------Manual input-------------------------------
/*	printf("Number of sites (N): ");
scanf("\n%d", &N);			//get number of sites
printf("Particle concentration (phi): ");
scanf("\n%f", &phi);			//get concentration
printf("Probability for tumbling (alpha): ");
scanf("\n%f", &alph);			//get tumbling probability
printf("Total Time for evolution (T): ");
scanf("\n%d",&tottime);
*/
//----------------------------------------------------------------

//----------mode----------
	int mode=2;
//-----------------------

	FILE *f;		//create file pointer		

	N = 60000;
	int m=1;
	alph =0.001;
	phi = 0.5;	
	int T = 10000;
	tottime=T;
	float M_ =(float)(N)*phi;	//M (number of Particles) --> if N*phi >= n.5 (with n natrual number) there is an error. This error is negligible for big N
	M=roundf(M_);
	switch (mode)
	{
	case 0:
		printf("no Mode");
		break;
	case 1:		//mean square distance calculation
		f = fopen("data_meandist_T.txt","w"); //open file stream
		fprintf(f,"meandist	T\n");
		for(i=0;i<1000;i++)
		{
			lattice = init_lat_2(N,i,phi);
			m_d=mean_dist_2(lattice,i,N,alph,T); 
	//		printf("mean dist: %lf \n",m_d);
			fprintf(f,"%lf	%d\n",m_d,i);
		}
		fclose(f);
		break;
	case 2:		//cluster-counting
		
		f = fopen("data_Fc_l.txt","w"); //open file stream
		fprintf(f,"l	Fc\n");
		double *Fc;
		Fc=(double *)calloc(M,sizeof(double));
		int *numofclusters;
		numofclusters=(int *)calloc(M,sizeof(int));
		int max_l=0;
		int iter=1000;
		lattice = init_lat_2(N,M,phi);
		for(i=0;i<T;i++)		//equilibration
		{
			lattice=timestep_2(lattice,N,M,alph);
			printf("i: %d\n",i);
		}
		for(ii=0;ii<iter;ii++)
		{
			printf("ii: %d\n",ii);
			for(i=0;i<100;i++)
			{
				lattice=timestep_2(lattice,N,M,alph);
			}
			matrix=transform_2d(lattice,M,m,N);
			int clusters = hoshen_kopelman(matrix,m,N);
			printf("clusters: %d\n",clusters);
			int ccount=0;
			int size=1;
			for(i=0;i<M;i++)
			{
				numofclusters[i]=0;
			}
			for(i=0;i<=N;i++)	//count the number of clusters of a specific size and write them into an array --> array[clustersize]=number of clusters with size clustersize
			{
				if(i==N)   //if the last site has been exeded
				{
					if(ccount==0)	//and there was no particle before
					{
						continue;	//do nothing
					}
					else 	//if there was one add one to the counter of its size
					{
						numofclusters[ccount]+=1;
						if(max_l<ccount)
						{
							max_l=ccount;
						}	 
						ccount=0;	//reset the counter
					}
					break;
				}
				if(matrix[m-1][i]==0) 	//if there is no particle
				{
					if(ccount==0)	//and there was no particle before
					{
						continue;	//do nothing
					}
					else 	//if there was one add one to the counter of its size
					{
						numofclusters[ccount]+=1;	 
						if(max_l<ccount)
						{
							max_l=ccount;
						}	 
						ccount=0;	//reset the counter
					}
				}
				else
				{
					ccount+=1;
				}
			}
			printf("max_l: %d\n",max_l);
			for(i=2;i<=max_l;i++)
			{	
				Fc[i]+=(double)numofclusters[i];
	//			Fc[i]=(double)numofclusters[i]/((double)M);
	//			fprintf(f,"%d	%lf\n",i,Fc[i]);
			}
		}
		for(i=2;i<=max_l;i++)
		{
			Fc[i]=Fc[i]/((double)M*(double)iter);
			fprintf(f,"%d	%lf\n",i,Fc[i]);
		}
		fclose(f);
	//	double l_m=(double)M/(double)clusters;	//mean cluster size (l_m)
	//	printf("mean cluster size: %lf\n",l_m);
		break;
	
	case 3: 	//calculating equilibration time (output average cluster size Lc over #timesteps T)
		
		f = fopen("data_Equilibrationtime.txt","w"); //open file stream
		fprintf(f,"T	Lc\n");
		lattice = init_lat_2(N,M,phi);		//initialize lattice
		T=1000000;
		int clusters;
		double Lc;
		double ti;
		int stepsize=200;
		for(i=stepsize;i<=T;i+=stepsize)
		{
			printf("i: %d\n",i);
			for(ii=0;ii<stepsize;ii++)
			{
				lattice=timestep_2(lattice,N,M,alph);		//make one timestep
			}
			matrix=transform_2d(lattice,M,m,N);		//create hk matrix
			clusters = hoshen_kopelman(matrix,m,N);		//calculate number of clusters
			for(ii=0;i<m;i++)
			{
				free(matrix[i]);
			}
			free(matrix);					//give space of matrix free (or I get a space problem)
			Lc=(double)M/(double)clusters;			//calculate mean cluster size
//			ti=(double)i/(double)T;
			fprintf(f,"%d	%lf\n",i,Lc);			//write data to file
		}
		fclose(f);	//close data stream
		break;
	}	
/*
	printf("M: %d\n",M);
	lattice = init_lat_2(N,M,phi);		//initalize lattice
	for(ii=0;ii<tottime;ii++)
	{
		r_lattice=transform(lattice,M,N);
		for(i=0;i<N;i++) 
		{
			//printf("%d",*(lattice+i));
			if(r_lattice[i]==0)
			{
				printf(" ");
			}
			else
			{		
				printf("|");
			}
		}
		printf("\n");
		lattice=timestep_2(lattice,N,M,alph);
	}
	printf("lattice index\n");
	for(i=0;i<M;i++)
	{
		printf("%d, ",lattice[i].ind);
	}
	printf("\n");
*/
//	int cl_count = cluster_counter(lattice,N);	
//	double Lc = (double)M/(double)cl_count;		//Lc is the average cluster size .... mistake: does not consider particles that are not in a cluster (all particles are in cluster ... wrong, but gets better for bigger numbers of particles)
//	printf("Cluster count: %d \n",cl_count);
//	printf("Lc: %lf \n",Lc);


}
//-----------------Functions-----------------------
//--init_lat--
//initialises the Lattice with cells on it
//input:int N (number of sites),float phi (particle concentration)
//output: int array lattice (lattice with particles at timestep 0)
int * init_lat(int N,int M,float phi) 
{	
	int * lattice;
	lattice=(int *)calloc(LATTICESIZE,sizeof(int));		//allocating 10000*sizeof(integer)bits space for the lattice array --> should be allocated dynamically, but didnt work till now	
	double interval=1/N;			//separate the space 0-1 into N pieces with length interval
	int i = 0;
	double rndnum;
	int len = N;
	int ind;
	for(i=0; i<N; i++)
	{
		lattice[i]=0;
	}
	for( i = 0 ; i < M ; i++ )  			//loop to find random indizes
	{
		ind=rand_index(len);
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
//--transform--
//input: struct particle * lattice
//output: int * r_lattice (transformed lattice to real lattice)
int * transform(struct particle * lattice,int M, int N)
{
	int i;
	static int  r_lattice[10000000];	//allocate memory for r_lattice
//	r_lattice=(int *)calloc(N,sizeof(int));
	for(i=0;i<N;i++)	//set all sites 0
	{
		r_lattice[i]=0;
	}
	for(i=0;i<M;i++)	//place particles on helper lattice
	{
		r_lattice[lattice[i].ind]=lattice[i].dir;
	}
	return r_lattice;
}
//--transform_2d--
//input: struct partilce * lattice, number of particles M, number of columns n, number of raws m)
//output: returns a ** consisting of integers -1, 1, 0
//important: after using this function the saving space for matrix AND all elements matrix[i] have to be freed by free(matrix) and free matrix[i] for all i
int ** transform_2d(struct particle * lattice, int M,int m, int n)
{
	int i,j;
	int ** matrix;
	matrix = (int **)calloc(m,sizeof(int *));
	for (i=0; i<m; i++)		//allocate memory 
	{
      		matrix[i] = (int *)calloc(n, sizeof(int));
		for(j=0;j<n;j++)
		{
			matrix[i][j]=0; //set all matrix elements 0	
		}
	}	
	for(i=0;i<M;i++)
	{
		int raw, column;
		raw=lattice[i].ind/n;
		column=lattice[i].ind-raw*n;
//		printf("raw: %d, column: %d \n",raw, column);
		matrix[raw][column]=lattice[i].dir;
	}
//	print_matrix(matrix,m,n);
	return matrix;
}
//--init_lat_2--
//initialises the Lattice with cells on it
//input:int N (number of sites),float phi (particle concentration)
//output: struct particle * lattice (lattice with particles at timestep 0)
struct particle * init_lat_2(int N,int M,float phi) 
{	
	struct particle * lattice;	//allocating 10000*sizeof(particle)bits space for the lattice array --> should be allocated dynamically, but didnt work till now	
	lattice = (struct particle *)calloc(M,sizeof(struct particle));
	int * r_lattice; 
	r_lattice=(int *)calloc(N,sizeof(int));		//helper lattice (real lattice)
	double interval=1/N;			//separate the space 0-1 into N pieces with length interval
	int i = 0;
	double rndnum;
	int ind;
	for(i=0; i<N; i++)
	{
		r_lattice[i]=0;
	}
	for(i=0;i<M;i++)  			//loop to find random indizes
	{
		ind=rand_index(N);			//find random index
		if(r_lattice[ind]==0) 			//if site is not occupied
		{	
			lattice[i].wallcount=0;
			lattice[i].type=0;
			lattice[i].ind=ind;		
			lattice[i].dir=rnddirection();	 
			r_lattice[ind]=lattice[i].dir;
		}
		else
		{
			i--;
		}
	}
	free(r_lattice);
	return lattice;
}
//--rand_index--
//returns a random index between 0 and arraylength
//input: double arraylength
//output: int random index
int rand_index(double arraylength) 
{
	double interval = 1.0/(arraylength-1);	//separate the space 0-1 into N pieces with length interval
//	printf("interval: %f\n",interval);
	int i = 0;
	double rndnum;			
	long seed = time(NULL);
	rndnum = ran3(&seed);			//pick a random number out of the spacing 0-1
	int index = roundf(rndnum/interval);	//evaluate the index of the lattice, this index is a randomly chosen one		
//	printf("rand_index, interval: %lf, rndnum: %lf, index: %d\n",interval,rndnum,index);
	return index;
}
//--rnddirection--
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
struct particle * timestep_alone(struct particle * lattice, int N, int M, double alph)
{	
	long seed = time(NULL);
	int i=0;
	int ind;
	int * r_lattice;
	r_lattice=transform(lattice,M,N);
	for(i=0;i<M;i++)
	{
		ind = rand_index(M);
		if(1) //if there is a particle 
		{
			lattice[ind].dir = tumble(lattice[ind].dir,alph);	//tumblin event
			if(lattice[ind].ind == N-1) //if upper periodic boundary 
			{
				if(lattice[ind].dir > 0) 		//find out the direction
				{
					lattice[ind].wallcount+=1;
					lattice[ind].ind=0; //move
					r_lattice[0]=lattice[ind].dir; //copy to helper lattice
					r_lattice[N-1]=0;	//clear old site
					continue;
				}
				else
				{
					lattice[ind].ind=N-2;		//move
					r_lattice[N-2]=lattice[ind].dir;	//copy to helper lattice
					r_lattice[N-1]=0;	//clear old site
					continue;
				}
			}
			else
			{
				if(lattice[ind].ind == 0)	//if lower periodic boundary
				{
					if(lattice[ind].dir > 0) 			//find out the direction
					{
						lattice[ind].ind=1;		//move
						r_lattice[1]=lattice[ind].dir;	//copy to helper lattice
						r_lattice[0]=0;	//clear old site
						continue;
					}
					else
					{	
						lattice[ind].ind=N-1;		//move
						lattice[ind].wallcount-=1;
						r_lattice[N-1]=lattice[ind].dir;	//copy to helper lattice
						r_lattice[0]=0;		//clear old site
						continue;
					}
				}
				else
				{
					if(lattice[ind].dir > 0) 			//find out the direction
					{
						lattice[ind].ind+=1;		//move
						r_lattice[lattice[ind].ind]=lattice[ind].dir;	//copy to helper lattice
						r_lattice[lattice[ind].ind-1]=0;	//clear old site
					}
					else	
					{
						lattice[ind].ind-=1;		//move
						r_lattice[lattice[ind].ind]=lattice[ind].dir;	//copy to helper lattice
						r_lattice[lattice[ind].ind+1]=0;	//clear old site
					}	
				}
			}
		}
		else
		{
			i--;
		}
	}
	return lattice;
}
//--timestep_2--
//input: struct particle * lattice (allocated in init)
//output: struct particle * lattice (evolution of lattice after one timestep)
struct particle * timestep_2(struct particle * lattice,int N,int M, double alph)
{
//	printf("timestep start\n");
	long seed = time(NULL);
	int i = 0;
	int ind;
	int * r_lattice;
	r_lattice=transform(lattice,M,N);
//--this algorithm has to be optimized, there must be a shorter way--
	for(i = 0; i < M; i++)
	{
		ind=rand_index(M);
//		printf("timestep .. lattice[ind]= %d\n",lattice[ind]);
//		printf("\n ind: %d\n\n",ind);
		if(r_lattice[lattice[ind].ind] != 0) //if there is a particle 
		{
			lattice[ind].dir = tumble(lattice[ind].dir,alph);	//tumbling event
			if(lattice[ind].ind == N-1) //if upper periodic boundary 
			{
				if(lattice[ind].dir > 0) 		//find out the direction
				{
					if(r_lattice[0] == 0) 	//find out if "the way if free"
					{
//						printf("upper -->\n");
//						printf("dir: %d\n",lattice[ind].dir);
//						printf("ind: %d\n",ind);
						lattice[ind].wallcount+=1;
						lattice[ind].ind=0; //move
						r_lattice[0]=lattice[ind].dir; //copy to helper lattice
						r_lattice[N-1]=0;	//clear old site
						continue;
					}
				}
				else
				{
					if(r_lattice[N-2] == 0)	//if the way is free
					{
//						printf("upper <--\n");
//						printf("dir: %d\n",lattice[ind].dir);
//						printf("ind: %d\n",ind);
						lattice[ind].ind=N-2;		//move
						r_lattice[N-2]=lattice[ind].dir;	//copy to helper lattice
						r_lattice[N-1]=0;	//clear old site
						continue;
					}	
				}
			}
			else
			{
				if(lattice[ind].ind == 0)	//if lower periodic boundary
				{
					if(lattice[ind].dir > 0) 			//find out the direction
					{
						if(r_lattice[1] == 0) 	//find out if "the way if free"
						{
//							printf("lower -->\n");
//							printf("dir: %d\n",lattice[ind].dir);
//							printf("ind: %d\n",ind);
							lattice[ind].ind=1;		//move
							r_lattice[1]=lattice[ind].dir;	//copy to helper lattice
							r_lattice[0]=0;	//clear old site
							continue;
						}
					}
					else
					{
						if(r_lattice[N-1] == 0) 		//if the way if free
						{
//							printf("lower <--\n");
//							printf("dir: %d\n",lattice[ind].dir);
//							printf("ind: %d\n",ind);
							lattice[ind].ind=N-1;		//move
							lattice[ind].wallcount-=1;
							r_lattice[N-1]=lattice[ind].dir;	//copy to helper lattice
							r_lattice[0]=0;		//clear old site
							continue;
						}	
					
					}
				}
				else
				{
					if(lattice[ind].dir > 0) 			//find out the direction
					{
						if(r_lattice[lattice[ind].ind+1] == 0) 	//find out if "the way if free"
						{
							lattice[ind].ind+=1;		//move
							r_lattice[lattice[ind].ind]=lattice[ind].dir;	//copy to helper lattice
							r_lattice[lattice[ind].ind-1]=0;	//clear old site
						}
					}
					else	
					{
						if(r_lattice[lattice[ind].ind-1] == 0)		//if the way is free
						{
							lattice[ind].ind-=1;		//move
							r_lattice[lattice[ind].ind]=lattice[ind].dir;	//copy to helper lattice
							r_lattice[lattice[ind].ind+1]=0;	//clear old site
						}
					}	
				}
			}
		}
		else
		{
			i--;
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
	if(alph==0)
	{
		return dir;	
	}
	if(rndnum < alph)
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
	int i,ii;
	int ind;
	double mean_distance=0;
	
	lattice = init_lat(N,M,phi);
//	printf("after init\n");
	//Declare nested function 
	//--find--
	//input: lattice
	//output: int ind (position of the particle)
	for(ii=0; ii<NumOfSweeps; ii++)		//run the whole counter NumOfSweeps times
	{
		
		ind = find(lattice,N); 
		for(i=0; i<NumOfTSteps; i++)	//let the particle evolve in time for one MCSweep 
		{
			lattice = timestep(lattice,N,M,alph);
		}
		mean_distance += (double)abs(ind - find(lattice,N));
	}
	mean_distance = mean_distance/(double)NumOfSweeps;	//calculate the mean distance over all runs
	printf("mean distance: %lf \n",mean_distance);
	return mean_distance;
}
//--mean_dist--
//input: lattice, number of particles M, Number of sites N, propability to tumble alph, number of timesteps to be made
//output: mean distance one particle moves

double mean_dist_2(struct particle * lattice,int M,int N,double alph,int T)
{
	int i;
	double dist=0;
	int * index;
	index=(int *)calloc(N,sizeof(int));
	for(i=0;i<M;i++)
	{
		index[i]=lattice[i].ind;
	}
	for(i=0;i<T;i++)
	{
		lattice=timestep_alone(lattice,N,M,alph);
	}
	for(i=0;i<M;i++)
	{
		dist += abs(abs(N*lattice[i].wallcount)-abs(lattice[i].ind-index[i]));
	}
	double m_d = dist/(double)M;
	free(index);
	return m_d;
}


int find(int * lattice,int N)
{
	int i;
//	printf("find\n");
	int store=-1;
	for(i=0; i<N; i++)		
	{
		if(lattice[i]!=0) 	
		{
			if(store==-1)
			{
//				printf("Find index: %d \n",i);
				store = i;
			}
			else
			{
				printf("meand_dist Error: there are two particles on the lattice, this function is only written for one!\n");
				exit(0);
			}
		}
		
	}
	return store;
}



