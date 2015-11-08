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
int rand_index(double arraylength);
int length(int * array);
int rnddirection(); 
//char * geturand();
int tumble(int dir,double alph);
int cluster_counter(int * lattice,int N);
int mean_dist(int NumOfTSteps, int NumOfSweeps,int N, int M, double phi,double alph);
int find(int * lattice,int N);
struct particle * init_lat_2(int N,int M,double phi,double alph,double tau); 
struct particle * timestep_2(struct particle * lattice,int N,int M, double alph);
double mean_dist_2(struct particle * lattice,int M,int N,double alph,int T);
struct particle * timestep_alone(struct particle * lattice, int N, int M, double alph);
int ** transform_2d(struct particle * lattice, int M,int m, int n);
int * transform(struct particle * lattice,int M, int N);
struct particle * timestep_2_2D(struct particle * lattice,int N,int m,int M, double alph);
int rnddirection_2D(); 
struct particle * init_lat_2_2D(int N,int M,double phi,double alph,double tau); 
int * get_index_2D(int ind,int N);
void waitFor (double secs);
double alph_change(double alph, double traptime,double tau);
int * cluster_counting(int ** matrix,int * numofclusters,int clusters,int N,int m);
int max_clustersize(int * numofclusters,int M);
//----------------Structures----------------------
struct particle
{
int ind;		//index on lattice
int dir;		//direction
int wallcount;		//counter how often it passed the end of lattice (right end is + left end is -)
int type;		//(optional) a type of particle (maybe not needed)
int wallcount_1;	//wallcount for 2D in 1 -1 direction
int wallcount_2;	//wallcount for 2D in 2 -2 direction
double alph;		//probabilitiy of tumbling
double traptime;	//time the particle was not able to move
double tau;		//rat with which the particle changes alpha (if trapped)
};

//-------------Global Variables--------------
int ALPH_CHANGE_TYPE=1;
//--------------Main Program---------------------
int main(int argc, char *argv[]) 
{
printf("\n--------------------\n");
printf("2D-Testing started\n");
printf("--------------------\n\n");
//----------------Argument info-------------------------
if(argv[1] == 0)
{
	printf("--------Help-------\n\
 Argument Info:\n\
	1. Mode\n\
		0: Testmode\n\
		1: mean square distance 1D\n\
		2:clustersizedistribution 1D\n\
		3: equilibration time 1D\n\
		4: visual 2D lattice output\n\
		5: clustersizedistribution 2D\n\
		6: data collapse plot 1D Lc over lc(other parameters are irrelevant, set 0!)\n\
		7:equilibration time 2D\n\
	2. Clustersize N (NxN in 2D)\n\
	3. alph\n\
	4. phi\n\
	5. initial time T\n\
	6. aditional name ([Name]_[aditional Name])\n\
	7. Tau (rate of changing alpha when trapped)\n");
	exit(0);
}

//-----------------Parameter input---------------------------
srand(time(0));
//double r = rand()/(double)RAND_MAX;
long t=time(0);
double r=ran3(&t);
double alph, phi;			//alph: propability for tumbling event; phi: particle concentration
int M,M_2D,N,tottime;			//M: total number of particles; N: total number of sites (or length of lattice array)
char word;
int i, ii;
struct particle * lattice;
int * r_lattice;
double m_d;
int ** matrix;
double Lc;
double mean_alph;

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
	int mode=atoi(argv[1]);
//-----------------------

	FILE *f;		//create file pointer		

	double tau=atof(argv[7]);
	N = atoi(argv[2]);
	int m=N;
	alph =atof(argv[3]);	//convert to double	
	phi = atof(argv[4]);	//convert to double
	int T = atoi(argv[5]);	//convert to int
	char output[100];
	sprintf(output,"Data_N=%d_alph=%1.3lf_phi=%1.3lf_T=%d_Tau=%1.3lf_%s",N,alph,phi,T,tau,argv[6]);
//	printf("%s\n",output);
//	printf("mode: %d\nN=%d\nalph=%.3f\nphi=%.3f\n",mode,N,alph,phi);
	tottime=T;
	float M_ =(float)(N)*phi;	//M (number of Particles) --> if N*phi >= n.5 (with n natrual number) there is an error. This error is negligible for big N
	float M_2D_=(float)(N)*(float)(N)*phi;
	M=roundf(M_);
	M_2D=roundf(M_2D_);
//	printf("before switch");
	switch (mode)
	{
	case 0:
		printf("Testmode");
		lattice = init_lat_2(100,100,0.5,0.5,100);
		int t=0;
		for(t=0;t<100;t++)
		{	
			printf("alph: %lf\n",alph_change(0.5,t,100));
		}
//		f = fopen("Empty_%d+%d_alph=%f_phi=%f","w",N,N,alph,phi);
//		fclose(f);
		break;
	case 1:		//mean square distance calculation
		printf("mode 1\n");
		f = fopen("data_meandist_T.txt","w"); //open file stream
		fprintf(f,"meandist	T\n");
		for(i=0;i<1;i++)
		{
	//		printf("in loopi\n");
			lattice = init_lat_2(N,M,phi,alph,tau);
	//		printf("after init\n");
			m_d=mean_dist_2(lattice,M,N,alph,T);
			 
			printf("mean dist: %lf \n",m_d);
			fprintf(f,"%lf	%d\n",m_d,M);
		}
		fclose(f);
		break;
	case 2:		//cluster-counting
		
		f = fopen(output,"w"); //open file stream
		fprintf(f,"l	Fc\n");
		double * Fc;
		Fc=(double *)calloc(M,sizeof(double));
			
		int * numofclusters = (int *)calloc(M,sizeof(int));	//the index stands for the clustersize and the number of Numofclusters[index] stands for the number of clusters with that size
		int max_l;
		int iter=1;
	
		for(ii=0;ii<iter;ii++)
		{
			lattice = init_lat_2(N,M,phi,alph,tau);
			for(i=0;i<T;i++)		//equilibration
			{
				lattice=timestep_2(lattice,N,M,alph);
		//		printf("i: %d\n",i);
			}	//	printf("ii: %d\n",ii);
	//		for(i=0;i<100;i++)
	//		{
	//			lattice=timestep_2(lattice,N,M,alph);
	//		}
			matrix=transform_2d(lattice,M,1,N);
			int clusters = hoshen_kopelman(matrix,1,N);
		//	printf("clusters: %d\n",clusters);
			int ccount=0;
			int size=1;
			numofclusters=cluster_counting(matrix,numofclusters,clusters,N,1);
			max_l=max_clustersize(numofclusters,M);
			

		}

		printf("max_l: %d\n",max_l);
		for(i=2;i<max_l;i++)
		{
			Fc[i]=(double)numofclusters[i]/((double)M*(double)iter);
			if (Fc[i]!=0)
			{
				fprintf(f,"%d	%1.10lf\n",i,Fc[i]);
			}
		}
		fclose(f);
	//	double l_m=(double)M/(double)clusters;	//mean cluster size (l_m)
	//	printf("mean cluster size: %lf\n",l_m);
		break;
	
	case 3: 	//calculating equilibration time (output average cluster size Lc over #timesteps T)
		
		//f = fopen("data_Equilibrationtime.txt","w"); //open file stream
		f = fopen(output,"w"); //open file stream
		fprintf(f,"T	Lc	<alph>\n");
		lattice = init_lat_2(N,M,phi,alph,tau);		//initialize lattice
		m=1;
		int clusters;
		double ti;
		int stepsize=10;
		for(i=stepsize;i<=T;i+=stepsize)
//>>>>>>> 6b0390817e69264e7757c8e463f1fcde443336d4
		{
//			printf("i: %d\n",i);
			for(ii=0;ii<stepsize;ii++)
			{
				lattice=timestep_2(lattice,N,M,alph);		//make one timestep
			}
			matrix=transform_2d(lattice,M,1,N);		//create hk matrix
			clusters = hoshen_kopelman(matrix,m,N);		//calculate number of clusters
		//	printf("clusters: %d\n", clusters);
			
			//--measuere mean tumbling rate---
			mean_alph=0;
			for(ii=0;ii<M;ii++)
			{
				mean_alph+=lattice[ii].alph;
			}
			mean_alph=mean_alph/(double)(M);
			//-----------------------

			for(ii=0;ii<m;ii++)
			{
				free(matrix[ii]);
			}
			free(matrix);					//give space of matrix free (or I get a space problem)
			Lc=(double)M/(double)clusters;			//calculate mean cluster size
		//	printf("Lc: %lf\n",Lc);
//			ti=(double)i/(double)T;
			fprintf(f,"%d	%lf	%lf\n",i,Lc,mean_alph);			//write data to file
		}
		fclose(f);	//close data stream
		break;
	case 4: 	//2D visualisation output is position of particles (snapshot)
		f = fopen(output,"w"); //open file stream
		fprintf(f,"x	y\n"); 		//A is the clustersize and Fc the distribution of it
		M_ =(float)(N)*(float)(N)*phi;	//M (number of Particles) --> if N*phi >= n.5 (with n natrual number) there is an error. This error is negligible for big N
		M=roundf(M_);
		lattice=init_lat_2_2D(N*N,M,phi,alph,tau);
		int ** matrix;
//		matrix=transform_2d(lattice,M,N,N);
//		print_matrix(matrix,N,N);
/*		double bla=0;
 calculating mean value of rnddirection_2D --> should be near 0
		for(i=0;i<100000000;i++)
		{
			bla+=rnddirection_2D();
		}
		printf("mean value: %lf",bla/100000000);
*/
		for(i=0;i<T;i++)
		{
//			waitFor(0.1);
			lattice=timestep_2_2D(lattice,N,N,M,alph);
//			system("clear");
//			matrix=transform_2d(lattice,M,N,N);
//			print_matrix(matrix,N,N);
		}
		matrix=transform_2d(lattice,M,N,N);
		clusters = hoshen_kopelman(matrix,N,N);
//		printf("----------------------------\n");
//		print_matrix(matrix,N,N);
		printf("clusters: %d\n",clusters);
		printf("alph: %lf\n",lattice[0].alph);
		
		for(i=0;i<M;i++)
		{
			int x=lattice[i].ind/N;		//calculate raw index in 2D matrix
			int y=lattice[i].ind-x*N;	//calculate column index in 2D matrix
			fprintf(f,"%d	%d\n",x,y);
			printf("%d	%d\n",x,y);
		}
		fclose(f);
		break;
	case 5:		//clustercounting to get clustersize distribution
		printf("mode: %d\n",mode);
		M_ =(float)(N)*(float)(N)*phi;	//M (number of Particles) --> if N*phi >= n.5 (with n natrual number) there is an error. This error is negligible for big N
		M=roundf(M_);
		int ccount,clsize=0;	//ccount is the number of the cluster looked at, clsize is the size of that cluster
		int * numofclusters_2D = (int *)calloc(M,sizeof(int));	//the index stands for the clustersize and the number of Numofclusters[index] stands for the number of clusters with that size
		int n,numofmeasure=1000;		//n ,number of measurements
//		printf("number of iterations: %d",numofmeasure);
		

	
		

		int T_step=100;
		f = fopen(output,"w"); //open file stream
		fprintf(f,"A	Fc\n"); 		//A is the clustersize and Fc the distribution of it

		int max_A;
		
		for(n=0;n<numofmeasure;n++)
		{
//			printf("n: %d",n);
			lattice=init_lat_2_2D(N*N,M,phi,alph,tau);	//initialize 2D lattice
			for(i=0;i<T;i++)	//do T timesteps to get equilibrium
			{
				lattice=timestep_2_2D(lattice,N,N,M,alph);
	//			printf("i: %d\n",i);
			}
		//	for(i=0;i<T_step;i++)	//do T_step timesteps
		//	{
		//		lattice=timestep_2_2D(lattice,N,N,M,alph);
		//	}
			matrix=transform_2d(lattice,M,N,N);
//			print_matrix(matrix,N,N);
			clusters = hoshen_kopelman(matrix,N,N);
			numofclusters_2D=cluster_counting(matrix,numofclusters_2D,clusters,N,N);
			max_A=max_clustersize(numofclusters,M);
		
		}

		double Fc_2D;
		for(i=1;i<=max_A;i++)
		{
			Fc_2D=(double)numofclusters_2D[i]/(((double)M)*(double)numofmeasure);		//normalisation of Fc
			if(Fc_2D!=0)
			{
				fprintf(f,"%d	%1.10lf\n",i,Fc_2D);
			}
		}
		fclose(f);
		break;
	
	case 6:		//data collapse average cluster size Lc as function of lengthscale lc 1D
		
		f = fopen(output,"w"); //open file stream
		fprintf(f,"lc	Lc	phi	alph\n");
		double lc;
		int j;
		int jj;
	//	for(i=1;i<=101;i+=10)
	//	{
//			printf("i: %d\n",i);
//			phi = (float)i/100.0;
			for(ii=1;ii<=1001*4;ii+=100)
			{
//				printf("ii: %d\n",ii);
		//		alph=-log((float)ii/1000.0)*0.999/log(0.001)+1;
		//		alph=-exp((float)ii/1000.0)*0.999/(exp(0.001)-exp(1))+1+exp(1)*0.999/(exp(0.001)-exp(1));
	//			alph=-exp((float)ii/1001.0)*0.999/(exp(0.001)-exp(1))+1+exp(1)*0.999/(exp(0.001)-exp(1));
				double k=(1.0-(0.001-pow(10.0,-0.999))/(1-pow(10.0,-0.999)))/10;
				double d=(0.001-pow(10.0,-0.999))/(1-pow(10.0,-0.999));
				alph=pow(10.0,(double)ii/1001.0/4)*k+d; 
				lc=sqrt(2*phi/((1-phi)*alph));
				float M_ =(float)(N)*phi;	//M (number of Particles) --> if N*phi >= n.5 (with n natrual number) there is an error. This error is negligible for big N
				M=roundf(M_);
				lattice = init_lat_2(N,M,phi,alph,tau);				
				for(j=0;j<T;j++)		//equilibration
				{
					lattice=timestep_2(lattice,N,M,alph);
			//		printf("i: %d\n",i);
				}
				matrix=transform_2d(lattice,M,1,N);
				int clusters = hoshen_kopelman(matrix,1,N);


/*				int count;
				int M_in_c=M;		//# particles in clusters
				int real_clusters=clusters;
				for(j=1;j<=clusters;j++)
				{
					count=0;
					for(jj=0;jj<N;jj++)
					{
						for(i=0;i<N;i++)
						{
							if(matrix[0][i]==j)
							{
								count+=1;
							}
						}
						if(count==1)
						{
							M_in_c-=1;
							real_clusters-=1;
						}
						else
						{
							
						}
					}
				}
				if(real_clusters!=0)
*/
				{
					Lc=(double)M/(double)clusters;
//					printf("rc: %d\n",clusters);
//					printf("c: %d\n",clusters);	
//					printf("mc: %d\n",M_in_c);
					fprintf(f,"%lf	%lf	%lf	%lf\n",lc,Lc,phi,alph);
				}
				free(lattice);
				free(matrix);
			}
	//	}
		break;	
	case 7: 	//calculating equilibration time 2D (output average cluster size Lc over #timesteps T)
		printf("mode 7\n");	
		f = fopen(output,"w"); //open file stream
		fprintf(f,"T	Lc	<alph>\n");
		lattice = init_lat_2_2D(N*N,M_2D,phi,alph,tau);		//initialize lattice
		int clusters_2;
		double ti_2;
		int stepsize_2=1;
		for(i=stepsize_2;i<=T;i+=stepsize_2)
		{
		//	printf("i: %d\n",i);
			for(ii=0;ii<stepsize_2;ii++)
			{
				lattice=timestep_2_2D(lattice,N,N,M_2D,alph);		//make one timestep
			}
			matrix=transform_2d(lattice,M_2D,N,N);		//create hk matrix
			clusters_2 = hoshen_kopelman(matrix,N,N);		//calculate number of clusters
		//	printf("clusters: %d\n",clusters_2);

			//--measuere mean tumbling rate---
			mean_alph=0;
			for(ii=0;ii<M_2D;ii++)
			{
				mean_alph+=lattice[ii].alph;
			}
			mean_alph=mean_alph/(double)(M_2D);
			//-----------------------

			for(ii=0;ii<m;ii++)
			{
				free(matrix[ii]);
			}
			free(matrix);					//give space of matrix free (or I get a space problem)
			Lc=(double)M_2D/(double)clusters_2;			//calculate mean cluster size
//			ti=(double)i/(double)T;
			fprintf(f,"%d	%lf	%lf\n",i,Lc,mean_alph);			//write data to file
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
//--transform_2d--
//input: struct partilce * lattice, number of particles M, number of columns n, number of raws m)
//output: returns a ** consisting of integers -1, 1, 0
//important: after using this function the saving space for matrix AND all elements matrix[i] have to be freed by free(matrix) and free matrix[i] for all i
int ** transform_2d(struct particle * lattice, int M,int m, int n)
{
	int i,j;
	int ** matrix;
	matrix = (int **)calloc(m,sizeof(int *));	//allocate memory
//	static int matrix[10000][10000];
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
		raw=lattice[i].ind/n;		//calculate raw index in 2D matrix
		column=lattice[i].ind-raw*n;	//calculate column index in 2D matrix
//		printf("raw: %d, column: %d \n",raw, column);
		matrix[raw][column]=lattice[i].dir;
	}
//	print_matrix(matrix,m,n);
	return matrix;
}
//--transform--
//input: struct particle * lattice
//output: int * r_lattice (transformed lattice to real lattice)
int * transform(struct particle * lattice,int M, int N)
{
	int i;
//	static int  r_lattice[10000000];	//allocate memory for r_lattice
	int *r_lattice=(int *)calloc(N,sizeof(int));
//	for(i=0;i<N;i++)	//set all sites 0
//	{
//		r_lattice[i]=0;
//	}
	for(i=0;i<M;i++)	//place particles on helper lattice
	{
		r_lattice[lattice[i].ind]=lattice[i].dir;
	}
	return r_lattice;
}
//--init_lat_2--
//initialises the Lattice with cells on it
//input:int N (number of sites),float phi (particle concentration)
//output: struct particle * lattice (lattice with particles at timestep 0)
struct particle * init_lat_2(int N,int M,double phi,double alph,double tau) 
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
			lattice[i].tau=tau;
			lattice[i].traptime=0;
			lattice[i].alph=alph;	
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
		if(r_lattice[lattice[ind].ind] != 0) //if there is a particle 
		{
			lattice[ind].alph=alph_change(lattice[ind].alph,lattice[ind].traptime,lattice[ind].tau);	//change alpha with rate tau
			lattice[ind].dir = tumble(lattice[ind].dir,lattice[ind].alph);	//tumbling event
			if(lattice[ind].ind == N-1) //if upper periodic boundary 
			{
				if(lattice[ind].dir > 0) 		//find out the direction
				{
					if(r_lattice[0] == 0) 	//find out if "the way if free"
					{
						lattice[ind].traptime=0;	//add one timeunit to the trap time
						lattice[ind].wallcount+=1;
						lattice[ind].ind=0; //move
						r_lattice[0]=lattice[ind].dir; //copy to helper lattice
						r_lattice[N-1]=0;	//clear old site
						continue;
					}

					else		//if the way is not free
					{
						lattice[ind].traptime+=1;	//add one timeunit to the trap time
						
					}
				}
				else
				{
					if(r_lattice[N-2] == 0)	//if the way is free
					{
						lattice[ind].traptime=0;	//add one timeunit to the trap time
						lattice[ind].ind=N-2;		//move
						r_lattice[N-2]=lattice[ind].dir;	//copy to helper lattice
						r_lattice[N-1]=0;	//clear old site
						continue;
					}
					else		//if the way is not free
					{
						lattice[ind].traptime+=1;	//add one timeunit to the trap time
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
							lattice[ind].traptime=0;	//add one timeunit to the trap time
							lattice[ind].ind=1;		//move
							r_lattice[1]=lattice[ind].dir;	//copy to helper lattice
							r_lattice[0]=0;	//clear old site
							continue;
						}
						else		//if the way is not free
						{
							lattice[ind].traptime+=1;	//add one timeunit to the trap time
						}
				}
					else
					{
						if(r_lattice[N-1] == 0) 		//if the way if free
						{
							lattice[ind].traptime=0;	//add one timeunit to the trap time
							lattice[ind].ind=N-1;		//move
							lattice[ind].wallcount-=1;
							r_lattice[N-1]=lattice[ind].dir;	//copy to helper lattice
							r_lattice[0]=0;		//clear old site
							continue;
						}	
						else		//if the way is not free
						{
							lattice[ind].traptime+=1;	//add one timeunit to the trap time
						}
	
					
					}
				}
				else
				{
					if(lattice[ind].dir > 0) 			//find out the direction
					{
						if(r_lattice[lattice[ind].ind+1] == 0) 	//find out if "the way if free"
						{
							lattice[ind].traptime=0;	//add one timeunit to the trap time
							lattice[ind].ind+=1;		//move
							r_lattice[lattice[ind].ind]=lattice[ind].dir;	//copy to helper lattice
							r_lattice[lattice[ind].ind-1]=0;	//clear old site
						}	
						else		//if the way is not free
						{
							lattice[ind].traptime+=1;	//add one timeunit to the trap time
						}

					}
					else	
					{
						if(r_lattice[lattice[ind].ind-1] == 0)		//if the way is free
						{
							lattice[ind].traptime=0;	//add one timeunit to the trap time
							lattice[ind].ind-=1;		//move
							r_lattice[lattice[ind].ind]=lattice[ind].dir;	//copy to helper lattice
							r_lattice[lattice[ind].ind+1]=0;	//clear old site
						}	
						else		//if the way is not free
						{
							lattice[ind].traptime+=1;	//add one timeunit to the trap time
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
	free(r_lattice);
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
//input: lattice, number of particles M, Number of sites N, propability to tumble alph, number of timesteps to be made
//output: mean distance one particle moves

double mean_dist_2(struct particle * lattice,int M,int N,double alph,int T)
{
//	printf("in md2");
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
		dist +=pow( (abs(N*lattice[i].wallcount)-abs(lattice[i].ind-index[i])),2);
//		dist += (lattice[i].ind-index[i])^2;
	}
	double m_d = dist/(double)M;
	free(index);
//	printf("m_d: %lf",m_d);
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
void waitFor (double secs)
{
    double retTime = time(0) + secs;     // Get finishing time.
    while (time(0) < retTime);    // Loop until it arrives.
}

//--------------------2D Expansion-------------------

//--rnddirection:--
//randomly chooses -1,1,-2 or 2
//input none
//output: -1,1,-2 or 2 as int
int rnddirection_2D() 
{		
//	printf("rnddirrection_2D\n");
	long seed = time(NULL);
	//printf("seed = %ld\n",seed);
	double rndnum = ran3(&seed);
	int sw;
	if(rndnum<0.25)		sw=1;
	else 	if (rndnum<0.5)		sw=2;
		else	if (rndnum<0.75)	sw=3;
			else	if(rndnum<1)		sw=4;
	switch (sw) 
	{
		case 1:
			return -1;
			break;
		case 2:
			return 1;
			break;
		case 3:
			return -2;
			break;
		case 4:
			return 2;
			break;
		
	}
}
//--init_lat_2_2D--
//initialises the Lattice with cells on it
//input:int N (number of sites),float phi (particle concentration)
//output: struct particle * lattice (lattice with particles at timestep 0)
struct particle * init_lat_2_2D(int N,int M,double phi,double alph,double tau) 
{	
//	printf("init lat\n");
	struct particle * lattice;	//allocating 10000*sizeof(particle)bits space for the lattice array --> should be allocated dynamically, but didnt work till now	
	lattice = (struct particle *)calloc(M,sizeof(struct particle));
	int * r_lattice; 
	r_lattice=(int *)calloc(N,sizeof(int));		//helper lattice (real lattice)
	int i = 0;
	int ind;


	for(i=0; i<N; i++)	//set all elements of real matrix 0
	{
		r_lattice[i]=0;
	}


	for(i=0;i<M;i++)  			//loop to find random indizes
	{
		ind=rand_index(N);			//find random index
		if(r_lattice[ind]==0) 			//if site is not occupied
		{	
			lattice[i].tau=tau;
			lattice[i].traptime=0;
			lattice[i].alph=alph;	
			lattice[i].wallcount=0;
			lattice[i].type=0;
			lattice[i].ind=ind;		
			lattice[i].dir=rnddirection_2D();	 
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

//--timestep_2_2D--
//input: struct particle * lattice (allocated in init)
//output: struct particle * lattice (evolution of lattice after one timestep)
struct particle * timestep_2_2D(struct particle * lattice,int N,int m,int M, double alph)
{
//	printf("timestep start\n");
	long seed = time(NULL);
	int i = 0;
	int ind;
	int ** r_lattice;
	int *ind2D;
	r_lattice=transform_2d(lattice,M,N,N);
	for(i = 0; i < M; i++)
	{
		ind=rand_index(M);	//generate random index
		ind2D=get_index_2D(lattice[ind].ind,N); //get index in 2D lattice ind2D[0]=row ind2D[1]=column
//		printf("timestep .. lattice[ind]= %d\n",lattice[ind]);
//		printf("\n ind: %d\n\n",ind);
		if(r_lattice[ind2D[0],ind2D[1]] != 0) //if there is a particle 
		{
			lattice[ind].alph=alph_change(lattice[ind].alph,lattice[ind].traptime,lattice[ind].tau);	//change alpha with rate tau
			lattice[ind].dir = tumble_2D(lattice[ind].dir,lattice[ind].alph);	//tumbling event
			switch (lattice[ind].dir)
			{
				case 1:		//moving right	
//					printf("case 1\n");
//					printf("z: %d,sp:%d N-1=%d\n",ind2D[0],ind2D[1],N-1);
					if(ind2D[1]==N-1)	//if particle is on right boundary
					{
						if(r_lattice[ind2D[0]][0]==0)	//if the way is free
						{
							lattice[ind].traptime=0;	//add one timeunit to the trap time
							lattice[ind].wallcount_1+=1;
							lattice[ind].ind=N*ind2D[0]+ind2D[1]-N+1; //move
							r_lattice[ind2D[0]][0]=lattice[ind].dir; //copy to helper lattice
							r_lattice[ind2D[0]][N-1]=0;	//clear old site
							continue;
						}
						
						else		//if the way is not free
						{
							lattice[ind].traptime+=1;	//add one timeunit to the trap time
							
						}
	
					}
					else // if not
					{
//						printf("else 1\n");
						if(r_lattice[ind2D[0]][ind2D[1]+1]==0)	//if the way is free
						{	
							lattice[ind].traptime=0;	//add one timeunit to the trap time
							lattice[ind].ind=(ind2D[0])*N+ind2D[1]+1; //move (lattice.ind=raw*N+column)
							r_lattice[ind2D[0]][ind2D[1]+1]=lattice[ind].dir; //copy to helper lattice
							r_lattice[ind2D[0]][ind2D[1]]=0;	//clear old site
							continue;
						}	
						else		//if the way is not free
						{
							lattice[ind].traptime+=1;	//add one timeunit to the trap time
							
						}
					}
					break;
				case -1:	//moving left
//					printf("case -1\n");
//					printf("index: %d,%d\n",ind2D[0],ind2D[1]);
					if(ind2D[1]==0)	//if particle is on left boundary
					{
						if(r_lattice[ind2D[0]][N-1]==0)	//if the way is free
						{
							lattice[ind].traptime=0;	//add one timeunit to the trap time
							lattice[ind].wallcount_1-=1;
							lattice[ind].ind=N*ind2D[0]+ind2D[1]+N-1; //move
							r_lattice[ind2D[0]][N-1]=lattice[ind].dir; //copy to helper lattice
							r_lattice[ind2D[0]][0]=0;	//clear old site
							continue;
						}	
						else		//if the way is not free
						{
							lattice[ind].traptime+=1;	//add one timeunit to the trap time
							
						}
					}
					else // if not
					{
						if(r_lattice[ind2D[0]][ind2D[1]-1]==0)	//if the way is free
						{	
							lattice[ind].traptime=0;	//add one timeunit to the trap time
							lattice[ind].ind=(ind2D[0])*N+ind2D[1]-1; //move (lattice.ind=raw*N+column)
							r_lattice[ind2D[0]][ind2D[1]-1]=lattice[ind].dir; //copy to helper lattice
							r_lattice[ind2D[0]][ind2D[1]]=0;	//clear old site
							continue;
						}	
						else		//if the way is not free
						{
							lattice[ind].traptime+=1;	//add one timeunit to the trap time
							
						}
					}
					break;
				case 2:		//moving down
//					printf("case 2\n");
					if(ind2D[0]==m-1)	//if particle is on upper boundary
					{
						if(r_lattice[0][ind2D[1]]==0)	//if the way is free
						{
							lattice[ind].traptime=0;	//add one timeunit to the trap time
							lattice[ind].wallcount_2+=1;
							lattice[ind].ind=/*(ind2D[0]+1)*N*/0*N+ind2D[1]; //move (lattice.ind=raw*N+column)
							r_lattice[0][ind2D[1]]=lattice[ind].dir; //copy to helper lattice
							r_lattice[m-1][ind2D[1]]=0;	//clear old site
							continue;
						}	
						else		//if the way is not free
						{
							lattice[ind].traptime+=1;	//add one timeunit to the trap time
							
						}
					}
					else // if not
					{
						if(r_lattice[ind2D[0]+1][ind2D[1]]==0)	//if the way is free
						{	
							lattice[ind].traptime=0;	//add one timeunit to the trap time
							lattice[ind].ind=(ind2D[0]+1)*N+ind2D[1]; //move (lattice.ind=raw*N+column)
							r_lattice[ind2D[0]+1][ind2D[1]]=lattice[ind].dir; //copy to helper lattice
							r_lattice[ind2D[0]][ind2D[1]]=0;	//clear old site
							continue;
						}	
						else		//if the way is not free
						{
							lattice[ind].traptime+=1;	//add one timeunit to the trap time
							
						}
					}
					break;
				case -2:	//moving up
//					printf("case -2\n");
					if(ind2D[0]==0)	//if particle is on lower boundary
					{
						if(r_lattice[N-1][ind2D[1]]==0)	//if the way is free
						{
							lattice[ind].traptime=0;	//add one timeunit to the trap time
							lattice[ind].wallcount_2-=1;
							lattice[ind].ind=/*(ind2D[0]+1)*N*/(N-1)*N+ind2D[1]; //move (lattice.ind=raw*N+column)
							r_lattice[N-1][ind2D[1]]=lattice[ind].dir; //copy to helper lattice
							r_lattice[0][ind2D[1]]=0;	//clear old site
							continue;
						}	
						else		//if the way is not free
						{
							lattice[ind].traptime+=1;	//add one timeunit to the trap time
							
						}
					}
					else // if not
					{
						if(r_lattice[ind2D[0]-1][ind2D[1]]==0)	//if the way is free
						{	
							lattice[ind].traptime=0;	//add one timeunit to the trap time
							lattice[ind].ind=(ind2D[0]-1)*N+ind2D[1]; //move (lattice.ind=raw*N+column)
							r_lattice[ind2D[0]-1][ind2D[1]]=lattice[ind].dir; //copy to helper lattice
							r_lattice[ind2D[0]][ind2D[1]]=0;	//clear old site
							continue;
						}	
						else		//if the way is not free
						{
							lattice[ind].traptime+=1;	//add one timeunit to the trap time
							
						}
					}
					break;
				
			}
		}
	}
	//----here was a memory leak, solved with this loop, forgot the inner pointer of matrix (or r_lattice)
	for(i=0;i<N;i++)
	{
		free(r_lattice[i]);
	}
	free(r_lattice);
	return lattice;
}

int * get_index_2D(int ind,int N)
{
//	printf("get_index_2D\n");	
	static int index[2];
	index[0]=0;
	index[1]=0;
	index [0]=ind/N;
	index[1]=ind-index[0]*N;
	return index;
}

//--tumble_2D()--
//input: dir (direction of the particle), alph (probability to tumble)
//output: direction of the particle (with probability alph changed)
int tumble_2D(int dir,double alph)
{
//	printf("tumble_2D\n");	
	int dir_ = dir;
	long seed = time(NULL);
	double rndnum = ran3(&seed);
	if(alph==0)
	{
//		printf("tumbe 1\n");
		return dir;	
	}
	if(rndnum < alph)
	{
//		printf("tumbe 2\n");
		return rnddirection_2D();
	}
	else
	{
		return dir;
	}
}
//--alph_change()--
//input: typ (type of changing), alph, traptime (time a particle cant move), tau (rate of poisson distribution)
//output: alph that may has changed according to the type of changing
double alph_change(double alph, double traptime,double tau)
{
	long seed = time(NULL);
	double rndnum = ran3(&seed);
	switch(ALPH_CHANGE_TYPE)
	{
		case 1: //random alpha change with prob poisson
			if (rndnum > exp(-traptime/tau))	//if random number is bigger than poisson dist at t = traptime with rate 1/tau
			{
//				return alph;
				return ran3(&seed);	//return a random alpha
			}
			else
			{
				return alph;		//return alpha like it was
			}
			break;
	}
}
//-----cluster counting--------
//input: matrix, clusters (number of clusters, N, m
//output: pointer numofclusters (numofclusters[1] is the number of clusters with size 1)
int * cluster_counting(int ** matrix,int * numofclusters,int clusters,int N,int m)
{
	int ccount;
	int i=0;
	int j=0;
	int clsize;
	for(ccount=1;ccount<=clusters;ccount++)
	{
		clsize=0;
		for(i=0;i<N;i++)	//search in whole matrix
		{
			for(j=0;j<m; j++)
			{
				if(matrix[j][i]==ccount)
				{
					clsize+=1;
				}
			}
	
		}
		numofclusters[clsize]+=1;
	}
	return numofclusters;
}
//----------maximal clustesize---------------
//input: pointer numberofclusters, M (number of particles)
//output: max_A (maximal clustersize) 
int max_clustersize(int * numofclusters,int M)
{
	int max_A; 	//maximal cluster size
	int i=0;
	for(i=1;i<M;i++)
	{
		if(numofclusters[i]!=0)
		{
			max_A=i;
		}
	}
	return max_A;
}
