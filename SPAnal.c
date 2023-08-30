#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <string>
#include <unistd.h>
#include <fftw3.h>
#include <time.h>


double randn_notrig();
int main(int argc, char **argv) 
{
  
	double PI=acos(-1.0);
	int t,k,j,i;

	int N=20000;
	double dt=0.0005;
	int t0=1;

	int N1=10000;
	double *x;

	int Overlap=1000;

	int KK;

	KK=(N-N1)/Overlap+1;

	double *f1, *f2;

	x = (double*) malloc(sizeof(double)*N);
	f1 = (double*) malloc(sizeof(double)*N);
	f2 = (double*) malloc(sizeof(double)*N1);

	char file_name[400];
	sprintf(file_name, "./Flow0_0.00e+00_8.67e-02_1.03e+00_dt_0.0005.dat");
	FILE *fp=fopen(file_name, "r");

	if(fp!=NULL) {
		int itmp;
		double tmp;
		double ff;
		int ii;
		for(i=0;i<N;i++)
		{
			fscanf(fp, "%d %le %le %le %le %le %le %le %le %le %le %le %le %le \n", &ii, &tmp, &tmp,  &tmp, &ff, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp);
			//fscanf(fp, "%le %le %le %le %le %le %le %le %le %le \n", &tmp, &tmp, &tmp,  &ff, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp);
			if (ii!=i+t0) i=i-1;
			f1[i]=ff;
		} 
		fclose(fp);
	}


	double f1_sum;
	for(i=0;i<N;i++)
	{
		f1_sum+=f1[i];
	}
	f1_sum/=(double)N;

	for(i=0;i<N;i++)
	{
		f1[i]-=f1_sum;
	}

	double uu=0.0;	
	for(i=0;i<N;i++)
	{
		uu+=f1[i]*f1[i];
	}
	uu/=N;

	FILE *f;
	char filen[80];
				
	sprintf(filen, "./In_u.dat");
	f = fopen(filen, "w");
	for(i=0;i<N1;i++)
	{
		fprintf(f, "%d %.7e \n", i, f1[i]);
	}
	fclose(f);

		
	fftw_complex *uin, *uin_1;
	fftw_complex *uout, *uout_1;

	uin = (fftw_complex*) fftw_malloc(N1*sizeof(fftw_complex));
	uout = (fftw_complex*) fftw_malloc(N1*sizeof(fftw_complex));
	uin_1 = (fftw_complex*) fftw_malloc(N1*sizeof(fftw_complex));
	uout_1 = (fftw_complex*) fftw_malloc(N1*sizeof(fftw_complex));

	fftw_plan pu;
	pu=fftw_plan_dft_1d(N1, uin, uout, FFTW_FORWARD, FFTW_ESTIMATE);

	fftw_plan pu_1;
	pu_1=fftw_plan_dft_1d(N1, uin_1, uout_1, FFTW_BACKWARD, FFTW_ESTIMATE);

	for(i=0;i<N1;i++) 
	{
		f2[i]=0.0;
	}

	for(k=0;k<KK;k++) 
	{

		int istart=k*Overlap;

		for(i=0;i<N1;i++) 
		{
			double t=(i-(N1-1)/2)/((N1-1)/2);
			double w=1-t*t;		
			uin[i][0]=f1[i+istart]*w;
			uin[i][1]=f1[i+istart]*w;
		}

		fftw_execute(pu); 

		for(i=0;i<N1;i++) 
		{
			uout[i][0]=uout[i][0]/((double)N1);
			uout[i][1]=uout[i][1]/((double)N1);
		}

		for(i=0;i<N1;i++)
		{
			double T=dt*N1;
			double rr=T*(uout[i][0]*uout[i][0]+uout[i][1]*uout[i][1])/uu;
			f2[i]+=rr;
		}

	}


	sprintf(filen, "./Out.dat");
	f = fopen(filen, "w");
	for(i=0;i<N1;i++)
	{
		fprintf(f, "%.7e %.7e \n", (double)i/(dt*(double)N1), f2[i]/((double)KK+1.e-19));
	}
	fclose(f);


	/*
	for(i=0;i<N1;i++) 
	{
		uin_1[i][0]=uout[i][0];
		uin_1[i][1]=uout[i][1];
	}

	fftw_execute(pu_1); 

	sprintf(filen, "./Out_1.dat");
	f = fopen(filen, "w");
	for(i=0;i<N1;i++)
	{
		fprintf(f, "%d %.7e %.7e\n", i, uout_1[i][0], uout_1[i][1]);
	}
	fclose(f);
	*/
	

	fftw_destroy_plan(pu);
	fftw_destroy_plan(pu_1);

	fftw_free(uin);
	fftw_free(uout);
	fftw_free(uin_1);
	fftw_free(uout_1);

	free(x);
	free(f1);
	free(f2);
}

double randn_notrig() {
	static bool deviateAvailable=false;   
	static float storedDeviate;
	double polar, rsquared, var1, var2;
	double mu=0.0; double sigma=1.0;
	if (!deviateAvailable) {
		do {
			var1=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			var2=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			rsquared=var1*var1+var2*var2;
		} while ( rsquared>=1.0 || rsquared == 0.0);

		polar=sqrt(-2.0*log(rsquared)/rsquared);
		storedDeviate=var1*polar;
		deviateAvailable=true;
		return var2*polar*sigma + mu;
	}
	else {
		deviateAvailable=false;
		return storedDeviate*sigma + mu;
	}
}


