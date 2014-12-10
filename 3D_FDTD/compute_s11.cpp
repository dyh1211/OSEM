
/******************* FDTD 3D - S11 Calculation************* *************/
/*                                                                      */
/*  After specifying the input files - will calculate the FFT of each   */
/*  and then the S parameters as specified...                           */
/*                                                                      */
/*  By :                                                                */
/*                                                                      */
/*  Additional files :                                                  */
/*                                                                      */
/************************************************************************/

#include <stdio.h>
#include <math.h>
#include "fdtd_complx.hpp"
#include <string.h>


void FFT(int dir, long m,double *x, double *y);

int main(int argc, char* argv[])
{
    
    FILE *fp_dt, *fp_nstop;     // Info read from files.
    FILE *obs_i, *obs_v;
    FILE *pulse_i, *pulse_v;
    FILE *fp_impi, *fp_impr, *fp_ffti, *fp_fftv;
    FILE *fp_s11;
    char fname[100];
    int n, nmax;
    double df;
    double fmax = 10e9;
    double freq =1e9;
    double dt = 0;
    double x_i[131072], y_i[131072];
    double x_v[131072], y_v[131072];
    double tmp1,tmp2;   
    
    short res;


    /* Open all the files */

    fp_dt = fopen("dt.fd","r");
    fp_nstop = fopen("nstop.fd","r");
    
    obs_v = fopen(argv[1],"r");

    fp_impr = fopen("imp_real.fd","w");
    fp_impi = fopen("imp_imag.fd","w");

    strcpy(fname,"s11_");
    strcat(fname,argv[1]);
    fp_s11 = fopen(fname,"w");

    fscanf(fp_dt,"%le\n",&dt);
    fscanf(fp_nstop,"%d\n",&nmax);

    printf("NMAX = %d\n",nmax);
    printf("DT   = %e\n",dt);

    
    /* read the dats in from the file */
    printf("Reading data file file - N = %d\n",nmax);
    for (n=0; n<nmax; n++)
    {
	fscanf(obs_v,"%le %le\n",&tmp2,&tmp1);

	/* get reflection only */
	x_i[n] = tmp1;
	x_v[n] = tmp2;

	y_i[n] = 0.0;
	y_v[n] = 0.0;
    }
    

    printf("Buffering the data..\n");
    
    for (n=(nmax-1); n<131072; n++)
      {
	x_i[n]=0.0;
	y_i[n]=0.0;
	x_v[n]=0.0;
	y_v[n]=0.0;
      }
    
    printf("FFT of current.. ");
    FFT(1,17,x_i,y_i);  
    printf("done\n");
	

    printf("FFT of voltage.. ");
    FFT(1,17,x_v,y_v);  
    printf("done\n");

    

    printf("Calculating impedance ...\n");

    df = 1.0 / (dt*131072);
    printf("df = %e\n",df);
    
    for (n=0; n<131072; n++)
    {

	freq = n*df;

//	tmp_i=sqrt(x_i[n]*x_i[n]+y_i[n]*y_i[n]);
//	tmp_v=sqrt(x_v[n]*x_v[n]+y_v[n]*y_v[n]);	

	complx fi;
	complx fv;
	complx imp;
	complx zo;
	complx tau;

	zo = complx(50,0);

	fi = complx(x_i[n],y_i[n]);
	fv = complx(x_v[n],y_v[n]);
	imp = fv / fi;
        tau = (imp - zo) / (imp + zo);

	double s11 = 20*log10(abs(tau));




	if ((freq > 100e6) && (freq < 5e9)) 
	{
	fprintf(fp_impr,"%e %e\n", freq, real(imp));
	fprintf(fp_impi,"%e %e\n", freq, imag(imp));
	fprintf(fp_s11,"%e %e\n", freq, s11);
	}

	/* if its at the max freq - break out */
	if (freq>fmax) break;

      }
    
   
    fclose(fp_dt);
    fclose(fp_nstop);

  //  fclose(obs_i);
    fclose(obs_v);

//    fclose(fp_ffti);
//    fclose(fp_fftv);

//    fclose(pulse_i);
//    fclose(pulse_v);


    fclose(fp_impi);
    fclose(fp_impr);

    fclose(fp_s11);

    return(0);

}  /******* end of main program ********/





/*
   This computes an in-place complex-to-complex FFT 
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform 
*/
void FFT(int dir, long m, double *x, double *y)

{
   long n,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   n = 1;
   for (i=0;i<m;i++) 
      n *= 2;

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0; 
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1; 
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1) 
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
    }
   
 }


