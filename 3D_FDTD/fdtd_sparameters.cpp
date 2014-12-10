#include <stdio.h>
#include <cmath>
#include "fdtd.hpp"
#include "fdtd_sparameters.hpp"
#include <string.h>



void computeSParameters(SIMULATION* simulation)
{
    int n, nmax;
    double df;
    double fmax = 10e9;
    double freq;

    double x_i[131072], y_i[131072];
    double x_v[131072], y_v[131072];

    double tmp1,tmp2;
    complx a1,b1;

    short res;

    double dt = simulation->dt;
    int nstop  = simulation->Nt;


    FILE* fp_imp = fopen("./sparameters/zin.fd","w");
    FILE* fp_s11  = fopen("./sparameters/s11.fd","w");
    FILE* fp_time  = fopen("./sparameters/excitation.fd","w");

    list <PORT>::iterator port = simulation->list_of_ports.begin();

    /* read the dats in from the file */
    for (n=0; n<nstop; n++)
    {

        x_v[n] = port->v_in[n];
        x_i[n] = port->i_in[n];
        fprintf(fp_time,"%e %e\n",port->v_in[n],port->i_in[n]);

        y_i[n] = 0.0;
        y_v[n] = 0.0;
    }
    fclose(fp_time);

    /* Pad Zeros */
    for (n=(nstop-1); n<131072; n++)
    {
        x_i[n]=0.0;
        y_i[n]=0.0;
        x_v[n]=0.0;
        y_v[n]=0.0;
    }

    FFT(1,17,x_i,y_i);
    FFT(1,17,x_v,y_v);
    df = 1.0 / (dt*131072);

    complx fi;
    complx fv;
    complx imp;
    complx zo;
    complx tau;
    for (n=0; n<131072; n++)
    {

        freq = n*df;
        zo = complx(50,0);

        fi = complx(x_i[n],y_i[n]);
        fv = complx(x_v[n],y_v[n]);
        imp = fv / fi;
        tau = (imp - zo) / (imp + zo);

        double s11 = 20*log10(abs(tau));

        if ((freq > 500e6) && (freq < 5e9))
        {
            fprintf(fp_imp,"%e %e %e\n", freq, real(imp),imag(imp));
            fprintf(fp_s11,"%e %e\n", freq, s11);
        }

        /* if its at the max freq - break out */
        if (freq>fmax) break;

    }

    fclose(fp_imp);
    fclose(fp_s11);


}


void FFT(int dir, long m, double *x, double *y)

{
    long n,i,i1,j,k,i2,l,l1,l2;
    double c1,c2,tx,ty,t1,t2,u1,u2,z;

    /* Calculate the number of points */
    n = 1;
    for (i=0; i<m; i++)
        n *= 2;

    /* Do the bit reversal */
    i2 = n >> 1;
    j = 0;
    for (i=0; i<n-1; i++) {
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
    for (l=0; l<m; l++) {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.0;
        u2 = 0.0;
        for (j=0; j<l1; j++) {
            for (i=j; i<n; i+=l2) {
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
        for (i=0; i<n; i++) {
            x[i] /= n;
            y[i] /= n;
        }
    }

}


