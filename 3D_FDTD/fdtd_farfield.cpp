#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fdtd.hpp"
//#include "fdtd_farfield.h"

void CLOSED_SURFACE::perform_DFT(GRID* grid,int n,double dt)
{
    double Exavg, Eyavg, Ezavg;
    double Hxavg, Hyavg, Hzavg;
    double omega = 2 * pi * frequency;

    double* Ex = grid->Ex;
    double* Ey = grid->Ey;
    double* Ez = grid->Ez;

    double* Hx = grid->Hx;
    double* Hy = grid->Hy;
    double* Hz = grid->Hz;

    /* E-Field Sizes */
    int exnx = grid->exnx;
    int exny = grid->exny;
    int exnz = grid->exnz;

    int eynx = grid->eynx;
    int eyny = grid->eyny;
    int eynz = grid->eynz;

    int eznx = grid->eznx;
    int ezny = grid->ezny;
    int eznz = grid->eznz;

    /* H-Field Sizes */
    int hxnx = grid->hxnx;
    int hxny = grid->hxny;
    int hxnz = grid->hxnz;

    int hynx = grid->hynx;
    int hyny = grid->hyny;
    int hynz = grid->hynz;

    int hznx = grid->hznx;
    int hzny = grid->hzny;
    int hznz = grid->hznz;


    int i,j,k;

    complx factor = complx(cos(-omega * n * dt),sin(-omega * n * dt)) * dt;

    for (int n = 0; n < 2; n++)
    {

        /* X Plates */
        if (n == 0) {
            i = is;
        }
        else {
            i = ie;
        }

        for (j = js; j < je; j++)
        {
            for (k = ks; k < ke; k++)
            {

                // Calc. Average E fields
                Eyavg =   (Ey[i * eyny * eynz + j * eynz + k] + Ey[i * eyny * eynz + j * eynz + (k+1)]) / 2.0;
                Ezavg =	  (Ez[i * ezny * eznz + j * eznz + k] + Ez[i * ezny * eznz + (j+1) * eznz + k]) / 2.0;

                Hyavg = (Hy[i * hyny * hynz + j * hynz + k]
                         + Hy[i * hyny * hynz + (j+1) * hynz + k]
                         + Hy[(i-1) * hyny * hynz + j * hynz + k]
                         + Hy[(i-1) * hyny * hynz + (j+1) * hynz + k]) / 4.0;

                Hzavg = (Hz[i * hzny * hznz + j * hznz + k]
                         + Hz[i * hzny * hznz + j * hznz + (k+1)]
                         + Hz[(i-1) * hzny * hznz + j * hznz + k]
                         + Hz[(i-1) * hzny * hznz + j * hznz + (k+1)]) / 4.0;

                x_surfaces[n].Ey[(k - ks) + (j - js) * kw] += factor * Eyavg;
                x_surfaces[n].Ez[(k - ks) + (j - js) * kw] += factor * Ezavg;
                x_surfaces[n].Hy[(k - ks) + (j - js) * kw] += factor * Hyavg;
                x_surfaces[n].Hz[(k - ks) + (j - js) * kw] += factor * Hzavg;
            }
        }


        /* Y Plates */
        if (n == 0) {
            j = js;
        }
        else {
            j = je;
        }
        for (i = is; i < ie; i++)
        {
            for (k = ks; k < ke; k++)
            {

                // Calc. Average E fields
                Exavg =   (Ex[i * exny * exnz + j * exnz + k] + Ex[i * exny * exnz + j * exnz + (k+1)]) / 2.0;
                Ezavg =	  (Ez[i * ezny * eznz + j * eznz + k] + Ez[(i+1) * ezny * eznz + j * eznz + k]) / 2.0;

                Hxavg = (Hx[i * hxny * hxnz + j * hxnz + k]
                         + Hx[(i+1) * hxny * hxnz + j * hxnz + k]
                         + Hx[i * hxny * hxnz + (j-1) * hxnz + k]
                         + Hx[(i+1) * hxny * hxnz + (j-1) * hxnz + k]) / 4.0;

                Hzavg = (Hz[i * hzny * hznz + j * hznz + k]
                         + Hz[i * hzny * hznz + j * hznz + (k+1)]
                         + Hz[i * hzny * hznz + (j-1) * hznz + k]
                         + Hz[i * hzny * hznz + (j-1) * hznz + (k+1)]) / 4.0;


                y_surfaces[n].Ex[(k - ks) + (i - is) * kw] += factor * Exavg;
                y_surfaces[n].Ez[(k - ks) + (i - is) * kw] += factor * Ezavg;
                y_surfaces[n].Hx[(k - ks) + (i - is) * kw] += factor * Hxavg;
                y_surfaces[n].Hz[(k - ks) + (i - is) * kw] += factor * Hzavg;
            }
        }

        /* Z Plates */
        if (n == 0) {
            k = ks;
        }
        else {
            k = ke;
        }
        for (i = is; i < ie; i++)
        {
            for (j = js; j < je; j++)
            {

                // Calc. Average E fields
                Exavg =   (Ex[i * exny * exnz + j * exnz + k] + Ex[i * exny * exnz + (j+1) * exnz + k]) / 2.0;
                Eyavg =   (Ey[i * eyny * eynz + j * eynz + k] + Ey[(i+1) * eyny * eynz + j * eynz + k]) / 2.0;

                Hxavg = ( Hx[i * hxny * hxnz + j * hxnz + k]
                          + Hx[ (i+1) * hxny * hxnz + j * hxnz + k]
                          + Hx[i * hxny * hxnz + j * hxnz + (k-1)]
                          + Hx[(i+1) * hxny * hxnz + j * hxnz + (k-1)]) / 4.0;
                Hyavg = ( Hy[i * hyny * hynz + j * hynz + k]
                          + Hy[i * hyny * hynz + (j+1) * hynz + k]
                          +  Hy[i * hyny * hynz + j * hynz + (k-1)]
                          + Hy[i * hyny * hynz + (j+1) * hynz + (k-1)]) / 4.0;



                z_surfaces[n].Ex[(j - js) + (i - is) * jw] += factor * Exavg;
                z_surfaces[n].Ey[(j - js) + (i - is) * jw] += factor * Eyavg;
                z_surfaces[n].Hx[(j - js) + (i - is) * jw] += factor * Hxavg;
                z_surfaces[n].Hy[(j - js) + (i - is) * jw] += factor * Hyavg;
            }
        }

    }

}

void CLOSED_SURFACE::compute_surface_currents()
{

    for (int j = 0; j < jw; j++)
    {
        for (int k = 0; k < kw; k++)
        {

            x_surfaces[0].Jy[k + j * kw] =  x_surfaces[0].Hz[k + j * kw];
            x_surfaces[0].Jz[k + j * kw] = -x_surfaces[0].Hy[k + j * kw];
            x_surfaces[0].My[k + j * kw] = -x_surfaces[0].Ez[k + j * kw];
            x_surfaces[0].Mz[k + j * kw] =  x_surfaces[0].Ey[k + j * kw];

            x_surfaces[1].Jy[k + j * kw] = -x_surfaces[1].Hz[k + j * kw];
            x_surfaces[1].Jz[k + j * kw] =  x_surfaces[1].Hy[k + j * kw];
            x_surfaces[1].My[k + j * kw] =  x_surfaces[1].Ez[k + j * kw];
            x_surfaces[1].Mz[k + j * kw] = -x_surfaces[1].Ey[k + j * kw];
        }
    }

    for (int i = 0; i < iw; i++)
    {
        for (int k = 0; k < kw; k++)
        {

            // Y = Ymin plate
            y_surfaces[0].Jx[k + i * kw] = -y_surfaces[0].Hz[k + i * kw];
            y_surfaces[0].Jz[k + i * kw] =  y_surfaces[0].Hx[k + i * kw];
            y_surfaces[0].Mx[k + i * kw] =  y_surfaces[0].Ez[k + i * kw];
            y_surfaces[0].Mz[k + i * kw] = -y_surfaces[0].Ex[k + i * kw];


            // Y = Ymax plate
            y_surfaces[1].Jx[k + i * kw] =  y_surfaces[1].Hz[k + i * kw];
            y_surfaces[1].Jz[k + i * kw] = -y_surfaces[1].Hx[k + i * kw];
            y_surfaces[1].Mx[k + i * kw] = -y_surfaces[1].Ez[k + i * kw];
            y_surfaces[1].Mz[k + i * kw] =  y_surfaces[1].Ex[k + i * kw];
        }
    }

    for (int i = 0; i < iw; i++)
    {
        for (int j = 0; j < jw; j++)
        {

            // Z = Zmin plate
            z_surfaces[0].Jx[j + i * jw] = z_surfaces[0].Hy[j + i * jw];
            z_surfaces[0].Jy[j + i * jw] = -z_surfaces[0].Hx[j + i * jw];
            z_surfaces[0].Mx[j + i * jw] = -z_surfaces[0].Ey[j + i * jw];
            z_surfaces[0].My[j + i * jw] = z_surfaces[0].Ex[j + i * jw];


            // Z = Zmax plate
            z_surfaces[1].Jx[j + i * jw] = -z_surfaces[1].Hy[j + i * jw];
            z_surfaces[1].Jy[j + i * jw] = z_surfaces[1].Hx[j + i * jw];
            z_surfaces[1].Mx[j + i * jw] = z_surfaces[1].Ey[j + i * jw];
            z_surfaces[1].My[j + i * jw] = -z_surfaces[1].Ex[j + i * jw];

        }
    }
}

void CLOSED_SURFACE::compute_radiated_power(list<PORT>::iterator port,double dx, double dy, double dz)
{

    /* Time Average Poynting Vectors */
    double Pvector_XPlates[jw][kw][2];
    double Pvector_YPlates[iw][kw][2];
    double Pvector_ZPlates[iw][jw][2];

    /* Time Average Poynting Vectors */
    double Pvector_Total_XPlates[2];
    double Pvector_Total_YPlates[2];
    double Pvector_Total_ZPlates[2];

    Pvector_Total_XPlates[0] = 0;
    Pvector_Total_XPlates[1] = 0;
    Pvector_Total_YPlates[0] = 0;
    Pvector_Total_YPlates[1] = 0;
    Pvector_Total_ZPlates[0] = 0;
    Pvector_Total_ZPlates[1] = 0;


    /* X Plates */
    for (int j = 0; j < jw; j++)
    {
        for (int k = 0; k < kw; k++)
        {
            Pvector_XPlates[j][k][0] = 0.5*real(conj(x_surfaces[0].Hz[k + kw * j])*x_surfaces[0].Ey[k + kw * j] - x_surfaces[0].Ez[k + kw * j]*conj(x_surfaces[0].Hy[k + kw * j]));
            Pvector_XPlates[j][k][1] = 0.5*real(conj(x_surfaces[1].Hy[k + kw * j])*x_surfaces[1].Ez[k + kw * j] - x_surfaces[1].Ey[k + kw * j]*conj(x_surfaces[1].Hz[k + kw * j]));

            Pvector_Total_XPlates[0] = Pvector_Total_XPlates[0] + Pvector_XPlates[j][k][0] * dy * dz;
            Pvector_Total_XPlates[1] = Pvector_Total_XPlates[1] + Pvector_XPlates[j][k][1] * dy * dz;
        }
    }

    /* Y Plates */
    for (int i = 0; i < iw; i++)
    {
        for (int k = 0; k < kw; k++)
        {
            Pvector_YPlates[i][k][0] = 0.5*real(conj(y_surfaces[0].Hx[k + kw * i])*y_surfaces[0].Ez[k + kw * i] - y_surfaces[0].Ex[k + kw * i]*conj(y_surfaces[0].Hz[k + kw * i]));
            Pvector_YPlates[i][k][1] = 0.5*real(conj(y_surfaces[1].Hz[k + kw * i])*y_surfaces[1].Ex[k + kw * i] - y_surfaces[1].Ez[k + kw * i]*conj(y_surfaces[1].Hx[k + kw * i]));

            Pvector_Total_YPlates[0] = Pvector_Total_YPlates[0] + Pvector_YPlates[i][k][0] * dx * dz;
            Pvector_Total_YPlates[1] = Pvector_Total_YPlates[1] + Pvector_YPlates[i][k][1] * dx * dz;
        }
    }

    /* Z Plates */
    for (int i = 0; i < iw; i++)
    {
        for (int j = 0; j < jw; j++)
        {
            Pvector_ZPlates[i][j][0] = 0.5*real(conj(z_surfaces[0].Hy[j + kw * i])*z_surfaces[0].Ex[j + kw * i] - z_surfaces[0].Ey[j + jw * i]*conj(z_surfaces[0].Hx[j + kw * i]));
            Pvector_ZPlates[i][j][1] = 0.5*real(conj(z_surfaces[1].Hx[j + kw * i])*z_surfaces[1].Ey[j + kw * i] - z_surfaces[1].Ex[j + jw * i]*conj(z_surfaces[1].Hy[j + kw * i]));

            Pvector_Total_ZPlates[0] = Pvector_Total_ZPlates[0] + Pvector_ZPlates[i][j][0] * dx * dy;
            Pvector_Total_ZPlates[1] = Pvector_Total_ZPlates[1] + Pvector_ZPlates[i][j][1] * dx * dy;

        }
    }

    prad = Pvector_Total_XPlates[0] + Pvector_Total_XPlates[1] + Pvector_Total_YPlates[0] + Pvector_Total_YPlates[1] +Pvector_Total_ZPlates[0] + Pvector_Total_ZPlates[1];
    prad = sqrt( pow(prad,2) );
    pin = 0.5*real(port->V_in[label] * conj(port->I_in[label]));
    pin = sqrt( pin * pin);
    efficiency = prad / pin;
}
void CLOSED_SURFACE::compute_farfield_xzplane(double dx,double dy, double dz,int label)
{
    double phi  = 0;
    double theta;
    double xf,yf,zf;
    double wavenum;
    double krcospsi;
    complx factor;
    double idis,jdis,kdis;

    complx Nth[360];
    complx Nph[360];
    complx Lth[360];
    complx Lph[360];

    double omega = 2 * pi * frequency;


    wavenum = omega  / c;


    for (int th = 0; th < 360 ; th++)
    {
        theta = (double)th;
        theta = theta * pi/180;
        phi = phi * pi/180;

        xf = sin(theta) * cos(phi);
        yf = sin(theta) * sin(phi);
        zf = cos(theta);

        /*********************************************************/
        /*              Intergration of the X Plates		 */
        /*********************************************************/
        for (int n = 0; n < 2 ; n++)
        {
            if (n == 0) {
                idis = is - xc;
            }
            else {
                idis = ie - xc;
            }
            for (int j = 0; j < jw; j++)
            {
                for (int k = 0; k < kw; k++)
                {
                    jdis = js + j - yc + 0.5;
                    kdis = ks + k - zc + 0.5;

                    krcospsi = wavenum * ((idis * dx * xf)	+ (jdis * dy * yf) + (kdis * dz * zf));
                    factor = complx(cos(krcospsi),  sin(krcospsi))*(complx)(dy*dz);

                    // --- <N> Far Field Vector Potential ---
                    Nth[th] += factor*(x_surfaces[n].Jy[k + kw * j] * cos(theta)*sin(phi) - x_surfaces[n].Jz[k + kw * j] * sin(theta));
                    Nph[th] += factor*(x_surfaces[n].Jy[k + kw * j] * cos(phi));

                    // --- <L> Far Field Vector Potential ---
                    Lth[th] += factor*(x_surfaces[n].My[k + kw * j] * cos(theta)*sin(phi) - x_surfaces[n].Mz[k + kw * j] * sin(theta));
                    Lph[th] += factor*(x_surfaces[n].My[k + kw * j] * cos(phi));


                }
            }
        }

        /*********************************************************/
        /*              Intergration of the Y Plates		 */
        /*********************************************************/

        for (int n = 0; n < 2 ; n++)
        {
            if (n == 0) {
                jdis = js - yc;
            }
            else {
                jdis = je- yc;
            }
            for (int i = 0; i < iw; i++)
            {
                for (int k = 0; k < kw; k++)
                {
                    idis = is + i - xc + 0.5;
                    kdis = ks + k - zc + 0.5;

                    krcospsi = wavenum * ((idis * dx * xf)	+ (jdis * dy * yf) + (kdis * dz * zf));
                    factor = complx(cos(krcospsi),  sin(krcospsi))*(complx)(dx*dz);

                    // --- <N> Far Field Vector Potential ---
                    Nth[th] += factor*(y_surfaces[n].Jx[k + kw * i] * cos(theta)*cos(phi) - y_surfaces[n].Jz[k + kw * i] * sin(theta));
                    Nph[th] += factor*(-y_surfaces[n].Jx[k + kw * i] * sin(phi));

                    // --- <L> Far Field Vector Potential ---
                    Lth[th] += factor*(y_surfaces[n].Mx[k + kw * i] * cos(theta)*cos(phi) - y_surfaces[n].Mz[k + kw * i] * sin(theta));
                    Lph[th] += factor*(-y_surfaces[n].Mx[k + kw * i] * sin(phi));


                }
            }
        }
        /*********************************************************/
        /*              Intergration of the Z Plates		 */
        /*********************************************************/

        for (int n = 0; n < 2 ; n++)
        {
            if (n == 0) {
                kdis = ks- zc;
            }
            else {
                kdis = ke- zc;
            }
            for (int i = 0; i < iw; i++)
            {
                for (int j = 0; j < jw; j++)
                {
                    idis = is + i - xc + 0.5;
                    jdis = js + j - yc + 0.5;

                    krcospsi = wavenum * ((idis * dx * xf)	+ (jdis * dy * yf) + (kdis * dz * zf));
                    factor = complx(cos(krcospsi),  sin(krcospsi))*(complx)(dx*dz);

                    // --- <N> Far Field Vector Potential ---
                    Nth[th] += factor*(z_surfaces[n].Jx[j + jw * i] * cos(theta)*cos(phi) + z_surfaces[n].Jy[j + kw * i] * cos(theta) * sin(phi));
                    Nph[th] += factor*(z_surfaces[n].Jx[j + jw * i] * sin(phi) + z_surfaces[n].Jy[j + jw * i] * cos(phi) );
                    // --- <L> Far Field Vector Potential ---
                    Lth[th] += factor*(z_surfaces[n].Mx[j + jw * i] * cos(theta)*cos(phi) + z_surfaces[n].My[j + kw * i] * cos(theta) *  sin(phi));
                    Lph[th] += factor*(z_surfaces[n].Mx[j + jw * i] * sin(phi) + z_surfaces[n].My[j + jw * i] * cos(phi) );
                }
            }
        }
    }

    double factor2 = (eta*omega * omega)/(8.0*c*c);
    char fname[100];
    sprintf(fname,"./farfield/xzRad_%i.fd",label);
    FILE* fp_out = fopen(fname,"w");
    for (int n = 0; n <360; n++)
    {
        double field =  factor2 * (norm(Nth[n] + Lph[n]/eta)+norm(Nph[n] - Lth[n]/eta));
        fprintf(fp_out,"%e \n",field);
    }
    fclose(fp_out);
}
void CLOSED_SURFACE::compute_farfield_xyplane(double dx,double dy, double dz,int label)
{
    double phi;
    double theta = 90;
    double xf,yf,zf;
    double wavenum;
    double krcospsi;
    complx factor;
    double idis,jdis,kdis;

    complx Nth[360];
    complx Nph[360];
    complx Lth[360];
    complx Lph[360];

    double omega = 2 * pi * frequency;


    wavenum = omega  / c;


    for (int ph = 0; ph < 360 ; ph++)
    {
        phi = (double)ph;
        theta = theta * pi/180;
        phi = phi * pi/180;

        xf = sin(theta) * cos(phi);
        yf = sin(theta) * sin(phi);
        zf = cos(theta);

        /*********************************************************/
        /*              Intergration of the X Plates		 */
        /*********************************************************/
        for (int n = 0; n < 2 ; n++)
        {
            if (n == 0) {
                idis = is - xc;
            }
            else {
                idis = ie - xc;
            }
            for (int j = 0; j < jw; j++)
            {
                for (int k = 0; k < kw; k++)
                {
                    jdis = js + j - yc + 0.5;
                    kdis = ks + k - zc + 0.5;

                    krcospsi = wavenum * ((idis * dx * xf)	+ (jdis * dy * yf) + (kdis * dz * zf));
                    factor = complx(cos(krcospsi),  sin(krcospsi))*(complx)(dy*dz);

                    // --- <N> Far Field Vector Potential ---
                    Nth[ph] += factor*(x_surfaces[n].Jy[k + kw * j] * cos(theta)*sin(phi) - x_surfaces[n].Jz[k + kw * j] * sin(theta));
                    Nph[ph] += factor*(x_surfaces[n].Jy[k + kw * j] * cos(phi));

                    // --- <L> Far Field Vector Potential ---
                    Lth[ph] += factor*(x_surfaces[n].My[k + kw * j] * cos(theta)*sin(phi) - x_surfaces[n].Mz[k + kw * j] * sin(theta));
                    Lph[ph] += factor*(x_surfaces[n].My[k + kw * j] * cos(phi));


                }
            }
        }

        /*********************************************************/
        /*              Intergration of the Y Plates		 */
        /*********************************************************/

        for (int n = 0; n < 2 ; n++)
        {
            if (n == 0) {
                jdis = js - yc;
            }
            else {
                jdis = je- yc;
            }
            for (int i = 0; i < iw; i++)
            {
                for (int k = 0; k < kw; k++)
                {
                    idis = is + i - xc + 0.5;
                    kdis = ks + k - zc + 0.5;

                    krcospsi = wavenum * ((idis * dx * xf)	+ (jdis * dy * yf) + (kdis * dz * zf));
                    factor = complx(cos(krcospsi),  sin(krcospsi))*(complx)(dx*dz);

                    // --- <N> Far Field Vector Potential ---
                    Nth[ph] += factor*(y_surfaces[n].Jx[k + kw * i] * cos(theta)*cos(phi) - y_surfaces[n].Jz[k + kw * i] * sin(theta));
                    Nph[ph] += factor*(-y_surfaces[n].Jx[k + kw * i] * sin(phi));

                    // --- <L> Far Field Vector Potential ---
                    Lth[ph] += factor*(y_surfaces[n].Mx[k + kw * i] * cos(theta)*cos(phi) - y_surfaces[n].Mz[k + kw * i] * sin(theta));
                    Lph[ph] += factor*(-y_surfaces[n].Mx[k + kw * i] * sin(phi));


                }
            }
        }
        /*********************************************************/
        /*              Intergration of the Z Plates		 */
        /*********************************************************/

        for (int n = 0; n < 2 ; n++)
        {
            if (n == 0) {
                kdis = ks- zc;
            }
            else {
                kdis = ke- zc;
            }
            for (int i = 0; i < iw; i++)
            {
                for (int j = 0; j < jw; j++)
                {
                    idis = is + i - xc + 0.5;
                    jdis = js + j - yc + 0.5;

                    krcospsi = wavenum * ((idis * dx * xf)	+ (jdis * dy * yf) + (kdis * dz * zf));
                    factor = complx(cos(krcospsi),  sin(krcospsi))*(complx)(dx*dz);

                    // --- <N> Far Field Vector Potential ---
                    Nth[ph] += factor*(z_surfaces[n].Jx[j + jw * i] * cos(theta)*cos(phi) + z_surfaces[n].Jy[j + kw * i] * cos(theta) * sin(phi));
                    Nph[ph] += factor*(z_surfaces[n].Jx[j + jw * i] * sin(phi) + z_surfaces[n].Jy[j + jw * i] * cos(phi) );
                    // --- <L> Far Field Vector Potential ---
                    Lth[ph] += factor*(z_surfaces[n].Mx[j + jw * i] * cos(theta)*cos(phi) + z_surfaces[n].My[j + kw * i] * cos(theta) *  sin(phi));
                    Lph[ph] += factor*(z_surfaces[n].Mx[j + jw * i] * sin(phi) + z_surfaces[n].My[j + jw * i] * cos(phi) );
                }
            }
        }
    }

    double factor2 = (eta*omega * omega)/(8.0*c*c);
    char fname[100];
    sprintf(fname,"./farfield/xyRad_%i.fd",label);
    FILE* fp_out = fopen(fname,"w");
    for (int n = 0; n <360; n++)
    {
        double field =  factor2 * (norm(Nth[n] + Lph[n]/eta)+norm(Nph[n] - Lth[n]/eta));
        fprintf(fp_out,"%e \n",field);
    }
    fclose(fp_out);
}


void CLOSED_SURFACE::allocate_surface_memory(int Nx, int Ny, int Nz, int index)
{

    is = 4;
    ie = Nx-4;
    js = 4;
    je = Ny-4;
    ks = 4;
    ke = Nz-4;

    iw = ie - is;
    jw = je - js;
    kw = ke - ks;
    label = index;



    /* The Centre of the Closed Surface */
    xc = Nx/2;
    yc = Ny/2;
    zc = Nz/2;

    for (int n = 0; n < 2 ; n++)
    {
        /* x Surfaces */
        x_surfaces[n].Ey = (complx*) malloc(iw * kw * sizeof(complx));
        x_surfaces[n].Ez = (complx*) malloc(iw * kw * sizeof(complx));

        x_surfaces[n].Hy = (complx*) malloc(iw * kw * sizeof(complx));
        x_surfaces[n].Hz = (complx*) malloc(iw * kw * sizeof(complx));

        x_surfaces[n].Jy = (complx*) malloc(iw * kw * sizeof(complx));
        x_surfaces[n].Jz = (complx*) malloc(iw * kw * sizeof(complx));

        x_surfaces[n].My = (complx*) malloc(iw * kw * sizeof(complx));
        x_surfaces[n].Mz = (complx*) malloc(iw * kw * sizeof(complx));

        /* y Surfaces */
        y_surfaces[n].Ex = (complx*) malloc(jw * kw * sizeof(complx));
        y_surfaces[n].Ez = (complx*) malloc(jw * kw * sizeof(complx));

        y_surfaces[n].Hx = (complx*) malloc(jw * kw * sizeof(complx));
        y_surfaces[n].Hz = (complx*) malloc(jw * kw * sizeof(complx));

        y_surfaces[n].Jx = (complx*) malloc(jw * kw * sizeof(complx));
        y_surfaces[n].Jz = (complx*) malloc(jw * kw * sizeof(complx));

        y_surfaces[n].Mx = (complx*) malloc(jw * kw * sizeof(complx));
        y_surfaces[n].Mz = (complx*) malloc(jw * kw * sizeof(complx));

        /* z Surfaces */
        z_surfaces[n].Ex = (complx*) malloc(iw * jw * sizeof(complx));
        z_surfaces[n].Ey = (complx*) malloc(iw * jw * sizeof(complx));

        z_surfaces[n].Hx = (complx*) malloc(iw * jw * sizeof(complx));
        z_surfaces[n].Hy = (complx*) malloc(iw * jw * sizeof(complx));

        z_surfaces[n].Jx = (complx*) malloc(iw * jw * sizeof(complx));
        z_surfaces[n].Jy = (complx*) malloc(iw * jw * sizeof(complx));

        z_surfaces[n].Mx = (complx*) malloc(iw * jw * sizeof(complx));
        z_surfaces[n].My = (complx*) malloc(iw * jw * sizeof(complx));

        complx null(0,0);
        for (int i=0; i<iw; i++)
        {
            for (int k=0; k<kw; k++)
            {
                x_surfaces[n].Ey[k + i * kw] = null;
                x_surfaces[n].Ez[k + i * kw] = null;
                x_surfaces[n].Hy[k + i * kw] = null;
                x_surfaces[n].Hz[k + i * kw] = null;


                y_surfaces[n].Ex[k + i * kw] = null;
                y_surfaces[n].Ez[k + i * kw] = null;
                y_surfaces[n].Hx[k + i * kw] = null;
                y_surfaces[n].Hz[k + i * kw] = null;

                z_surfaces[n].Ex[k + i * kw] = null;
                z_surfaces[n].Ey[k + i * kw] = null;
                z_surfaces[n].Hx[k + i * kw] = null;
                z_surfaces[n].Hy[k + i * kw] = null;
            }
        }

    }

}

