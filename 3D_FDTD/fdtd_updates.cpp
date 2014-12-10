#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fdtd.hpp"
#include "fdtd_updates.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <vector>
using namespace std;


void updateExcitation(GRID* grid, SIMULATION* simulation)
{
    int f;
    double 	ca, cbx, cby, cv;	// Constants for field update

    double* hz = grid->Hz;
    double* hy = grid->Hy;
    double* hx = grid->Hx;
    double* ex = grid->Ex;
    double* ey = grid->Ey;
    double* ez = grid->Ez;

    double  w_pulse = 40e-12;             /* Width of the excitation pulse */
    double  t0 = 4.0 * w_pulse;            /* time of the peak of the pulse */
    double  t = simulation->n * simulation->dt - (0.5*simulation->dt); // -.5 dt per. p.459 Taflove.


    int ezny = grid->ezny;
    int eznz = grid->eznz;

    int hxny = grid->hxny;
    int hxnz = grid->hxnz;

    int hyny = grid->hyny;
    int hynz = grid->hynz;


    for (list<PORT>::iterator it=simulation->list_of_ports.begin(); it!=simulation->list_of_ports.end(); ++it)
    {

        int i = it->i;
        int j = it->j;
        int k = it->k;
        double zo = it->zo;
        double exc;

        ca = (1 - (simulation->dt * simulation->dz) / (2*zo*eps0*simulation->dx*simulation->dy))/
             (1 + (simulation->dt * simulation->dz) / (2*zo*eps0*simulation->dx*simulation->dy));

        cbx = (simulation->dt/eps0)/(1+(simulation->dt*simulation->dz) / (2*zo*eps0*simulation->dx*simulation->dy));
        cby = (simulation->dt/eps0)/(1+(simulation->dt*simulation->dz) / (2*zo*eps0*simulation->dx*simulation->dy));
        cv =  (simulation->dt / (zo*eps0*simulation->dx*simulation->dy)) / (1+( (simulation->dt*simulation->dz) / (2*50*eps0*simulation->dx*simulation->dy) ));

        exc = it->signal(t,it->parameters[0],it->parameters[1],it->parameters[2]);

        ez[i * ezny * eznz + j * eznz +k] = ca * it->Ezn1+ cbx * ((hy[i * hyny * hynz + j * hynz +k] - hy[(i-1) * hyny * hynz + j * hynz +k]) / simulation->dx)- cby * ((hx[i * hxny * hxnz + j * hxnz +k] - hx[i * hxny * hxnz + (j-1) * hxnz +k]) / simulation->dy)+ cv * -exc;

        it->Ezn1 = ez[i * ezny * eznz + j * eznz +k];


        double current = (hx[i * hxny * hxnz + (j-1) * hxnz +k] - hx[i * hxny * hxnz + j * hxnz +k]) * simulation->dx + (hy[i * hyny * hynz + j * hynz +k] - hy[(i-1) * hyny * hynz + j * hynz +k])*simulation->dy;
        double voltage =  ez[i * ezny * eznz + j * eznz +k] * -simulation->dz;

        it->v_in[simulation->n]  = voltage;
        it->i_in[simulation->n]  = current;


        int f = 0;
        for  (list<double>::iterator freq=simulation->frequencies.begin(); freq!=simulation->frequencies.end(); ++freq)
        {
            double omega = 2 * pi * *freq;
            complx factor = complx(cos(-omega * simulation->n * simulation->dt),sin(-omega * simulation->n * simulation->dt)) * simulation->dt;
            it->V_in[f] += complx(voltage,0) * factor;
            it->I_in[f] += complx(current,0) * factor;
            f++;
        }
    }
}

void updateWires(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants)
{
    double ca,cbx,cby,cv;

    double* ex = grid->Ex;
    double* ey = grid->Ey;
    double* ez = grid->Ez;

    double* hx = grid->Hx;
    double* hy = grid->Hy;
    double* hz = grid->Hz;

    int exny = grid->exny;
    int exnz = grid->exnz;

    int eyny = grid->eyny;
    int eynz = grid->eynz;

    int ezny = grid->ezny;
    int eznz = grid->eznz;

    int hxny = grid->hxny;
    int hxnz = grid->hxnz;

    int hyny = grid->hyny;
    int hynz = grid->hynz;

    int hzny = grid->hzny;
    int hznz = grid->hznz;

    int* xe = grid->xe;
    int* ye = grid->ye;
    int* ze = grid->ze;

    int* xm = grid->xm;
    int* ym = grid->ym;
    int* zm = grid->zm;

    double* da = constants->da;
    double* dbx = constants->dbx;
    double* dbz = constants->dbz;


    for (list <WIRE>::iterator wir = simulation->list_of_wires.begin(); wir != simulation->list_of_wires.end(); wir++)
    {
        int i = wir->is;
        int j = wir->js;
        int ks = wir->ks;
        int ke = ks + wir->length - 1;
        double radius = wir->radius;

        /* << Hx Updates >> */
        /* Do Front of wire */
        for(int k = ks; k < ke; k++)
        {
            int material_idx = xm[i * hxny * hxnz + j * hxnz +k];
            double Ur = constants->Ur[material_idx];
            double da = 1;
            double db2 = (2 * simulation->dt) / (mu0*Ur * simulation->dy *log(simulation->dx / wir->radius) );

            hx[i * hxny * hxnz + j * hxnz +k] = da * wir->hx_buf[0][k-ks]
                                                + dbz[material_idx]  * ( ey[i * eyny * eynz + j * eynz + (k+1)] - ey[i * eyny * eynz + j * eynz +k])
                                                - db2  * ( ez[i * ezny * eznz + (j+1) * eznz + k]);
            wir->hx_buf[0][k-ks] = hx[i * hxny * hxnz + j * hxnz +k];
        }


        j = wir->js - 1;
        /* Do Back of wire */
        for(int k = ks; k < ke; k++)
        {
            int material_idx = xm[i * hxny * hxnz + (j-1) * hxnz +k];
            double Ur = constants->Ur[material_idx];
            double da = 1;
            double db2 = (2 * simulation->dt) / (mu0*Ur * simulation->dy *log(simulation->dy / wir->radius) );


            hx[i * hxny * hxnz + j * hxnz +k] = da * wir->hx_buf[1][k-ks]
                                                + dbz[material_idx]  * ( ey[i * eyny * eynz + j * eynz + (k+1)] - ey[i * eyny * eynz + j * eynz +k])
                                                + db2  * ( ez[i * ezny * eznz + j * eznz + k] );
            wir->hx_buf[1][k-ks] = hx[i * hxny * hxnz + j * hxnz +k];

        }

        /* Reset Position */
        i = wir->is;
        j = wir->js;
        ks = wir->ks;
        ke = ks + wir->length - 1;

        /* << Hy Updates >> */
        /* Do Front of wire */
        for(int k = ks; k < ke; k++)
        {
            int material_idx = ym[i * hyny * hynz + j * hynz +k];
            double Ur = constants->Ur[material_idx];
            double da = 1;
            double db2 = (2 * simulation->dt) / (mu0*Ur * simulation->dx *log(simulation->dx / wir->radius) );

            hy[i * hyny * hynz + j * hynz +k] = da * wir->hy_buf[0][k-ks]
                                                - dbz[material_idx]    * ( ex[i * exny * exnz + j * exnz +(k+1)] - ex[i * exny * exnz + j * exnz +k]);
            + db2  * ( ez[(i+1) * ezny * eznz + j * eznz +k]);
            wir->hy_buf[0][k-ks] = hy[i * hyny * hynz + j * hynz +k];

        }

        j = wir->js - 1;
        /* Do Back of wire */
        for(int k = ks; k < ke; k++)
        {
            int material_idx = ym[i * hyny * hynz + (j-1) * hynz +k];
            double Ur = constants->Ur[material_idx];

            double da = 1;
            double db2 = (2 * simulation->dt) / (mu0*Ur * simulation->dx *log(simulation->dx / wir->radius) );


            hy[i * hyny * hynz + j * hynz +k] = da * wir->hy_buf[1][k-ks]
                                                - dbz[material_idx] * ( ex[i * exny * exnz + j * exnz +(k+1)] - ex[i * exny * exnz + j * exnz +k])
                                                - db2  * ( ez[i * ezny * eznz + j * eznz +k]);

            wir->hy_buf[1][k-ks] = hy[i * hyny * hynz + j * hynz +k];

        }


    }
}

void updateInductors(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants)
{
    double ca,cbx,cby,cv;

    double* ez = grid->Ez;
    double* hx = grid->Hx;
    double* hy = grid->Hy;

    int ezny = grid->ezny;
    int eznz = grid->eznz;

    int hxny = grid->hxny;
    int hxnz = grid->hxnz;

    int hyny = grid->hyny;
    int hynz = grid->hynz;

    int* xe = grid->xe;
    int* ye = grid->ye;
    int* ze = grid->ze;


    for (list <INDUCTOR>::iterator ind = simulation->list_of_inductors.begin(); ind != simulation->list_of_inductors.end(); ind++)
    {
        int i = ind->i;
        int j = ind->j;
        int k = ind->k;
        double value = ind->value;
        int material_idx = ze[i * ezny * eznz + j * eznz +k];
        double Er = constants->Er[material_idx];
        double eCond = 0;

        ca = (2 * eps0*Er - simulation->dt * eCond)/(2 * eps0*Er - simulation->dt * eCond);
        cbx = ( (2 * simulation->dt) / (2 * eps0*Er + simulation->dt * eCond) ) * simulation->dx;
        cby = -( (2 * simulation->dt) / (2 * eps0*Er + simulation->dt * eCond) ) * simulation->dy;
        cv =  -( (2 * simulation->dt) / (2 * eps0*Er + simulation->dt * eCond) );

        ind->Jz = ind->Jz + ( (simulation->dt * simulation->dz) / (value * simulation->dx * simulation->dy) ) * ind->Ez_Old;

        ez[eznz * i * ezny + j * eznz +k] = ca * ind->Ez_Old
                                            + cbx * ((hy[i * hyny * hynz + j * hynz +k]
                                                    - hy[(i-1) * hyny * hynz + j * hynz +k]) / simulation->dx)
                                            - cby * ((hx[i * hxny * hxnz + j * hxnz +k]
                                                    - hx[i * hxny * hxnz + (j-1) * hxnz +k]) / simulation->dy) + cv * ind->Jz;

        ind->Ez_Old = ez[eznz * i * ezny + j * eznz +k];

    }
}

void updateCapacitors(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants)
{
    double ca,cbx,cby;

    double* ez = grid->Ez;
    double* hx = grid->Hx;
    double* hy = grid->Hy;

    int ezny = grid->ezny;
    int eznz = grid->eznz;

    int hxny = grid->hxny;
    int hxnz = grid->hxnz;

    int hyny = grid->hyny;
    int hynz = grid->hynz;

    int* xe = grid->xe;
    int* ye = grid->ye;
    int* ze = grid->ze;


    for (list <CAPACITOR>::iterator cap = simulation->list_of_capacitors.begin(); cap != simulation->list_of_capacitors.end(); cap++)
    {
        int i = cap->i;
        int j = cap->j;
        int k = cap->k;
        double value = cap->value;

        int material_idx = ze[i * ezny * eznz + j * eznz +k];
        double Er = constants->Er[material_idx];

        ca = (2 * Er*eps0  + ((2 * value * simulation->dz)/(simulation->dx * simulation->dy))) /
             (2 * Er*eps0  + ((2 * value * simulation->dz)/(simulation->dx * simulation->dy)));

        cbx = (2 * simulation->dt) / ( (2 * Er*eps0 + ((2 * value * simulation->dz)/(simulation->dx * simulation->dy)) * simulation->dx));
        cby = (2 * simulation->dt) / ( (2 * Er*eps0 + ((2 * value * simulation->dz)/(simulation->dx * simulation->dy) ) * simulation->dy));


        ez[eznz * i * ezny + j * eznz +k] = ca * cap->Ez_Old
                                            + cbx * ((hy[i * hyny * hynz + j * hynz +k]
                                                    - hy[(i-1) * hyny * hynz + j * hynz +k]) / simulation->dx)
                                            - cby * ((hx[i * hxny * hxnz + j * hxnz +k]
                                                    - hx[i * hxny * hxnz + (j-1) * hxnz +k]) / simulation->dy);

        cap->Ez_Old = ez[i * ezny * eznz + j * eznz +k];

    }
}
void updateResistors(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants)
{
    double ca,cbx,cby;

    double* ez = grid->Ez;
    double* hx = grid->Hx;
    double* hy = grid->Hy;

    int* xe = grid->xe;
    int* ye = grid->ye;
    int* ze = grid->ze;


    int ezny = grid->ezny;
    int eznz = grid->eznz;

    int hxny = grid->hxny;
    int hxnz = grid->hxnz;

    int hyny = grid->hyny;
    int hynz = grid->hynz;

    for (list <RESISTOR>::iterator res = simulation->list_of_resistors.begin(); res != simulation->list_of_resistors.end(); res++)
    {
        int i = res->i;
        int j = res->j;
        int k = res->k;
        double value = res->value;



        int material_idx = ze[i * ezny * eznz + j * eznz +k];
        double Er = constants->Er[material_idx];

        ca = (1 - (simulation->dt * simulation->dz) / (2*value*eps0*simulation->dx*simulation->dy))/
             (1 + (simulation->dt * simulation->dz) / (2*value*eps0*simulation->dx*simulation->dy));

        cbx = (simulation->dt/Er*eps0)/(1+(simulation->dt*simulation->dz) / (2*value*Er*eps0*simulation->dx*simulation->dy));
        cby = (simulation->dt/Er*eps0)/(1+(simulation->dt*simulation->dz) / (2*value*Er*eps0*simulation->dx*simulation->dy));

        ez[eznz * i * ezny + j * eznz +k] = ca * res->Ez_Old
                                            + cbx * ((hy[i * hyny * hynz + j * hynz +k]
                                                    - hy[(i-1) * hyny * hynz + j * hynz +k]) / simulation->dx)
                                            - cby * ((hx[i * hxny * hxnz + j * hxnz +k]
                                                    - hx[i * hxny * hxnz + (j-1) * hxnz +k]) / simulation->dy);

        res->Ez_Old = ez[i * ezny * eznz + j * eznz +k];

    }
}





void Ex_Update(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants)
{

    int i,j,k,material_idx;

    int exnx = grid->exnx;
    int exny = grid->exny;
    int exnz = grid->exnz;

    int hznx = grid->hznx;
    int hzny = grid->hzny;
    int hznz = grid->hznz;

    int hynx = grid->hynx;
    int hyny = grid->hyny;
    int hynz = grid->hynz;

    int Nx = simulation->nx;
    int Ny = simulation->ny;
    int Nz = simulation->nz;

    double* ex = grid->Ex;
    double* hy = grid->Hy;
    double* hz = grid->Hz;
    double* ca = constants->ca;
    double* cby = constants->cby;
    double* cbz = constants->cbz;

    for (i=0; i<Nx ; i++)
    {
        for (j=0; j<Ny ; j++)
        {
            for (k=0; k<Nz ; k++)
            {
                material_idx = grid->xe[i * exny * exnz + j * exnz +k];

                ex[i * exny * exnz + j * exnz +k] = ca[material_idx] * ex[i * exny * exnz + j * exnz +k]
                                                    + cby[material_idx]  * ( hz[i * hzny * hznz + j * hznz +k] - hz[i * hzny * hznz + (j-1) * hznz +k])
                                                    - cbz[material_idx]  * ( hy[i * hyny * hynz + j * hynz +k] - hy[i * hyny * hynz + j * hynz + (k-1)]);
            }
        }
    }
}

void Ey_Update(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants)
{

    int i,j,k,material_idx;
    int eynx = grid->eynx;
    int eyny = grid->eyny;
    int eynz = grid->eynz;

    int hzny = grid->hzny;
    int hznz = grid->hznz;

    int hxny = grid->hxny;
    int hxnz = grid->hxnz;

    int Nx = simulation->nx;
    int Ny = simulation->ny;
    int Nz = simulation->nz;

    double* ey = grid->Ey;
    double* hx = grid->Hx;
    double* hz = grid->Hz;
    double* ca = constants->ca;
    double* cbx = constants->cbx;
    double* cbz = constants->cbz;

    for (i=0; i<Nx ; i++)
    {
        for (j=0; j<Ny ; j++)
        {
            for (k=0; k<Nz ; k++)
            {
                material_idx = grid->ye[i * eyny * eynz + j * eynz +k];

                ey[i * eyny * eynz + j * eynz +k] = ca[material_idx] * ey[i * eyny * eynz + j * eynz +k]
                                                    + cbz[material_idx]  * ( hx[i * hxny * hxnz + j * hxnz +k] - hx[i * hxny * hxnz + j * hxnz + (k-1)])
                                                    - cbx[material_idx]  * ( hz[i * hzny * hznz + j * hznz +k] - hz[(i-1) * hzny * hznz + j * hznz +k]) ;
            }
        }
    }

}

void Ez_Update(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants)
{
    int i,j,k,material_idx;

    int eznx = grid->eznx;
    int ezny = grid->ezny;
    int eznz = grid->eznz;

    int hxny = grid->hxny;
    int hxnz = grid->hxnz;

    int hyny = grid->hyny;
    int hynz = grid->hynz;

    double* ez = grid->Ez;
    double* hx = grid->Hx;
    double* hy = grid->Hy;
    double* ca = constants->ca;
    double* cbx = constants->cbx;
    double* cby = constants->cby;

    int Nx = simulation->nx;
    int Ny = simulation->ny;
    int Nz = simulation->nz;

    for (i=0; i<Nx ; i++)
    {
        for (j=0; j<Ny ; j++)
        {
            for (k=0; k<Nz ; k++)
            {

                material_idx = grid->ze[i * ezny * eznz + j * eznz +k];
                ez[i * ezny * eznz + j * eznz +k] = ca[material_idx] * ez[i * ezny * eznz + j * eznz +k]
                                                    + cbx[material_idx]  * ( hy[i * hyny * hynz + j * hynz +k] - hy[(i-1) * hyny * hynz + j * hynz +k])
                                                    - cby[material_idx]  * ( hx[i * hxny * hxnz + j * hxnz +k] - hx[i * hxny * hxnz + (j-1) * hxnz +k]);
            }
        }
    }

}

void Hx_Update(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants)
{

    int i,j,k,material_idx;

    int hxnx = grid->hxnx;
    int hxny = grid->hxny;
    int hxnz = grid->hxnz;

    int eyny = grid->eyny;
    int eynz = grid->eynz;

    int ezny = grid->ezny;
    int eznz = grid->eznz;


    double* hx = grid->Hx;
    double* ey = grid->Ey;
    double* ez = grid->Ez;
    double* da = constants->da;
    double* dby = constants->dby;
    double* dbz = constants->dbz;
    int Nx = simulation->nx;
    int Ny = simulation->ny;
    int Nz = simulation->nz;



    for (i=0; i<Nx ; i++)
    {
        for (j=0; j<Ny ; j++)
        {
            for (k=0; k<Nz ; k++)
            {
                material_idx = grid->xm[i * hxny * hxnz + j * hxnz +k];

                hx[i * hxny * hxnz + j * hxnz +k] = da[material_idx] * hx[i * hxny * hxnz + j * hxnz +k]
                                                    + dbz[material_idx]  * ( ey[i * eyny * eynz + j * eynz +(k+1)] - ey[i * eyny * eynz + j * eynz +k])
                                                    - dby[material_idx]  * ( ez[i * ezny * eznz + (j+1) * eznz +k] - ez[i * ezny * eznz + j * eznz +k]);
            }
        }
    }

}

void Hy_Update(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants)
{
    int i,j,k,material_idx;

    int hynx = grid->hynx;
    int hyny = grid->hyny;
    int hynz = grid->hynz;

    int exnx = grid->exnx;
    int exny = grid->exny;
    int exnz = grid->exnz;

    int eznx = grid->eznx;
    int ezny = grid->ezny;
    int eznz = grid->eznz;

    double* hy = grid->Hy;
    double* ex = grid->Ex;
    double* ez = grid->Ez;
    double* da = constants->da;
    double* dbx = constants->dbx;
    double* dbz = constants->dbz;

    int Nx = simulation->nx;
    int Ny = simulation->ny;
    int Nz = simulation->nz;


    for (i=0; i<Nx ; i++)
    {
        for (j=0; j<Ny ; j++)
        {
            for (k=0; k<Nz ; k++)
            {
                material_idx = grid->ym[i * hyny * hynz + j * hynz +k];
                hy[i * hyny * hynz + j * hynz +k] = da[material_idx] * hy[i * hyny * hynz + j * hynz +k]
                                                    + dbx[material_idx]  * ( ez[(i+1) * ezny * eznz + j * eznz +k] - ez[i * ezny * eznz + j * eznz +k])
                                                    - dbz[material_idx]  * ( ex[i * exny * exnz + j * exnz +(k+1)] - ex[i * exny * exnz + j * exnz +k]);
            }
        }
    }

}

void Hz_Update(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants)
{

    int i,j,k,material_idx;


    int hznx = grid->hznx;
    int hzny = grid->hzny;
    int hznz = grid->hznz;

    int exnx = grid->exnx;
    int exny = grid->exny;
    int exnz = grid->exnz;

    int eynx = grid->eynx;
    int eyny = grid->eyny;
    int eynz = grid->eynz;


    double* hz = grid->Hz;
    double* ex = grid->Ex;
    double* ey = grid->Ey;
    double* da = constants->da;
    double* dbx = constants->dbx;
    double* dby = constants->dby;

    int Nx = simulation->nx;
    int Ny = simulation->ny;
    int Nz = simulation->nz;


    for (i=0; i<Nx ; i++)
    {
        for (j=0; j<Ny ; j++)
        {
            for (k=0; k<Nz ; k++)
            {
                material_idx = grid->zm[i * hzny * hznz + j * hznz +k];
                hz[i * hzny * hznz + j * hznz +k] = da[material_idx] * hz[i * hzny * hznz + j * hznz +k]
                                                    + dby[material_idx]  * ( ex[i * exny * exnz + (j+1) * exnz +k] - ex[i * exny * exnz + j * exnz +k])
                                                    - dbx[material_idx]  * ( ey[(i+1) * eyny * eynz + j * eynz +k] - ey[i * eyny * eynz + j * eynz +k]);
            }
        }
    }

}


