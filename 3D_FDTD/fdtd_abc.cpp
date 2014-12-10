#include <stdio.h>
#include <stdlib.h>
#include "fdtd.hpp"
#include <cstring>
#include "fdtd_abc.hpp"



void record_previous_fields(GRID* grid, SIMULATION* simulation,PREVIOUS_EFIELDS* previous_Efields)
{
    int Nx = simulation->nx;
    int Ny = simulation->ny;
    int Nz = simulation->nz;

    int exnx = grid->exnx;
    int exny = grid->exny;
    int exnz = grid->exnz;

    int eynx = grid->eynx;
    int eyny = grid->eyny;
    int eynz = grid->eynz;

    int eznx = grid->eznx;
    int ezny = grid->ezny;
    int eznz = grid->eznz;

    double* ex = grid->Ex;
    double* ey = grid->Ey;
    double* ez = grid->Ez;


    int simSize = exnx*exnz*exnz * sizeof(double);
    // Rotate buffers via the pointers
    double *Extemp = previous_Efields->Ex_old_old_old;
    previous_Efields->Ex_old_old_old = previous_Efields->Ex_old_old;
    previous_Efields->Ex_old_old = previous_Efields->Ex_old;
    previous_Efields->Ex_old = Extemp;
    memcpy(previous_Efields->Ex_old, ex, simSize);

    double *Eytemp = previous_Efields->Ey_old_old_old;
    previous_Efields->Ey_old_old_old = previous_Efields->Ey_old_old;
    previous_Efields->Ey_old_old = previous_Efields->Ey_old;
    previous_Efields->Ey_old = Eytemp;
    memcpy(previous_Efields->Ey_old, ey, simSize);

    double *Eztemp = previous_Efields->Ez_old_old_old;
    previous_Efields->Ez_old_old_old = previous_Efields->Ez_old_old;
    previous_Efields->Ez_old_old = previous_Efields->Ez_old;
    previous_Efields->Ez_old = Eztemp;
    memcpy(previous_Efields->Ez_old, ez, simSize);

}

void abc_liao(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants,PREVIOUS_EFIELDS* previous_Efields)
{

    double* ex = grid->Ex;
    double* ey = grid->Ey;
    double* ez = grid->Ez;

    double* exn = previous_Efields->Ex_old;
    double* eyn = previous_Efields->Ey_old;
    double* ezn = previous_Efields->Ez_old;

    double* exn1 = previous_Efields->Ex_old_old;
    double* eyn1 = previous_Efields->Ey_old_old;
    double* ezn1 = previous_Efields->Ez_old_old;

    double u1,u2,u3,du1,du2,ddu1;

    int Nx = simulation->nx;
    int Ny = simulation->ny;
    int Nz = simulation->nz;

    int exnx = grid->exnx;
    int exny = grid->exny;
    int exnz = grid->exnz;

    int eynx = grid->eynx;
    int eyny = grid->eyny;
    int eynz = grid->eynz;

    int eznx = grid->eznx;
    int ezny = grid->ezny;
    int eznz = grid->eznz;

    /* Y - Plate */
    for (int j = 0 ; j <= Nx; j++)
    {
        for (int k = 0 ; k <= Ny; k++)
        {

            /* Left Side */

            if  (j < Nx) {
                ey[0 * eyny * eynz + j * eynz +k] = 2 * eyn[1 * eyny * eynz + j * eynz +k] - eyn1[2 * eyny * eynz + j * eynz +k];
            }

            if (k < Nx) {
                ez[0 * ezny * eznz + j * eznz +k] = 2 * ezn[1 * ezny * eznz + j * eznz +k] - ezn1[2 * ezny * eznz + j * eznz +k];
            }

            /* Right Side */

            if (j < Nx) {
                ey[Nx * eyny * eynz + j * eynz +k] = 2 * eyn[(Nx-1) * eyny * eynz + j * eynz +k] - eyn1[(Nx-2) * eyny * eynz + j * eynz +k];
            }

            if (k < Nx) {
                ez[Nx * ezny * eznz + j * eznz +k] = 2 * ezn[(Nx-1) * ezny * eznz + j * eznz +k] - ezn1[(Nx-2) * ezny * eznz + j * eznz +k];
            }
        }
    }

    /* X-Plate */
    for (int i = 0 ; i <= Nx; i++)
    {
        for (int k = 0 ; k <= Ny; k++)
        {

            /* Front Side */
            if (i < Nx) {
                ex[i * exny * exnz + (Ny) * exnz +k] = 2 * exn[i * exny * exnz + (Ny-1) * exnz +k] - exn1[i * exny * exnz + (Ny-2) * exnz +k];
            }



            if (k < Nz) {
                ez[i * ezny * eznz + (Ny) * eznz +k] = 2 * ezn[i * ezny * eznz + (Ny-1) * eznz +k] - ezn1[i * ezny * eznz + (Ny-2) * eznz +k];
            }


            /* Back Side */
            if (i < Nx) {
                ex[i * exny * exnz + 0 * exnz +k] = 2 * exn[i * exny * exnz + 1 * exnz +k] - exn1[i * exny * exnz + 2 * exnz +k];
            }



            if (k < Nz) {
                ez[i * ezny * eznz + 0 * eznz +k] = 2 * ezn[i * ezny * eznz + 1 * eznz +k] - ezn1[i * ezny * eznz + 2 * eznz +k];
            }
        }
    }

    /* Z-Plate */
    for (int i = 0 ; i <= Nx; i++)
    {
        for (int j = 0 ; j <= Ny; j++)
        {

            /* Top Size */
            if (i < Nx) {
                ex[i * exny * exnz + j * exnz + (Nz)] = 2 * exn[i * exny * exnz + j * exnz + (Nz-1)] - exn1[i * exny * exnz + j * exnz + (Nz-2)];
            }

            if (j < Ny) {
                ey[i * eyny * eynz + j * eynz + (Nz)] = 2 * eyn[i * eyny * eynz + j * eynz + (Nz-1)] - eyn1[i * eyny * eynz + j * eynz + (Nz-2)];
            }


            /* Bottom Size */
            if (i < Nx ) {
                ex[i * exny * exnz + j * exnz + 0] = 2 * exn[i * exny * exnz + j * exnz + 1] - exn1[i * exny * exnz + j * exnz + 2];
            }

            if (j < Ny) {
                ey[i * eyny * eynz + j * eynz + 0] = 2 * eyn[i * eyny * eynz + j * eynz + 1] - eyn1[i * eyny * eynz + j * eynz + 2];
            }
        }
    }


}


void abc_murfirst(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants,PREVIOUS_EFIELDS* previous_Efields)
{
    double* ex = grid->Ex;
    double* ey = grid->Ey;
    double* ez = grid->Ez;

    int Nx = simulation->nx;
    int Ny = simulation->ny;
    int Nz = simulation->nz;

    int j,k;
    int exnx = grid->exnx;
    int exny = grid->exny;
    int exnz = grid->exnz;

    int eynx = grid->eynx;
    int eyny = grid->eyny;
    int eynz = grid->eynz;

    int eznx = grid->eznx;
    int ezny = grid->ezny;
    int eznz = grid->eznz;


    double k1 = (c * simulation->dt - simulation->dx) / (c * simulation->dt + simulation->dx);
    for (j = 0 ; j <= simulation->nx; j++)
    {
        for (k = 0 ; k <= simulation->ny; k++)
        {


            /* Left Side */
            if  (j < 120) {
                ey[0 * eyny * eynz + j * eynz + k] = previous_Efields->Ey_old[1 * eyny * eynz + j * eynz + k] + k1 * (ey[1 * eyny * eynz + j * eynz + k] - previous_Efields->Ey_old[0 * eyny * eynz + j * eynz + k]);
            }

            if (k < 120) {
                ez[0 * ezny * eznz + j * eznz + k] = previous_Efields->Ez_old[1 * ezny * eznz + j * eznz + k] + k1 * (ez[1 * ezny * eznz + j * eznz + k] - previous_Efields->Ez_old[0 * ezny * eznz + j * eznz + k]);
            }

            /* Right Side */
            if (j < 120) {
                ey[Nx * eyny * eynz + j * eynz + k] = previous_Efields->Ey_old[Nx-1 * eyny * eynz + j * eynz + k] + k1 * (ey[Nx-1 * eyny * eynz + j * eynz + k] - previous_Efields->Ey_old[Nx * eyny * eynz + j * eynz + k]);
            }

            if (k < 120) {
                ez[Nx * ezny * eznz + j * eznz + k] = previous_Efields->Ez_old[Nx-1 * ezny * eznz + j * eznz + k] + k1 * (ez[Nx-1 * ezny * eznz + j * eznz + k] - previous_Efields->Ez_old[Nx * ezny * eznz + j * eznz + k]);
            }



            /* Front Side */
            if (j < 120) {
                ex[j * exny * exnz + Ny * exnz + k] = previous_Efields->Ex_old[j * exny * exnz + (Nx-1) * exnz + k] + k1 * (ex[j * exny * exnz + (Nx-1) * exnz + k] - previous_Efields->Ex_old[j * exny * exnz + Nx * exnz + k]);
            }


            if (k < 120) {
                ez[j * ezny * eznz + Ny * eznz + k] = previous_Efields->Ez_old[j * ezny * eznz + (Nx-1) * eznz + k] + k1 * (ez[j * ezny * eznz + (Nx-1) * eznz + k] - previous_Efields->Ez_old[j * ezny * eznz + Nx * eznz + k]);
            }


            /* Back Side */
            if (j < 120)  {
                ex[j * exny * exnz + 0 * exnz + k] = previous_Efields->Ex_old[j * exny * exnz + 1 * exnz + k] + k1 * (ex[j * exny * exnz + 1 * exnz + k] - previous_Efields->Ex_old[j * exny * exnz + 0 * exnz + k]);
            }


            if (k < 120)  {
                ez[j * ezny * eznz + 0 * eznz + k] = previous_Efields->Ez_old[j * ezny * eznz + 1 * eznz + k] + k1 * (ez[j * ezny * eznz + 1 * eznz + k] - previous_Efields->Ez_old[j * ezny * eznz + 0 * eznz + k]);
            }

            /* Top Size */
            if (j < 120) {
                ex[j * exny * exnz + k * exnz + Nz] = previous_Efields->Ex_old[j * exny * exnz + k * exnz + (Nz-1)] + k1 * (ex[j * exny * exnz + k * exnz + (Nx-1)] - previous_Efields->Ex_old[j * exny * exnz + k * exnz + Nz]);
            }

            if (k < 120) {
                ey[j * eyny * eynz + k * eynz + Nz] = previous_Efields->Ey_old[j * eyny * eynz + k * eynz + (Nz-1)] + k1 * (ey[j * eyny * eynz + k * eynz + (Nx-1)] - previous_Efields->Ey_old[j * eyny * eynz + k * eynz + Nz]);
            }


            /* Bottom Size */
            if (j < 120 ) {
                ex[j * exny * exnz + k * exnz + 0] = previous_Efields->Ex_old[j * exny * exnz + k * exnz + 1] + k1 * (ex[j * exny * exnz + k * exnz + 1] - previous_Efields->Ex_old[j * exny * exnz + k * exnz + 0]);
            }

            if (k < 120) {
                ey[j * eyny * eynz + k * eynz + 0] = previous_Efields->Ey_old[j * eyny * eynz + k * eynz + 1] + k1 * (ey[j * eyny * eynz + k * eynz + 1] - previous_Efields->Ey_old[j * eyny * eynz + k * eynz + 0]);
            }

        }

    }

}





