#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <list>
#include "fdtd.hpp"
#include "fdtd_sparameters.hpp"
#include "fdtd_farfield.hpp"
#include "fdtd_abc.hpp"
#include "fdtd_mesh.hpp"
#include "fdtd_updates.hpp"

/* ---------------------------------<BEGIN> Function Declarations --------------------------------------*/
void initilise_fdconstants(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants,MATERIAL_LIST &list);
void allocate_memory_fdconstants(FD_CONSTANTS* constants,GRID* grid, MATERIAL_LIST &matl);
void allocate_memory_previousEfields(GRID* grid,SIMULATION* simulation,PREVIOUS_EFIELDS* previous_Efields);
void allocate_memory(GRID* grid,SIMULATION* simulation,FD_CONSTANTS* constants);
/* ---------------------------------<END> Function Declarations --------------------------------------*/

int main (void)
{
    int f;
    /* Time-Domain Parameters */
    GRID grid;
    SIMULATION simulation;
    FD_CONSTANTS fd_constants;
    PREVIOUS_EFIELDS previous_Efields;
    MATERIAL_LIST list_of_materials;

    /* Frequency Domain Parameters */
    RADIATION_POWER radiation_power;

    /* Mesh Size and Resolution */
    simulation.set_mesh_size(120,120,120);
    simulation.set_mesh_resolution(1e-3,1e-3,1e-3);
    simulation.set_simulation_length(5000);

    /* A List of Frequencies Must Be Defined First */
    for (double fr = 1e9; fr < 10e9; fr+=50e6) {
	    simulation.addFrequency(fr);
    }
    simulation.initClosedSurface();


    /* Intitilise list of dielectric/magnetic materials */
    // No special initialization needed.

    /* Allocate Memory fo the FD equations */
    allocate_memory(&grid,&simulation,&fd_constants);
    allocate_memory_previousEfields(&grid,&simulation,&previous_Efields);

    /* Construct the Mesh of Voxels */
    construct_mesh(list_of_materials, &grid, &simulation);

    /* Set the FD Constants */
    allocate_memory_fdconstants(&fd_constants, &grid, list_of_materials);
    initilise_fdconstants(&grid,&simulation, &fd_constants, list_of_materials);


    printf("=== Time Constant....\n");
    printf("\t dT %e \n",simulation.dt);

    int nstop = simulation.Nt;
    char data[100];
    for (int n = 1; n <= nstop ; n++)
    {
        sprintf(data,"Time Index: [%i / %i]",n,nstop);
        cout << data << "\r" ;
        fflush(stdout);

        simulation.n = n;

        Ex_Update(&grid,&simulation,&fd_constants);
        Ey_Update(&grid,&simulation,&fd_constants);
        Ez_Update(&grid,&simulation,&fd_constants);
        abc_liao(&grid,&simulation,&fd_constants,&previous_Efields);
        record_previous_fields(&grid, &simulation,&previous_Efields);

        simulation.updateClosedSurfaces(&grid);

        updateExcitation(&grid, &simulation);

        updateResistors(&grid,&simulation,&fd_constants);
        updateCapacitors(&grid,&simulation,&fd_constants);
        updateInductors(&grid,&simulation,&fd_constants);

        Hx_Update(&grid,&simulation,&fd_constants);
        Hy_Update(&grid,&simulation,&fd_constants);
        Hz_Update(&grid,&simulation,&fd_constants);

        updateWires(&grid,&simulation,&fd_constants);

    }
    computeSParameters(&simulation);

    cout << endl;
    cout << endl;
    cout << endl;

    FILE* fp = fopen("./radiation/efficiency.fd","w");
    int index = 0;
    for (list<CLOSED_SURFACE>::iterator it = simulation.list_of_closed_surfaces.begin(); it != simulation.list_of_closed_surfaces.end(); it++)
    {
        list <PORT>::iterator port  = simulation.list_of_ports.begin();

        cout << "**** Computing Near-Far Field at: "  << it->frequency / 1e9 << "GHz *****" << endl;
        it->compute_surface_currents();
        it->compute_radiated_power(port,simulation.dx,simulation.dy,simulation.dz);
        it->compute_farfield_xyplane(simulation.dx,simulation.dy,simulation.dz,index);
        it->compute_farfield_xzplane(simulation.dx,simulation.dy,simulation.dz,index);
 	fprintf(fp,"%lg %lg %lg %lg\n",it->frequency, it->pin, it->prad, it->efficiency*100);
        cout << "Radiated Power: " << it->prad << " Watts" << endl;
        cout << "Port Power: " << it->pin << " Watts" << endl;
        cout << "Efficiency: " << it->efficiency*100 << "%" << endl << endl;
        index++;
    }
    fclose(fp);
}

double getGaussian(double t,double amplitude,double width,double t0)
{
    return amplitude * exp ( (-1.0 * pow((t - t0),2.0)) / (width * width) );
}

double getSinusoidal(double t,double ampltiude,double frequency,double tmp)
{
    return ampltiude * sin(2 * pi * frequency * t);
}

double getNull(double t,double val1, double val2,double val3)
{
    return 0;
}



void SIMULATION::addPort(int i, int j, int k, double zo, double val1, double val2, double val3, string Type)
{
    PORT newPort;
    newPort.i = i;
    newPort.j = j;
    newPort.k = k;
    newPort.zo = zo;
    newPort.parameters[0] = val1;
    newPort.parameters[1] = val2;
    newPort.parameters[2] = val3;

    if (Type.compare("Sin") == 0 )
        newPort.signal = &getSinusoidal;
    else if (Type.compare("Gauss") == 0)
        newPort.signal = &getGaussian;
    else if (Type.compare("Null") == 0)
        newPort.signal = &getNull;

    else {
        cout << "No Port of that Type Exists !! "<< endl;
        exit(-1);
    }
    newPort.v_in = (double*)calloc(Nt,sizeof(double));
    newPort.i_in = (double*)calloc(Nt,sizeof(double));
    newPort.V_in = (complx*)malloc(getFrequencyLength() * sizeof(complx));
    newPort.I_in = (complx*)malloc(getFrequencyLength() * sizeof(complx));

    for (int i = 0; i < getFrequencyLength(); i++)
    {
        newPort.V_in[i] = (0,0);
        newPort.I_in[i] = (0,0);
    }

    list_of_ports.push_back(newPort);
}


void SIMULATION::set_simulation_length(int N)
{
    Nt = N;
}


void SIMULATION::set_mesh_size(int Nx,int Ny, int Nz)
{
    nx = Nx;
    ny = Ny;
    nz = Nz;
}
void SIMULATION::set_mesh_resolution(double delx, double dely, double delz)
{
    dx = delx;
    dy = dely;
    dz = delz;
    dt = 1.0 / (c * sqrt( ((1.0/(dx*dx)) + (1.0/(dy*dy)) +(1.0/(dz*dz)))));
}

int SIMULATION::getFrequencyLength()
{
    return frequencies.size();
}

void SIMULATION::addResistor(int i,int j,int k,double value)
{
    RESISTOR newResistor;
    newResistor.i = i;
    newResistor.j = j;
    newResistor.k = k;
    newResistor.value = value;
    list_of_resistors.push_back(newResistor);
}

void SIMULATION::addCapacitor(int i,int j,int k,double value)
{
    CAPACITOR newCapacitor;
    newCapacitor.i = i;
    newCapacitor.j = j;
    newCapacitor.k = k;
    newCapacitor.value = value;
    list_of_capacitors.push_back(newCapacitor);
}

void SIMULATION::addInductor(int i,int j,int k,double value)
{
    INDUCTOR newInductor;
    newInductor.i = i;
    newInductor.j = j;
    newInductor.k = k;
    newInductor.value = value;
    list_of_inductors.push_back(newInductor);
}

void SIMULATION::addWire(int i,int j,int k,int len,double rad)
{
    WIRE newWire;
    newWire.is = i;
    newWire.js = j;
    newWire.ks = k;
    newWire.length = len;
    newWire.radius = rad;
    list_of_wires.push_back(newWire);
}


void SIMULATION::addFrequency(double f)
{
    frequencies.push_back(f);
}



void SIMULATION::initClosedSurface()
{
    int index = 0;
    for (list<double>::iterator it=frequencies.begin(); it != frequencies.end(); it++)
    {
        CLOSED_SURFACE newClosedSurface;
        newClosedSurface.frequency = *it;
        newClosedSurface.allocate_surface_memory(nx,ny,nz,index);
        list_of_closed_surfaces.push_back(newClosedSurface);
        index++;
    }
}


void SIMULATION::updateClosedSurfaces(GRID* grid)
{
    for (list<CLOSED_SURFACE>::iterator it=list_of_closed_surfaces.begin(); it != list_of_closed_surfaces.end(); it++)
        it->perform_DFT(grid,n,dt);
}

void initilise_fdconstants(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants,MATERIAL_LIST &list)
{
    double top1, top2, bot1, bot2;
    int NoMaterials = list.mats.size();
    int i;
    int exnz = grid->exnz;
    int eynz = grid->eynz;
    int eznz = grid->eznz;

    int hxnz = grid->hxnz;
    int hynz = grid->hynz;
    int hznz = grid->hxnz;

    for (int n = 0; n < NoMaterials; n++)
    {
        MATERIAL *mat = &list.mats[n];

        constants->Er[n] = mat->Er;
        constants->Ur[n] = mat->Ur;


        /* Electric Field Updates Constants*/
        top1 = mat->econd * simulation->dt;
        top2 = 2.0 * eps0 * mat->Er;
        bot1 = mat->econd * simulation->dt;
        bot2 = 2.0 * eps0 * mat->Er;
        constants->ca[n] = (1.0 - top1 / top2) / (1.0 + bot1 / bot2);

        top1 = simulation->dt;
        top2 = mat->Er * eps0 * simulation->dx;
        bot1 = mat->econd * simulation->dt;
        bot2 = 2.0 * eps0 * mat->Er;
        constants->cbx[n] = (top1 / top2) / (1.0 + bot1 / bot2);

        top1 = simulation->dt;
        top2 = mat->Er * eps0 * simulation->dy;
        bot1 = mat->econd * simulation->dt;
        bot2 = 2.0 * eps0 * mat->Er;
        constants->cby[n] = (top1 / top2) / (1.0 + bot1 / bot2);

        top1 = simulation->dt;
        top2 = mat->Er * eps0 * simulation->dz;
        bot1 = mat->econd * simulation->dt;
        bot2 = 2.0 * eps0 * mat->Er;
        constants->cbz[n] = (top1 / top2) / (1.0 + bot1 / bot2);

        /* Magnetic Field Updates Constants*/
        top1 = mat->mcond * simulation->dt;
        top2 = 2.0 * mu0 * mat->Ur ;
        bot1 = mat->mcond * simulation->dt;
        bot2 = 2.0 * mu0 * mat->Ur;
        constants->da[n] = (1.0 - top1 / top2) / (1.0 + bot1 / bot2);

        top1 = simulation->dt;
        top2 = mat->Ur * mu0 * simulation->dx;
        bot1 = mat->mcond * simulation->dt;
        bot2 = 2.0 * mu0 * mat->Ur;
        constants->dbx[n] = (top1 / top2) / (1.0 + bot1 / bot2);

        top1 = simulation->dt;
        top2 = mat->Ur * mu0 * simulation->dy;
        bot1 = mat->mcond * simulation->dt;
        bot2 = 2.0 * mu0 * mat->Ur;
        constants->dby[n] = (top1 / top2) / (1.0 + bot1 / bot2);

        top1 = simulation->dt;
        top2 = mat->Ur * mu0 * simulation->dz;
        bot1 = mat->mcond * simulation->dt;
        bot2 = 2.0 * mu0 * mat->Ur;
        constants->dbz[n] = (top1 / top2) / (1.0 + bot1 / bot2);

        /* PEC */
        if (mat->name == "PEC")
        {
            constants->ca[n] = 0;
            constants->cbx[n] = 0;
            constants->cby[n] = 0;
            constants->cbz[n] = 0;

            constants->da[n] = 0;
            constants->dbx[n] = 0;
            constants->dby[n] = 0;
            constants->dbz[n] = 0;
        }
    }
}

void allocate_memory_previousEfields(GRID* grid,SIMULATION* simulation,PREVIOUS_EFIELDS* previous_Efields)
{
    int exnx = simulation->nx;
    int exny = simulation->ny + 1;
    int exnz = simulation->nz + 1;

    int eynx = simulation->nx + 1;
    int eyny = simulation->ny;
    int eynz = simulation->nz + 1;

    int eznx = simulation->nx + 1;
    int ezny = simulation->ny + 1;
    int eznz = simulation->nz;


    int hxnx = simulation->nx + 1;
    int hxny = simulation->ny + 1;
    int hxnz = simulation->nz + 1;

    int hynx = simulation->nx + 1;
    int hyny = simulation->ny + 1;
    int hynz = simulation->nz + 1;

    int hznx = simulation->nx + 1;
    int hzny = simulation->ny + 1;
    int hznz = simulation->nz + 1;


    previous_Efields->Ex_old = (double*)calloc(exnx * exny * exnz,sizeof(double));
    previous_Efields->Ey_old = (double*)calloc(eynx * eyny * eynz,sizeof(double));
    previous_Efields->Ez_old = (double*)calloc(eznx * ezny * eznz,sizeof(double));

    previous_Efields->Ex_old_old = (double*)calloc(exnx * exny * exnz,sizeof(double));
    previous_Efields->Ey_old_old = (double*)calloc(eynx * eyny * eynz,sizeof(double));
    previous_Efields->Ez_old_old = (double*)calloc(eznx * ezny * eznz,sizeof(double));

    previous_Efields->Ex_old_old_old = (double*)calloc(exnx * exny * exnz,sizeof(double));
    previous_Efields->Ey_old_old_old = (double*)calloc(eynx * eyny * eynz,sizeof(double));
    previous_Efields->Ez_old_old_old = (double*)calloc(eznx * ezny * eznz,sizeof(double));

}



void allocate_memory(GRID* grid,SIMULATION* simulation,FD_CONSTANTS* constants)
{
    int exnx = simulation->nx;
    int exny = simulation->ny + 1;
    int exnz = simulation->nz + 1;

    int eynx = simulation->nx + 1;
    int eyny = simulation->ny;
    int eynz = simulation->nz + 1;

    int eznx = simulation->nx + 1;
    int ezny = simulation->ny + 1;
    int eznz = simulation->nz;


    int hxnx = simulation->nx + 1;
    int hxny = simulation->ny + 1;
    int hxnz = simulation->nz + 1;

    int hynx = simulation->nx + 1;
    int hyny = simulation->ny + 1;
    int hynz = simulation->nz + 1;

    int hznx = simulation->nx + 1;
    int hzny = simulation->ny + 1;
    int hznz = simulation->nz + 1;

    /* E-Field Size */
    grid->exnx = exnx;
    grid->exny = exny;
    grid->exnz = exnz;

    grid->eynx = eynx;
    grid->eyny = eyny;
    grid->eynz = eynz;

    grid->eznx = eznx;
    grid->ezny = ezny;
    grid->eznz = eznz;

    /* H-Field Size */
    grid->hxnx = hxnx;
    grid->hxny = hxny;
    grid->hxnz = hxnz;

    grid->hynx = hynx;
    grid->hyny = hyny;
    grid->hynz = hynz;

    grid->hznx = hznx;
    grid->hzny = hzny;
    grid->hznz = hznz;

    grid->Ex = (double*)calloc(exnx * exny * exnz,sizeof(double));
    grid->Ey = (double*)calloc(eynx * eyny * eynz,sizeof(double));
    grid->Ez = (double*)calloc(eznx * ezny * eznz,sizeof(double));

    grid->Hx = (double*)calloc(hxnx * hxny * hxnz,sizeof(double));
    grid->Hy = (double*)calloc(hynx * hyny * hynz,sizeof(double));
    grid->Hz = (double*)calloc(hznx * hzny * hznz,sizeof(double));

    grid->xe = (int*)calloc(exnx * exny * exnz,sizeof(int));
    grid->ye = (int*)calloc(eynx * eyny * eynz,sizeof(int));
    grid->ze = (int*)calloc(eznx * ezny * eznz,sizeof(int));

    grid->xm = (int*)calloc(hxnx * hxny * hxnz,sizeof(int));
    grid->ym = (int*)calloc(hynx * hyny * hynz,sizeof(int));
    grid->zm = (int*)calloc(hznx * hzny * hznz,sizeof(int));
}

void allocate_memory_fdconstants(FD_CONSTANTS* constants,GRID* grid, MATERIAL_LIST &matl)
{
    int NoMaterials = matl.mats.size();
    int exnz = grid->exnz;
    int eynz = grid->eynz;
    int eznz = grid->eznz;

    int hxnz = grid->hxnz;
    int hynz = grid->hynz;
    int hznz = grid->hxnz;

    constants->Er =  (double*)calloc(NoMaterials,sizeof(double));
    constants->Ur =  (double*)calloc(NoMaterials,sizeof(double));

    constants->ca =  (double*)calloc(NoMaterials,sizeof(double));
    constants->cbx = (double*)calloc(NoMaterials,sizeof(double));
    constants->cby = (double*)calloc(NoMaterials,sizeof(double));
    constants->cbz = (double*)calloc(NoMaterials,sizeof(double));

    constants->da =  (double*)calloc(NoMaterials,sizeof(double));
    constants->dbx = (double*)calloc(NoMaterials,sizeof(double));
    constants->dby = (double*)calloc(NoMaterials,sizeof(double));
    constants->dbz = (double*)calloc(NoMaterials,sizeof(double));

}








