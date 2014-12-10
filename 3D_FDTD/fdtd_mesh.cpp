#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "fdtd.hpp"
#include "fdtd_mesh.hpp"


void readMaterials(MATERIAL_LIST &material_data) {

    std::ifstream matfile("input/materials.txt");

    while (!matfile.eof()) {
        char line[1000];
        matfile.getline(line, 1000);
        if (matfile.eof()) break;

        std::istringstream istr(line);
        std::string name;
        double Er, Ur, econd, mcond;
        istr >> name;
        if (name == "-") // useful for first line ignore
            continue;
        istr >> Er >> Ur >> econd >> mcond;
        material_data.add_material(material(name, Er, Ur, econd, mcond));

        printf("Material [%i]; Name: %s\n", material_data.mats.back().label, name.c_str());
    }

    matfile.close();
}
void construct_mesh(MATERIAL_LIST &material_data, GRID* grid, SIMULATION* simulation)
{
    //readMaterials(material_data);

    add_material(material_data, "Freespace", 1, 1, 0, 0);
    add_material(material_data, "PEC", 1, 1, 0, 0);
    add_material(material_data, "FR4", 4, 1, 1e-4, 0);

    /* Add New Ports */
    // simulation->addPort(60,60,60,50,1,1e9,0,"Sin");
    simulation->addPort(60,60,60,50,1,40e-12,160e-12,"Gauss");

    int is = 60;
    int js = 60;
    int zs = 60;

    di_voxels(grid,simulation,material_data,60-2,60,60-25,4,0,25,"PEC");
    di_voxels(grid,simulation,material_data,60-2,60,61,4,0,25,"PEC");
}



void di_voxels(GRID* grid, SIMULATION* simulation, MATERIAL_LIST &matl,
               int is, int js, int ks,            /* Start Coordinates      */
               int iw, int jw, int kw,            /* Widths in all axes     */
               const char *name)
{

    int i, j, k;

    int exnx = grid->exnx;
    int exny = grid->exny;
    int exnz = grid->exnz;

    int eynx = grid->eynx;
    int eyny = grid->eyny;
    int eynz = grid->eynz;

    int eznx = grid->eznx;
    int ezny = grid->ezny;
    int eznz = grid->eznz;

    int mlabel = get_mlabel(matl, name);

    if (mlabel  < 0)
    {
        printf("=== Error: Undefined Material %s \n",name);
        exit(-1);
    }
    /* If X width = 0, then build a thin plate in the Y/Z plane */
    if(iw==0)
    {
        for(j=js; j<js+jw; j++)
        {
            for(k=ks; k<ks+kw; k++)
            {
                grid->ye[is * eyny * eynz + j * eynz +k] = mlabel;
                grid->ye[is * eyny * eynz + j * eynz + (k+1)] = mlabel;
                grid->ze[is * ezny * eznz + j * eznz + k] = mlabel;
                grid->ze[is * ezny * eznz + (j+1) * eznz + k] = mlabel;
            }
        }
    }

    /* If Y width = 0, then build a thin plate in the X/Z plane */
    if (jw==0)
    {
        for(i=is; i<is+iw; i++)
        {
            for(k=ks; k<ks+kw; k++)
            {
                grid->xe[i * exny * exnz + js * exnz +k]=mlabel;
                grid->xe[i * exny * exnz + js * exnz +(k+1)]=mlabel;
                grid->ze[i * ezny * eznz + js * eznz +k]=mlabel;
                grid->ze[(i+1) * ezny * eznz + js * eznz +k]=mlabel;
            }
        }
    }

    /* If Z width = 0, then build a thin plate in the X/Y plane */
    if (kw==0)
    {
        for(i=is; i<is+iw; i++)
        {
            for(j=js; j<js+jw; j++)
            {
                grid->xe[i * exny * exnz + j * exnz +ks]=mlabel;
                grid->xe[i * exny * exnz + (j+1) * exnz +ks]=mlabel;
                grid->ye[i * eyny * eynz + j * eynz +ks]=mlabel;
                grid->ye[(i+1) * eyny * eynz + j * eynz +ks]=mlabel;

            }
        }
    }
    /* end of thin plate loops */

    for (i=is; i<is+iw; i++)
    {
        for (j=js; j<js+jw; j++)
        {
            for (k=ks; k<ks+kw; k++)
            {
                grid->xe[i * exny * exnz + j * exnz +k]=mlabel;
                grid->xe[i * exny * exnz + j * exnz +(k+1)]=mlabel;
                grid->xe[(i+1) * exny * exnz + (j+1) * exnz +k]=mlabel;
                grid->xe[i * exny * exnz + (j+1) * exnz +k]=mlabel;

                grid->ye[i * eyny * eynz + j * eynz +k]=mlabel;
                grid->ye[(i+1) * eyny * eynz + j * eynz +k]=mlabel;
                grid->ye[(i+1) * eyny * eynz + j * eynz +(k+1)]=mlabel;
                grid->ye[i * eyny * eynz + j * eynz +(k+1)]=mlabel;

                grid->ze[i * ezny * eznz + j * eznz +k]=mlabel;
                grid->ze[(i+1) * ezny * eznz + j * eznz +k]=mlabel;
                grid->ze[(i+1) * ezny * eznz + (j+1) * eznz +k]=mlabel;
                grid->ze[i * ezny * eznz + (j+1) * eznz +k]=mlabel;
            }
        }
    }
}


int get_mlabel(MATERIAL_LIST &matl, const char *name)
{
    if (matl.labelmap.find(name) != matl.labelmap.end()) // if exists
        return(matl.labelmap[name]);
    // not exists
    return(-1);
}

void add_material (MATERIAL_LIST &matl, const char *name, double Er, double Ur, double econd,double mcond)
{
    MATERIAL newMaterial(name, Er, Ur, econd, mcond);
    matl.add_material(newMaterial);

    printf("Material [%i]; Name: %s\n", matl.mats.back().label, name);
}


