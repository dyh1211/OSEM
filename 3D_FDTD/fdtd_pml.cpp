#include <stdio.h>
#include <stdlib.h>
#include "fdtd.hpp"
#include "fdtd_pml.hpp"

//  Specify the CPML Thickness in Each Direction (Value of Zero
//  Corresponds to No PML, and the Grid is Terminated with a PEC)
// PML thickness in each direction
int nxPML_1, nxPML_2, nyPML_1;
int nyPML_2, nzPML_1, nzPML_2;

//  Specify the CPML Order and Other Parameters:
int m = 3, ma = 1;

double sig_x_max;
double sig_y_max;
double sig_z_max;
double alpha_x_max;
double alpha_y_max;
double alpha_z_max;
double kappa_x_max;
double kappa_y_max;
double kappa_z_max;

//  CPML components (Taflove 3rd Edition, Chapter 7)
/* 3 - D Arrays */
double ***psi_Ezx_1;
double ***psi_Ezx_2;
double ***psi_Hyx_1;
double ***psi_Hyx_2;
double ***psi_Ezy_1;
double ***psi_Ezy_2;
double ***psi_Hxy_1;
double ***psi_Hxy_2;
double ***psi_Hxz_1;
double ***psi_Hxz_2;
double ***psi_Hyz_1;
double ***psi_Hyz_2;
double ***psi_Exz_1;
double ***psi_Exz_2;
double ***psi_Eyz_1;
double ***psi_Eyz_2;
double ***psi_Hzx_1;
double ***psi_Eyx_1;
double ***psi_Hzx_2;
double ***psi_Eyx_2;
double ***psi_Hzy_1;
double ***psi_Exy_1;
double ***psi_Hzy_2;
double ***psi_Exy_2;

double *be_x_1, *ce_x_1, *alphae_x_PML_1, *sige_x_PML_1, *kappae_x_PML_1;
double *bh_x_1, *ch_x_1, *alphah_x_PML_1, *sigh_x_PML_1, *kappah_x_PML_1;
double *be_x_2, *ce_x_2, *alphae_x_PML_2, *sige_x_PML_2, *kappae_x_PML_2;
double *bh_x_2, *ch_x_2, *alphah_x_PML_2, *sigh_x_PML_2, *kappah_x_PML_2;
double *be_y_1, *ce_y_1, *alphae_y_PML_1, *sige_y_PML_1, *kappae_y_PML_1;
double *bh_y_1, *ch_y_1, *alphah_y_PML_1, *sigh_y_PML_1, *kappah_y_PML_1;
double *be_y_2, *ce_y_2, *alphae_y_PML_2, *sige_y_PML_2, *kappae_y_PML_2;
double *bh_y_2, *ch_y_2, *alphah_y_PML_2, *sigh_y_PML_2, *kappah_y_PML_2;
double *be_z_1, *ce_z_1, *alphae_z_PML_1, *sige_z_PML_1, *kappae_z_PML_1;
double *bh_z_1, *ch_z_1, *alphah_z_PML_1, *sigh_z_PML_1, *kappah_z_PML_1;
double *be_z_2, *ce_z_2, *alphae_z_PML_2, *sige_z_PML_2, *kappae_z_PML_2;
double *bh_z_2, *ch_z_2, *alphah_z_PML_2, *sigh_z_PML_2, *kappah_z_PML_2;



void allocate_pml_memory(GRID* grid, SIMULATION* simulation,PML_PARAMETERS* pml_parameters)
{
    //PML Layers (10 layers)
    nxPML_1 = 11;
    nxPML_2 = 11;
    nyPML_1 = 11;
    nyPML_2 = 11;
    nzPML_1 = 11;
    nzPML_2 = 11;
    int i,j,k;
    int Imax = simulation->nx;
    int Jmax = simulation->ny;
    int Kmax = simulation->nz;

    psi_Ezx_1 = (double ***)malloc(nxPML_1 * sizeof(double **));
    for(i = 0; i < nxPML_1; i++) {
        psi_Ezx_1[i] = (double **)malloc(Jmax * sizeof(double *));
        for(j = 0; j < Jmax; j++) {
            psi_Ezx_1[i][j] = (double *)malloc(Kmax * sizeof(double));
            for(k = 0; k < Kmax; k++) {
                psi_Ezx_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Ezx_2 = (double ***)malloc(nxPML_2 * sizeof(double **));
    for(i = 0; i < nxPML_2; i++) {
        psi_Ezx_2[i] = (double **)malloc(Jmax * sizeof(double *));
        for(j = 0; j < Jmax; j++) {
            psi_Ezx_2[i][j] = (double *)malloc(Kmax * sizeof(double));
            for(k = 0; k < Kmax; k++) {
                psi_Ezx_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Hyx_1 = (double ***)malloc((nxPML_1-1) * sizeof(double **));
    for(i = 0; i < nxPML_1-1; i++) {
        psi_Hyx_1[i] = (double **)malloc(Jmax * sizeof(double *));
        for(j = 0; j < Jmax; j++) {
            psi_Hyx_1[i][j] = (double *)malloc(Kmax * sizeof(double));
            for(k = 0; k < Kmax; k++) {
                psi_Hyx_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Hyx_2 = (double ***)malloc((nxPML_2-1) * sizeof(double **));

    for(i = 0; i < nxPML_1-1; i++) {

        psi_Hyx_2[i] = (double **)malloc(Jmax * sizeof(double *));

        for(j = 0; j < Jmax; j++) {

            psi_Hyx_2[i][j] = (double *)malloc(Kmax * sizeof(double));

            for(k = 0; k < Kmax; k++) {

                psi_Hyx_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Ezy_1 = (double ***)malloc(Imax * sizeof(double **));

    for(i = 0; i < Imax; i++) {

        psi_Ezy_1[i] = (double **)malloc(nyPML_1 * sizeof(double *));

        for(j = 0; j < nyPML_1; j++) {

            psi_Ezy_1[i][j] = (double *)malloc(Kmax * sizeof(double));

            for(k = 0; k < Kmax; k++) {

                psi_Ezy_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Ezy_2 = (double ***)malloc(Imax * sizeof(double **));

    for(i = 0; i < Imax; i++) {

        psi_Ezy_2[i] = (double **)malloc(nyPML_2 * sizeof(double *));

        for(j = 0; j < nyPML_2; j++) {

            psi_Ezy_2[i][j] = (double *)malloc(Kmax * sizeof(double));

            for(k = 0; k < Kmax; k++) {

                psi_Ezy_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Hxy_1 = (double ***)malloc(Imax * sizeof(double **));

    for(i = 0; i < Imax; i++) {

        psi_Hxy_1[i] = (double **)malloc((nyPML_1-1) * sizeof(double *));

        for(j = 0; j < nyPML_1-1; j++) {

            psi_Hxy_1[i][j] = (double *)malloc(Kmax * sizeof(double));

            for(k = 0; k < Kmax; k++) {

                psi_Hxy_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Hxy_2 = (double ***)malloc(Imax * sizeof(double **));

    for(i = 0; i < Imax; i++) {

        psi_Hxy_2[i] = (double **)malloc((nyPML_2-1) * sizeof(double *));

        for(j = 0; j < nyPML_2-1; j++) {

            psi_Hxy_2[i][j] = (double *)malloc(Kmax * sizeof(double));

            for(k = 0; k < Kmax; k++) {

                psi_Hxy_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Hxz_1 = (double ***)malloc(Imax * sizeof(double **));

    for(i = 0; i < Imax; i++) {

        psi_Hxz_1[i] = (double **)malloc((Jmax-1) * sizeof(double *));

        for(j = 0; j < Jmax; j++) {

            psi_Hxz_1[i][j] = (double *)malloc((nzPML_1-1) * sizeof(double));

            for(k = 0; k < nzPML_1-1; k++) {

                psi_Hxz_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Hxz_2 = (double ***)malloc(Imax * sizeof(double **));

    for(i = 0; i < Imax; i++) {

        psi_Hxz_2[i] = (double **)malloc((Jmax-1) * sizeof(double *));

        for(j = 0; j < Jmax; j++) {

            psi_Hxz_2[i][j] = (double *)malloc((nzPML_2-1) * sizeof(double));

            for(k = 0; k < nzPML_2-1; k++) {

                psi_Hxz_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Hyz_1 = (double ***)malloc((Imax-1) * sizeof(double **));

    for(i = 0; i < Imax-1; i++) {

        psi_Hyz_1[i] = (double **)malloc(Jmax * sizeof(double *));

        for(j = 0; j < Jmax; j++) {

            psi_Hyz_1[i][j] = (double *)malloc((nzPML_1-1) * sizeof(double));

            for(k = 0; k < nzPML_1-1; k++) {

                psi_Hyz_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Hyz_2 = (double ***)malloc((Imax-1) * sizeof(double **));

    for(i = 0; i < Imax-1; i++) {

        psi_Hyz_2[i] = (double **)malloc(Jmax * sizeof(double *));

        for(j = 0; j < Jmax; j++) {

            psi_Hyz_2[i][j] = (double *)malloc((nzPML_2-1) * sizeof(double));

            for(k = 0; k < nzPML_2-1; k++) {

                psi_Hyz_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Exz_1 = (double ***)malloc((Imax-1) * sizeof(double **));

    for(i = 0; i < Imax-1; i++) {

        psi_Exz_1[i] = (double **)malloc(Jmax * sizeof(double *));

        for(j = 0; j < Jmax; j++) {

            psi_Exz_1[i][j] = (double *)malloc(nzPML_1 * sizeof(double));

            for(k = 0; k < nzPML_1; k++) {

                psi_Exz_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Exz_2 = (double ***)malloc((Imax-1) * sizeof(double **));

    for(i = 0; i < Imax-1; i++) {

        psi_Exz_2[i] = (double **)malloc(Jmax * sizeof(double *));

        for(j = 0; j < Jmax; j++) {

            psi_Exz_2[i][j] = (double *)malloc(nzPML_2 * sizeof(double));

            for(k = 0; k < nzPML_2; k++) {

                psi_Exz_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Eyz_1 = (double ***)malloc((Imax-1) * sizeof(double **));

    for(i = 0; i < Imax; i++) {

        psi_Eyz_1[i] = (double **)malloc((Jmax-1) * sizeof(double *));

        for(j = 0; j < Jmax-1; j++) {

            psi_Eyz_1[i][j] = (double *)malloc(nzPML_1 * sizeof(double));

            for(k = 0; k < nzPML_1; k++) {

                psi_Eyz_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Eyz_2 = (double ***)malloc((Imax-1) * sizeof(double **));

    for(i = 0; i < Imax; i++) {

        psi_Eyz_2[i] = (double **)malloc((Jmax-1) * sizeof(double *));

        for(j = 0; j < Jmax-1; j++) {

            psi_Eyz_2[i][j] = (double *)malloc(nzPML_2 * sizeof(double));

            for(k = 0; k < nzPML_2; k++) {

                psi_Eyz_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Hzx_1 = (double ***)malloc((nxPML_1-1) * sizeof(double **));

    for(i = 0; i < nxPML_1-1; i++) {

        psi_Hzx_1[i] = (double **)malloc((Jmax-1) * sizeof(double *));

        for(j = 0; j < Jmax-1; j++) {

            psi_Hzx_1[i][j] = (double *)malloc((Kmax-1) * sizeof(double));

            for(k = 0; k < Kmax-1; k++) {

                psi_Hzx_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Hzx_2 = (double ***)malloc((nxPML_2-1) * sizeof(double **));

    for(i = 0; i < nxPML_2-1; i++) {

        psi_Hzx_2[i] = (double **)malloc((Jmax-1) * sizeof(double *));

        for(j = 0; j < Jmax-1; j++) {

            psi_Hzx_2[i][j] = (double *)malloc((Kmax-1) * sizeof(double));

            for(k = 0; k < Kmax-1; k++) {

                psi_Hzx_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Eyx_1 = (double ***)malloc(nxPML_1 * sizeof(double **));

    for(i = 0; i < nxPML_1; i++) {

        psi_Eyx_1[i] = (double **)malloc((Jmax-1) * sizeof(double *));

        for(j = 0; j < Jmax-1; j++) {

            psi_Eyx_1[i][j] = (double *)malloc((Kmax-1) * sizeof(double));

            for(k = 0; k < Kmax-1; k++) {

                psi_Eyx_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Eyx_2 = (double ***)malloc(nxPML_2 * sizeof(double **));

    for(i = 0; i < nxPML_2; i++) {

        psi_Eyx_2[i] = (double **)malloc((Jmax-1) * sizeof(double *));

        for(j = 0; j < Jmax-1; j++) {

            psi_Eyx_2[i][j] = (double *)malloc((Kmax-1) * sizeof(double));

            for(k = 0; k < Kmax-1; k++) {

                psi_Eyx_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Hzy_1 = (double ***)malloc((Imax-1) * sizeof(double **));

    for(i = 0; i < Imax-1; i++) {

        psi_Hzy_1[i] = (double **)malloc((nyPML_1-1) * sizeof(double *));

        for(j = 0; j < nyPML_1-1; j++) {

            psi_Hzy_1[i][j] = (double *)malloc((Kmax-1) * sizeof(double));

            for(k = 0; k < Kmax-1; k++) {

                psi_Hzy_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Hzy_2 = (double ***)malloc((Imax-1) * sizeof(double **));

    for(i = 0; i < Imax-1; i++) {

        psi_Hzy_2[i] = (double **)malloc((nyPML_2-1) * sizeof(double *));

        for(j = 0; j < nyPML_2-1; j++) {

            psi_Hzy_2[i][j] = (double *)malloc((Kmax-1) * sizeof(double));

            for(k = 0; k < Kmax-1; k++) {

                psi_Hzy_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Exy_1 = (double ***)malloc((Imax-1) * sizeof(double **));

    for(i = 0; i < Imax-1; i++) {

        psi_Exy_1[i] = (double **)malloc(nyPML_1 * sizeof(double *));

        for(j = 0; j < nyPML_1; j++) {

            psi_Exy_1[i][j] = (double *)malloc((Kmax-1) * sizeof(double));

            for(k = 0; k < Kmax-1; k++) {

                psi_Exy_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Exy_2 = (double ***)malloc((Imax-1) * sizeof(double **));

    for(i = 0; i < Imax-1; i++) {

        psi_Exy_2[i] = (double **)malloc(nyPML_2 * sizeof(double *));

        for(j = 0; j < nyPML_2; j++) {

            psi_Exy_2[i][j] = (double *)malloc((Kmax-1) * sizeof(double));

            for(k = 0; k < Kmax-1; k++) {

                psi_Exy_2[i][j][k] = 0.0;
            }
        }
    }

    be_x_1 = (double *)malloc((nxPML_1) * sizeof(double));
    for(i = 0; i < nxPML_1; i++) {

        be_x_1[i] = 0.0;
    }

    ce_x_1 = (double *)malloc((nxPML_1) * sizeof(double));
    for(i = 0; i < nxPML_1; i++) {

        ce_x_1[i] = 0.0;
    }

    alphae_x_PML_1 = (double *)malloc((nxPML_1) * sizeof(double));
    for(i = 0; i < nxPML_1; i++) {

        alphae_x_PML_1[i] = 0.0;
    }

    sige_x_PML_1 = (double *)malloc((nxPML_1) * sizeof(double));
    for(i = 0; i < nxPML_1; i++) {

        sige_x_PML_1[i] = 0.0;
    }

    kappae_x_PML_1 = (double *)malloc((nxPML_1) * sizeof(double));
    for(i = 0; i < nxPML_1; i++) {

        kappae_x_PML_1[i] = 0.0;
    }

    bh_x_1 = (double *)malloc((nxPML_1-1) * sizeof(double));
    for(i = 0; i < nxPML_1-1; i++) {

        bh_x_1[i] = 0.0;
    }

    ch_x_1 = (double *)malloc((nxPML_1-1) * sizeof(double));
    for(i = 0; i < nxPML_1-1; i++) {

        ch_x_1[i] = 0.0;
    }

    alphah_x_PML_1 = (double *)malloc((nxPML_1-1) * sizeof(double));
    for(i = 0; i < nxPML_1-1; i++) {

        alphah_x_PML_1[i] = 0.0;
    }

    sigh_x_PML_1 = (double *)malloc((nxPML_1-1) * sizeof(double));
    for(i = 0; i < nxPML_1-1; i++) {

        sigh_x_PML_1[i] = 0.0;
    }

    kappah_x_PML_1 = (double *)malloc((nxPML_1-1) * sizeof(double));
    for(i = 0; i < nxPML_1-1; i++) {

        kappah_x_PML_1[i] = 0.0;
    }

    be_x_2 = (double *)malloc((nxPML_2) * sizeof(double));
    for(i = 0; i < nxPML_2; i++) {

        be_x_2[i] = 0.0;
    }

    ce_x_2 = (double *)malloc((nxPML_2) * sizeof(double));
    for(i = 0; i < nxPML_2; i++) {

        ce_x_2[i] = 0.0;
    }

    alphae_x_PML_2 = (double *)malloc((nxPML_2) * sizeof(double));
    for(i = 0; i < nxPML_2; i++) {

        alphae_x_PML_2[i] = 0.0;
    }


    sige_x_PML_2 = (double *)malloc((nxPML_2) * sizeof(double));
    for(i = 0; i < nxPML_2; i++) {

        sige_x_PML_2[i] = 0.0;
    }


    kappae_x_PML_2 = (double *)malloc((nxPML_2) * sizeof(double));
    for(i = 0; i < nxPML_2; i++) {

        kappae_x_PML_2[i] = 0.0;
    }


    bh_x_2 = (double *)malloc((nxPML_2-1) * sizeof(double));
    for(i = 0; i < nxPML_2-1; i++) {

        bh_x_2[i] = 0.0;
    }


    ch_x_2 = (double *)malloc((nxPML_2-1) * sizeof(double));
    for(i = 0; i < nxPML_2-1; i++) {

        ch_x_2[i] = 0.0;
    }

    alphah_x_PML_2 = (double *)malloc((nxPML_2-1) * sizeof(double));
    for(i = 0; i < nxPML_2-1; i++) {

        alphah_x_PML_2[i] = 0.0;
    }

    sigh_x_PML_2 = (double *)malloc((nxPML_2-1) * sizeof(double));
    for(i = 0; i < nxPML_2-1; i++) {

        sigh_x_PML_2[i] = 0.0;
    }

    kappah_x_PML_2 = (double *)malloc((nxPML_2-1) * sizeof(double));
    for(i = 0; i < nxPML_1-1; i++) {

        kappah_x_PML_2[i] = 0.0;
    }

    be_y_1 = (double *)malloc((nyPML_1) * sizeof(double));
    for(i = 0; i < nyPML_1; i++) {

        be_y_1[i] = 0.0;
    }

    ce_y_1 = (double *)malloc((nyPML_1) * sizeof(double));
    for(i = 0; i < nyPML_1; i++) {

        ce_y_1[i] = 0.0;
    }

    alphae_y_PML_1 = (double *)malloc((nyPML_1) * sizeof(double));
    for(i = 0; i < nyPML_1; i++) {

        alphae_y_PML_1[i] = 0.0;
    }

    sige_y_PML_1 = (double *)malloc((nyPML_1) * sizeof(double));
    for(i = 0; i < nyPML_1; i++) {

        sige_y_PML_1[i] = 0.0;
    }

    kappae_y_PML_1 = (double *)malloc((nyPML_1) * sizeof(double));
    for(i = 0; i < nyPML_1; i++) {

        kappae_y_PML_1[i] = 0.0;
    }

    bh_y_1 = (double *)malloc((nyPML_1-1) * sizeof(double));
    for(i = 0; i < nyPML_1-1; i++) {

        bh_y_1[i] = 0.0;
    }

    ch_y_1 = (double *)malloc((nyPML_1-1) * sizeof(double));
    for(i = 0; i < nyPML_1-1; i++) {

        ch_y_1[i] = 0.0;
    }

    alphah_y_PML_1 = (double *)malloc((nyPML_1-1) * sizeof(double));
    for(i = 0; i < nyPML_1-1; i++) {

        alphah_y_PML_1[i] = 0.0;
    }

    sigh_y_PML_1 = (double *)malloc((nyPML_1-1) * sizeof(double));
    for(i = 0; i < nyPML_1-1; i++) {

        sigh_y_PML_1[i] = 0.0;
    }

    kappah_y_PML_1 = (double *)malloc((nyPML_1-1) * sizeof(double));
    for(i = 0; i < nyPML_1-1; i++) {

        kappah_y_PML_1[i] = 0.0;
    }

    be_y_2 = (double *)malloc((nyPML_2) * sizeof(double));
    for(i = 0; i < nyPML_2; i++) {

        be_y_2[i] = 0.0;
    }

    ce_y_2 = (double *)malloc((nyPML_2) * sizeof(double));
    for(i = 0; i < nyPML_2; i++) {

        ce_y_2[i] = 0.0;
    }

    alphae_y_PML_2 = (double *)malloc((nyPML_2) * sizeof(double));
    for(i = 0; i < nyPML_2; i++) {

        alphae_y_PML_2[i] = 0.0;
    }

    sige_y_PML_2 = (double *)malloc((nyPML_2) * sizeof(double));
    for(i = 0; i < nyPML_2; i++) {

        sige_y_PML_2[i] = 0.0;
    }

    kappae_y_PML_2 = (double *)malloc((nyPML_2) * sizeof(double));
    for(i = 0; i < nyPML_2; i++) {

        kappae_y_PML_2[i] = 0.0;
    }

    bh_y_2 = (double *)malloc((nyPML_2-1) * sizeof(double));
    for(i = 0; i < nyPML_2-1; i++) {

        bh_y_2[i] = 0.0;
    }

    ch_y_2 = (double *)malloc((nyPML_2-1) * sizeof(double));
    for(i = 0; i < nyPML_2-1; i++) {

        ch_y_2[i] = 0.0;
    }

    alphah_y_PML_2 = (double *)malloc((nyPML_2-1) * sizeof(double));
    for(i = 0; i < nyPML_2-1; i++) {

        alphah_y_PML_2[i] = 0.0;
    }

    sigh_y_PML_2 = (double *)malloc((nyPML_2-1) * sizeof(double));
    for(i = 0; i < nyPML_2-1; i++) {

        sigh_y_PML_2[i] = 0.0;
    }

    kappah_y_PML_2 = (double *)malloc((nyPML_2-1) * sizeof(double));
    for(i = 0; i < nyPML_1-1; i++) {

        kappah_y_PML_2[i] = 0.0;
    }

    be_z_1 = (double *)malloc((nzPML_1) * sizeof(double));
    for(i = 0; i < nzPML_1; i++) {

        be_z_1[i] = 0.0;
    }

    ce_z_1 = (double *)malloc((nzPML_1) * sizeof(double));
    for(i = 0; i < nzPML_1; i++) {

        ce_z_1[i] = 0.0;
    }

    alphae_z_PML_1 = (double *)malloc((nzPML_1) * sizeof(double));
    for(i = 0; i < nzPML_1; i++) {

        alphae_z_PML_1[i] = 0.0;
    }

    sige_z_PML_1 = (double *)malloc((nzPML_1) * sizeof(double));
    for(i = 0; i < nzPML_1; i++) {

        sige_z_PML_1[i] = 0.0;
    }

    kappae_z_PML_1 = (double *)malloc((nzPML_1) * sizeof(double));
    for(i = 0; i < nzPML_1; i++) {

        kappae_z_PML_1[i] = 0.0;
    }

    bh_z_1 = (double *)malloc((nzPML_1-1) * sizeof(double));
    for(i = 0; i < nzPML_1-1; i++) {

        bh_z_1[i] = 0.0;
    }

    ch_z_1 = (double *)malloc((nzPML_1-1) * sizeof(double));
    for(i = 0; i < nzPML_1-1; i++) {

        ch_z_1[i] = 0.0;
    }

    alphah_z_PML_1 = (double *)malloc((nzPML_1-1) * sizeof(double));
    for(i = 0; i < nzPML_1-1; i++) {

        alphah_z_PML_1[i] = 0.0;
    }

    sigh_z_PML_1 = (double *)malloc((nzPML_1-1) * sizeof(double));
    for(i = 0; i < nzPML_1-1; i++) {

        sigh_z_PML_1[i] = 0.0;
    }

    kappah_z_PML_1 = (double *)malloc((nzPML_1-1) * sizeof(double));
    for(i = 0; i < nzPML_1-1; i++) {

        kappah_z_PML_1[i] = 0.0;
    }

    be_z_2 = (double *)malloc((nzPML_2) * sizeof(double));
    for(i = 0; i < nzPML_2; i++) {

        be_z_2[i] = 0.0;
    }

    ce_z_2 = (double *)malloc((nzPML_2) * sizeof(double));
    for(i = 0; i < nzPML_2; i++) {

        ce_z_2[i] = 0.0;
    }

    alphae_z_PML_2 = (double *)malloc((nzPML_2) * sizeof(double));
    for(i = 0; i < nzPML_2; i++) {

        alphae_z_PML_2[i] = 0.0;
    }

    sige_z_PML_2 = (double *)malloc((nzPML_2) * sizeof(double));
    for(i = 0; i < nzPML_2; i++) {

        sige_z_PML_2[i] = 0.0;
    }

    kappae_z_PML_2 = (double *)malloc((nzPML_2) * sizeof(double));
    for(i = 0; i < nzPML_2; i++) {

        kappae_z_PML_2[i] = 0.0;
    }

    bh_z_2 = (double *)malloc((nzPML_2-1) * sizeof(double));
    for(i = 0; i < nzPML_2-1; i++) {

        bh_z_2[i] = 0.0;
    }

    ch_z_2 = (double *)malloc((nzPML_2-1) * sizeof(double));
    for(i = 0; i < nzPML_2-1; i++) {

        ch_z_2[i] = 0.0;
    }

    alphah_z_PML_2 = (double *)malloc((nzPML_2-1) * sizeof(double));
    for(i = 0; i < nzPML_2-1; i++) {

        alphah_z_PML_2[i] = 0.0;
    }

    sigh_z_PML_2 = (double *)malloc((nzPML_2-1) * sizeof(double));
    for(i = 0; i < nzPML_2-1; i++) {

        sigh_z_PML_2[i] = 0.0;
    }


    kappah_z_PML_2 = (double *)malloc((nzPML_2-1) * sizeof(double));
    for(i = 0; i < nzPML_1-1; i++) {

        kappah_z_PML_2[i] = 0.0;
    }


}




