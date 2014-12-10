#include "fdtd_complx.hpp"
#include <list>
using namespace std;
const double    c = 2.99795637e+8;
const double	pi = 3.14159265358979323846;
const double	eps0 = 8.854e-12;
const double	mu0 = 1.2566306e-6;
const double	eta = 376.733340625715;

typedef struct wire
{
    int is,ie;
    int js,je;
    int ks;

    double hy_buf[2][100];
    double hx_buf[2][100];
    int length;
    double radius;
} WIRE;


typedef struct resistor
{
    int i;
    int j;
    int k;

    double value;
    double* v_in;
    double* i_in;
    double  Ez_Old;
    complx* V_in;
    complx* I_in;;

} RESISTOR;

typedef struct capacitor
{
    int i;
    int j;
    int k;

    double value;
    double* v_in;
    double* i_in;
    double  Ez_Old;
    complx* V_in;
    complx* I_in;;

} CAPACITOR;

typedef struct inductor
{
    int i;
    int j;
    int k;

    double value;
    double* v_in;
    double* i_in;
    double  Ez_Old;
    complx* V_in;
    complx* I_in;
    double Jx,Jy,Jz;

} INDUCTOR;


class PORT
{
public:
    int i;
    int j;
    int k;

    double zo;
    double* v_in;
    double* i_in;
    double Ezn1;

    /* Parameters for Excitation */

    // Gaussian Pulse
    // parameters[0] - amplitude ; parameters[1] - t0; parameters[2] - pulse width

    // Sinusoidal Pulse
    // parameters[0] - amplitude ; parameters[1] - frequency; parameters[2] - not used

    double parameters[3];

    complx* V_in;
    complx* I_in;

    double (*signal)(double,double,double,double);

    int polarisation; /* 0 - Ex; 1 - Ey; 2 - Ez */
    int type; /* 0 - Active ; 1 - Passive */
};



class vPROBE
{
public:
    int i;
    int j;
    int k;

    double* v_in;
    complx* V_in;

    int orientation; /* 0 - Ex; 1 - Ey; 2 - Ez */
};




typedef struct grid
{
    double *Ex;
    double *Ey;
    double *Ez;

//	double *Dx;
//	double *Dy;
//	double *Dz;

    double *Hx;
    double *Hy;
    double *Hz;

    int *xe;
    int *ye;
    int *ze;

    int *xm;
    int *ym;
    int *zm;

    int exnx;
    int exny;
    int exnz;

    int eynx;
    int eyny;
    int eynz;

    int eznx;
    int ezny;
    int eznz;

    int hxnx;
    int hxny;
    int hxnz;

    int hynx;
    int hyny;
    int hynz;

    int hznx;
    int hzny;
    int hznz;
} GRID;

typedef struct xsurfaces
{
    complx* Ey;
    complx* Ez;

    complx* Hy;
    complx* Hz;

    complx* My;
    complx* Mz;

    complx* Jy;
    complx* Jz;


} XSURFACES;

typedef struct ysurfaces
{
    complx* Ex;
    complx* Ez;
    complx* Hx;
    complx* Hz;

    complx* Mx;
    complx* Mz;

    complx* Jx;
    complx* Jz;

} YSURFACES;

typedef struct zsurfaces
{
    complx* Ex;
    complx* Ey;
    complx* Hx;
    complx* Hy;

    complx* Mx;
    complx* My;

    complx* Jx;
    complx* Jy;

} ZSURFACES;


class CLOSED_SURFACE
{
    int is,js,ks;
    int ie,je,ke;
    int iw,jw,kw;

    int xc,yc,zc; /* Reference point for Far-Field XTransform */
    XSURFACES x_surfaces[2];
    YSURFACES y_surfaces[2];
    ZSURFACES z_surfaces[2];


public:
    double frequency;
    double prad;
    double pin;
    double efficiency;
    int label;

    /* Methods */
    void allocate_surface_memory(int,int,int,int);
    void perform_DFT(GRID *grid,int n, double dt);
    void compute_surface_currents();
    void compute_farfield_xzplane(double dx,double dy, double dz,int);
    void compute_farfield_xyplane(double dx,double dy, double dz,int);
    void compute_radiated_power(list<PORT>::iterator,double,double,double);
};


class SIMULATION
{
public:
    int nx,ny,nz,n,Nt;
    double dx,dy,dz;
    double dt;
    int number_of_freqs;
    list<double> frequencies;

    list<PORT> list_of_ports;
    list<CLOSED_SURFACE> list_of_closed_surfaces;
    list<RESISTOR> list_of_resistors;
    list<INDUCTOR> list_of_inductors;
    list<CAPACITOR> list_of_capacitors;
    list<WIRE> list_of_wires;

    /* Methods */
    int getFrequencyLength();
    void set_simulation_length(int N);
    void set_mesh_size(int Nx,int Ny, int Nz);
    void set_mesh_resolution(double dx, double dy, double dz);
    void initClosedSurface();
    void updateClosedSurfaces(GRID* grid);
    void addResistor(int i,int j,int k,double value);
    void addCapacitor(int i,int j,int k,double value);
    void addInductor(int i,int j,int k,double value);
    void addWire(int,int,int,int,double);

    /* Ports */
    void addFrequency(double);
    void addPort(int i, int j, int k, double zo,double,double,double, string);



};

typedef struct constants
{
    double* ca;
    double* cbx;
    double* cby;
    double* cbz;


    double* da;
    double* dbx;
    double* dby;
    double* dbz;

    double* Er;
    double* Ur;
} FD_CONSTANTS;







