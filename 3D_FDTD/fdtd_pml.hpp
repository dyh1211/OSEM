typedef struct pml_parameters
{
    double* Psi_eyx_xn;
    double* Psi_ezx_xn;
    double* Psi_hyx_xn;
    double* Psi_hzx_xn;

    double* CPsi_eyx_xn;
    double* CPsi_ezx_xn;
    double* CPsi_hyx_xn;
    double* CPsi_hzx_xn;

    double* cpml_b_ex_xn;
    double* cpml_a_ex_xn;
    double* cpml_b_mx_xn;
    double* cpml_a_mx_xn;
} PML_PARAMETERS;

