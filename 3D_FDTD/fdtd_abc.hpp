

typedef struct previous_Efields
{
    double*	Ex_old;
    double*	Ey_old;
    double*	Ez_old;

    double*	Ex_old_old;
    double*	Ey_old_old;
    double*	Ez_old_old;

    double*	Ex_old_old_old;
    double*	Ey_old_old_old;
    double*	Ez_old_old_old;

} PREVIOUS_EFIELDS;

void record_previous_fields(GRID*, SIMULATION*,PREVIOUS_EFIELDS*);
void abc_liao(GRID*, SIMULATION*,FD_CONSTANTS* ,PREVIOUS_EFIELDS*);
void abc_murfirst(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants,PREVIOUS_EFIELDS* previous_Efields);


