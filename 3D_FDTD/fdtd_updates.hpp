void Ex_Update(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants);
void Ey_Update(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants);
void Ez_Update(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants);

void Hx_Update(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants);
void Hy_Update(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants);
void Hz_Update(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants);

void updateExcitation(GRID* grid, SIMULATION* simulation);
void updateResistors(GRID*, SIMULATION*, FD_CONSTANTS*);
void updateCapacitors(GRID*, SIMULATION*, FD_CONSTANTS*);
void updateInductors(GRID*, SIMULATION*, FD_CONSTANTS*);
void updateWires(GRID* grid, SIMULATION* simulation,FD_CONSTANTS* constants);

