#include <stdio.h>

// координаты (массивы)
double *coordx =NULL;
double *coordy =NULL;
double *coordz =NULL;

// скорости (массивы)
double *vx = NULL;
double *vy = NULL;
double *vz = NULL;

// потенциальная, кин., кин. тепл., полная, внутр. энергии (скаляры)
double U = 0.0;
double Ekin = 0.0;
double Eterm = 0.0;
double E = 0.0;
double Eint = 0.0;

// силы взаимодействия 
double *Fx = NULL;
double *Fy = NULL;
double *Fz = NULL;
double F = 0.0;


