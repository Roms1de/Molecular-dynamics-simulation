#pragma once
#include <stdio.h>



// координаты (массивы)
double *coordx =NULL;
double *coordy =NULL;
double *coordz =NULL;

// скорости (массивы)
double *vx = NULL;
double *vy = NULL;
double *vz = NULL;
static double v_mx = 0.0;
static double v_my = 0.0;
static double v_mz = 0.0;

// потенциальная, кин., кин. тепл., полная, внутр. энергии (скаляры)
double U = 0.0;
double U_full;
double Ekin = 0.0;
double Eterm = 0.0;
double E = 0.0;
double Eint = 0.0;

// силы взаимодействия 
double *Fx = NULL;
double *Fy = NULL;
double *Fz = NULL;
double F = 0.0;

// Температура
double T = 0.0;

// Давление и его компоненты
double P = 0.0;
double P1 = 0.0; // Кинетическая часть
double P2 = 0.0; // Вириальная (силовая) часть

