#include <stdio.h>

#include "global_var.h"


//Создать функцию для задания н.у. двух частиц:
void start_cond_two_particles() {
    coordx[0] = 0.75; 
    coordy[0] = 0.75; 
    coordz[0] = 0.5;
    
    coordx[1] = 1.25; 
    coordy[1] = 0.75; 
    coordz[1] = 0.5;
    
    vx[0] = 1.0; 
    vy[0] = 1.0; 
    vz[0] = 0.0;

    vx[1] =-1.0;
    vy[1] = 1.0; 
    vz[1] = 0.0;
}

void start_cond_four_particles() {
    coordx[0] = 0.25; 
    coordy[0] = 0.75; 
    coordz[0] = 0.5;
    
    coordx[1] = 0.75; 
    coordy[1] = 0.75; 
    coordz[1] = 0.5;

    coordx[2] = 1.35; 
    coordy[2] = 1.75; 
    coordz[2] = 0.5;

    coordx[3] = 1.68; 
    coordy[3] = 1.75; 
    coordz[3] = 0.5;
    
    vx[0] = 1.0; 
    vy[0] = 1.0; 
    vz[0] = 0.0;

    vx[1] = -1.0;
    vy[1] = -1.0; 
    vz[1] = 0.0;

    vx[2] = 1.0;
    vy[2] = 2.0; 
    vz[2] = 0.0;

    vx[3] = -1.0;
    vy[3] = -2.0; 
    vz[3] = 0.0;
}