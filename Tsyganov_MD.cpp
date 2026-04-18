#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "params.h"
#include "constans.h"
#include "start_cond.h"
#include "global_var.h"

using namespace std;

// Добавим глобальные переменные для записи
double rx_ij, ry_ij, rz_ij, r_ij, r_ij2;
double ex_ij, ey_ij, ez_ij;

// Выделение памяти
void memoryAllocation() {
    coordx = (double*) malloc(NUMBERPARTICLES * sizeof(double));
    coordy = (double*) malloc(NUMBERPARTICLES * sizeof(double));
    coordz = (double*) malloc(NUMBERPARTICLES * sizeof(double));

    vx = (double*) malloc(NUMBERPARTICLES * sizeof(double));
    vy = (double*) malloc(NUMBERPARTICLES * sizeof(double));
    vz = (double*) malloc(NUMBERPARTICLES * sizeof(double));

    Fx = (double*) malloc(NUMBERPARTICLES * sizeof(double));
    Fy = (double*) malloc(NUMBERPARTICLES * sizeof(double));
    Fz = (double*) malloc(NUMBERPARTICLES * sizeof(double));
}

// очищаем память
void memoryFree() {

    free(coordx); coordx = NULL;
    free(coordy); coordy = NULL;
    free(coordz); coordz = NULL;

    free(vx); vx = NULL;
    free(vy); vy = NULL;
    free(vz); vz = NULL;

    free(Fx); Fx = NULL;
    free(Fy); Fy = NULL;
    free(Fz); Fz = NULL;
}


void calculate_potential_U_F() {

    U = 0.0;
    F = 0.0;
    for (int k = 0; k < NUMBERPARTICLES; ++k) {
        Fx[k] = 0.0;
        Fy[k] = 0.0;
        Fz[k] = 0.0;
    }

    for (int i = 0; i < NUMBERPARTICLES; ++i) {
        for (int j = i + 1; j < NUMBERPARTICLES; ++j) {

            rx_ij = coordx[i] - coordx[j];
            ry_ij = coordy[i] - coordy[j];
            rz_ij = coordz[i] - coordz[j];

            // виртуальные частицы
            if (rx_ij >  LX / 2.0) rx_ij -= LX;
            if (rx_ij < -LX / 2.0) rx_ij += LX;

            if (ry_ij >  LY / 2.0) ry_ij -= LY;
            if (ry_ij < -LY / 2.0) ry_ij += LY;

            if (rz_ij >  LZ / 2.0) rz_ij -= LZ;
            if (rz_ij < -LZ / 2.0) rz_ij += LZ;
            
            // находим квадрат расстояния
            r_ij2 = rx_ij * rx_ij + ry_ij * ry_ij + rz_ij * rz_ij;

            
            if (r_ij2 <= RCUT2) { 

                
                // растояние между частицами
                r_ij = sqrt(r_ij2);

                // вводим доп переменные
                double sr = SIGMA / r_ij;
                double sr2 = sr * sr;
                double sr6 = sr2 * sr2 * sr2;

                // вспомогательные единичные векторы
                ex_ij = rx_ij / r_ij;
                ey_ij = ry_ij / r_ij;
                ez_ij = rz_ij / r_ij;

                U += 4.0 * EPS * (sr6 * sr6 - sr6) - UCUT;
                F = 24.0 * EPS / r_ij * (2.0 * sr6 * sr6 - sr6);

                // силы (3 закон Ньютона)
                Fx[i] += F * ex_ij;
                Fy[i] += F * ey_ij;
                Fz[i] += F * ez_ij;

                Fx[j] -= F * ex_ij;
                Fy[j] -= F * ey_ij;
                Fz[j] -= F * ez_ij;
            } 
        }
    }
}

// Velocity Verlet (половинный шаг скоростей)
void velocity_Verlet_half() {
    for (int i = 0; i < NUMBERPARTICLES; ++i) {
        // v(t + dt/2) = v(t) + F(t) * dt / (2m)
        vx[i] += 0.5 * STEP * Fx[i] / MASS;
        vy[i] += 0.5 * STEP * Fy[i] / MASS;
        vz[i] += 0.5 * STEP * Fz[i] / MASS;
    }
}

// ПГУ, периодические граничные условия
void PBC() {
    for (int i = 0; i < NUMBERPARTICLES; ++i) {

        if (coordx[i] >= LX) coordx[i] -= LX;
        if (coordx[i] < 0.0) coordx[i] += LX;

        if (coordy[i] >= LY) coordy[i] -= LY;
        if (coordy[i] < 0.0) coordy[i] += LY;

        if (coordz[i] >= LZ) coordz[i] -= LZ;
        if (coordz[i] < 0.0) coordz[i] += LZ;
    }
}


// вычисление координат
void coord_Verlet() {
    for (int i = 0; i < NUMBERPARTICLES; ++i) {
        // r(t + dt) = r(t) + v(t + dt/2) * dt
        coordx[i] += STEP * vx[i];
        coordy[i] += STEP * vy[i];
        coordz[i] += STEP * vz[i];
    }
}

// запись в файл
void writeToFile(int step) {

    FILE *file = fopen("Tsyganov_MD_8.txt", "a");

    

    fprintf(file, "Step = %d\n", step);

    fprintf(file, "r1 = (rx1; ry1; rz1) = (%0.8f; %0.8f; %0.8f)\n", coordx[0], coordy[0], coordz[0]);
    fprintf(file, "r2 = (rx2; ry2; rz2) = (%0.8f; %0.8f; %0.8f)\n", coordx[1], coordy[1], coordz[1]);

    fprintf(file, "r12 = (rx12; ry12; rz12) = (%0.8f; %0.8f; %0.8f)\n", rx_ij, ry_ij, rz_ij);
    fprintf(file, "r12_abs = %0.8f\n", r_ij);

    fprintf(file, "(rx12; ry12; rz12) / r12_abs = (%0.8f; %0.8f; %0.8f)\n", ex_ij, ey_ij, ez_ij);

    fprintf(file, "U12 = %0.8f\n", U);
    fprintf(file, "F12 = %0.8f\n", F);
    fprintf(file, "F1 = (Fx1; Fy1; Fz1) = (%0.8f; %0.8f; %0.8f)\n", Fx[0], Fy[0], Fz[0]);

    fprintf(file, "v1 = (vx1; vy1; vz1) = (%0.8f; %0.8f; %0.8f)\n", vx[0], vy[0], vz[0]);
    fprintf(file, "v2 = (vx2; vy2; vz2) = (%0.8f; %0.8f; %0.8f)\n", vx[1], vy[1], vz[1]);

    fprintf(file, "\n");
    
    fclose(file);
}

// основной алгоритм 
void algorithm_MD_problem() {

    // задаем начальные условия для частиц
    start_cond_two_particles();

    // считаем нулевой шаг
    calculate_potential_U_F();

    writeToFile(0);

    
    for (int step = 1; step <= LASTSTEP; step++){

        velocity_Verlet_half();

        coord_Verlet();

        // ПГУ, смещение координат при необходимости
        PBC();
        
        calculate_potential_U_F();

        velocity_Verlet_half();

        writeToFile(step);
    }
}

int main() {
    
    // веделяем память
    memoryAllocation();

    // запускаем алгоритм
    algorithm_MD_problem();

    // очищаем память
    memoryFree();

    return 0;
}