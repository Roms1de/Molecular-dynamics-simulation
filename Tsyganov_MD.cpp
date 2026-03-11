#include <iostream>
#include <fstream>
#include <iomanip>

#include "params.h"
#include "constans.h"
#include "start_cond.h"
#include "global_var.h"


using namespace std;

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

// очищаем память и задаем нул, чтобы не было крашей
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


// запись в файл
void writeToFile() {


    int currentStep = 0;

    double rx1, ry1, rz1, rx2, ry2, rz2;
    double rx12, ry12, rz12, rx21, ry21, rz21, r12_abs, r21_abs;
    double ex12, ey12, ez12, ex21, ey21, ez21;

    double vx1, vy1, vz1, vx2, vy2, vz2;

    double Fx12, Fy12, Fz12, Fx21, Fy21, Fz21, Fs;
    double U12, sr_abs, sr_abs2, sr_abs6;
    

    // Координаты первой частицы 
    rx1 = coordx[0];
    ry1 = coordy[0];
    rz1 = coordz[0];

    // Координаты второй частицы
    rx2 = coordx[1];
    ry2 = coordy[1];
    rz2 = coordz[1];

    // Скорости первой частицы
    vx1 = vx[0];
    vy1 = vy[0];
    vz1 = vz[0];

    // Скорости второй частицы
    vx2 = vx[1];
    vy2 = vy[1];
    vz2 = vz[1];

    // Радиус векторы
    rx12 = rx1 - rx2;
    ry12 = ry1 - ry2;
    rz12 = rz1 - rz2;

    rx21 = rx2 - rx1;
    ry21 = ry2 - ry1;
    rz21 = rz2 - rz1;

    // Длина вектора r12
    r12_abs = sqrt(rx12 * rx12 + ry12 * ry12 + rz12 * rz12);

    // Длина вектора r21
    r21_abs = sqrt(rx21 * rx21 + ry21 * ry21 + rz21 * rz21);

    //вычисляем единичный вектор 
    ex12 = rx12 / r12_abs;
    ey12 = ry12 / r12_abs;
    ez12 = rz12 / r12_abs;

    ex21 = rx21 / r21_abs;
    ey21 = ry21 / r21_abs;
    ez21 = rz21 / r21_abs;


    // Вспомогательные переменные
    sr_abs = SIGMA / r12_abs;
    sr_abs2 = sr_abs * sr_abs;
    sr_abs6 = sr_abs2 * sr_abs2 * sr_abs2;

    // Потенциал взаимодействия
    U12 = 4.0 * EPS * (sr_abs6 * sr_abs6 - sr_abs6);

    // Сила взаимодействия 
    Fs = 24.0 * EPS / r12_abs * (2.0 * sr_abs6 * sr_abs6 - sr_abs6);

    //проекции силы, действующей на первую частицу
    Fx12 = Fs * ex12;
    Fy12 = Fs * ey12;
    Fz12 = Fs * ez12;

    //сила действующая на вторую частицу
    Fx21 = -Fx12;
    Fy21 = -Fy12;
    Fz21 = -Fz12;



    // Создаем новый файл или перезаписываем
    FILE *file = fopen("Tsyganov_MD_3.txt", "w");

    fprintf(file, "Step = %d\n", currentStep);
    fprintf(file, "r1 = (rx1; ry1; rz1) = (%0.8f; %0.8f; %0.8f)\n", rx1, ry1, rz1);
    fprintf(file, "r2 = (rx2; ry2; rz2) = (%0.8f; %0.8f; %0.8f)\n", rx2, ry2, rz2);
    fprintf(file, "r12 = (rx12; ry12; rz12) = (%0.8f; %0.8f; %0.8f)\n", rx12, ry12, rz12);
    fprintf(file, "r12_abs = %0.8f\n", r12_abs);
    fprintf(file, "(rx12; ry12; rz12) / r12_abs = (%0.8f; %0.8f; %0.8f)\n", ex12, ey12, ez12);
    fprintf(file, "(rx21; ry21; rz21) / r21_abs = (%0.8f; %0.8f; %0.8f)\n", ex21, ey21, ez21);
    fprintf(file, "v1 = (vx1; vy1, vz1) = (%0.8f; %0.8f; %0.8f)\n", vx1, vy1, vz1);
    fprintf(file, "v2 = (vx2; vy2, vz2) = (%0.8f; %0.8f; %0.8f)\n", vx2, vy2, vz2);
    fprintf(file, "U12 = %0.8f\n", U12);
    fprintf(file, "F = %0.8f\n", Fs);
    fprintf(file, "F1 = (Fx1; Fy1; Fz1) = (%0.8f; %0.8f; %0.8f)\n", Fx12, Fy12, Fz12);
    fprintf(file, "F2 = (Fx2; Fy2; Fz2) = (%0.8f; %0.8f; %0.8f)\n", Fx21, Fy21, Fz21);

    
    fclose(file);
}

// Основная функция-алгоритм 
void algorithm_MD_problem() {
    start_cond_two_particles();
    writeToFile();
}

int main() {
    
    memoryAllocation();

    algorithm_MD_problem();

    memoryFree();

    return 0;
}

