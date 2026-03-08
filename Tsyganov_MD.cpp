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


    ofstream file("Tsyganov_MD_2.txt");

    // огругление вывода вещественных чисел до 8 знаков после запятой
    file << fixed << setprecision(8);

    file << endl;
    file << "Number of Particles = " << NUMBERPARTICLES << endl;
    file << "Number of Steps = " << NSTEPS << endl;
    file << "Size of System on X = " << LX << endl;
    file << "Boltzmann Constant = " << K_B << endl;
    file << "Volume = " << VOLUME << endl;

    file.close();
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

