#include <stdio.h>


#include "params.h"


// масса одной частицы (атом аргона)
const double MASS = 66.335;

//постоянная Больцмана
const double K_B = 1.380648528; 

//параметры потенциала Леннарда-Джонса для аргона
const double EPS = 1.712; //параметр энергии (глубина потенциальной ямы)
const double SIGMA = 0.3418; //параметр длины взаимодействия
const double RCUT = 2.5 * SIGMA; // радиус обрезания потенциала
const double RCUT2 = RCUT * RCUT;
const double SR = SIGMA / RCUT; 
const double SR2 = SR * SR;
const double SR6 = SR2 * SR2 * SR2;
const double UCUT = 4 * EPS * (SR6 * SR6 - SR6);

const double ACRIST = 1.0; //длина ребра элементарной ячейки


//размеры системы по осям координат
const double LX = NUMCRIST_X * ACRIST; 
const double LY = NUMCRIST_Y * ACRIST;
const double LZ = NUMCRIST_Z * ACRIST;

// LX = LY = LZ = <ЗНАЧЕНИЕ>;
const double VOLUME = LX * LY * LZ; //объем системы

