#include <stdio.h>

#pragma once

//Число шагов
const int NSTEPS = 2000;
const int LASTSTEP = NSTEPS - 1;

//число элементарных ячеек (кристаллов) по осям координат
const int NUMCRIST_X = 4; 
const int NUMCRIST_Y = 4; 
const int NUMCRIST_Z = 4;

// число частиц, первый конкретный расчет для двух частиц
const int NUMBERPARTICLES = 4;

//шаг интегрирования разностной схемы
const double STEP = 0.002;