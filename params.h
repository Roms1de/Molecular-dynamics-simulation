#include <stdio.h>

#pragma once

//Число шагов
const int NSTEPS = 10;
const int LASTSTEP = NSTEPS - 1;

//число элементарных ячеек (кристаллов) по осям координат
const int NUMCRIST_X = 2; 
const int NUMCRIST_Y = 2; 
const int NUMCRIST_Z = 2;

// число частиц, первый конкретный расчет для двух частиц
const int NUMBERPARTICLES = 2;

//шаг интегрирования разностной схемы
const double STEP = 0.002;