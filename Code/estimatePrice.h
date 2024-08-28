#pragma once
#include"stdafx.h"
void caltotMin(map<int, float> *minAvg, float distance, int *mintype, float *minpri);
void pri_calMinForEvery(map<int, multiset<Worker>> *avaWorker, map<int, float> *minAvg, float priepsilon);
void calMinForEvery(map<int, multiset<Worker>> *avaWorker, map<int, float> *minAvg);
//void calAverForW(map<int, multiset<Worker>> *avaWorker, map<int, float> *avgForEP);
void calAvgMin(map<int, float> *avgForEP, int *mintype, float *minpri);