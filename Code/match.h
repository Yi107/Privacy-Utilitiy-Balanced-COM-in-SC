#pragma once
#include"stdafx.h";
void basedMatch(map<int, multiset<Worker>*>* workers, vector<Task>* tasks, vector<tuple<int, int, float> >* M, map<int, vector<int>>* matchNum, map<int, vector<double>>* rev, float epsilon, float geoepsilon);