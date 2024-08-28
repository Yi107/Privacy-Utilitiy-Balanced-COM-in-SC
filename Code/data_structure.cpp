#include "stdafx.h"

extern map <float, int> tmp_p;
extern float estimate_value_sp;
extern int total2;


int cmp_worker_arrive(const Worker &a, const Worker &b)
{
	return a.arrive < b.arrive;
}

int cmp_task_arrive(const Task &a, const Task &b)
{
	return a.arrive < b.arrive;
}

int cmp_task_value(const Task &a, const Task &b)
{
	return a.value > b.value;
}

int cmp_task_id(const Task &a, const Task &b)
{
	return a.id < b.id;
}

