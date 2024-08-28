#pragma once
#ifndef DATA_STRUCTURE_H
#define DATA_STRUCTURE_H

#include<vector>
#include<set>
using namespace std;

struct Worker
{
	int id;
	float x, y, radius;
	float per_x, per_y;
	int arrive; 
	int available; 
	vector<float> history_value;
	map<string, int> num_avgPrice;
	int type;
	

	bool operator == (const Worker b) const
	{
		return id == b.id;
	}

	bool operator < (const Worker b) const
	{
		return arrive < b.arrive;
	}

};

struct Task
{
	int id; 
	float end_x, end_y;
	float x, y;
	float value;
	float distance;
	int arrive;
	int endtime;
	int type;
};

const float alpha = 0.1, eta = 0.2, delta = 0.01;

multiset<Worker>::iterator find_my_worker(multiset<Worker> *workers, Task *task);
multiset<Worker>::iterator find_other_worker(multiset<Worker> *workers, Task *task, float value);
multiset<Worker>::iterator find_other_worker_with_bound(multiset<Worker> *workers, Task *task, float value);
int cmp_worker_arrive(const Worker &a, const Worker &b);
int cmp_task_arrive(const Task &a, const Task &b);
int cmp_task_value(const Task &a, const Task &b);
int cmp_task_id(const Task &a, const Task &b);
void find_my_workers(multiset<Worker> *workers, Task *task, vector<multiset<Worker>::iterator> *ans_it);
void find_my_workers(multiset<Worker> *workers, Task *task, vector<Worker> *ans_workers);
void find_my_workers(vector<Worker> *workers, Task *task, vector<Worker> *ans_workers);
int sub_size_serve_task(vector<Worker> *sub_worker, Task *task, float value);
pair<float, int> find_ac_ratio(vector<Worker> *workers, Task *task, float value);

#endif // !DATA_STRUCTURE_H
