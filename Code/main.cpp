#include "stdafx.h"

using namespace std;

int maxtime = 0;

void SplitString(string& s, vector<float>& v, const string& c)
{
	std::string::size_type pos1, pos2;
	pos2 = s.find(c);
	pos1 = 0;
	while (std::string::npos != pos2)
	{
		v.push_back(atof((s.substr(pos1, pos2 - pos1)).c_str()));

		pos1 = pos2 + c.size();
		pos2 = s.find(c, pos1);
	}
	if (pos1 != s.length())
		v.push_back(atof((s.substr(pos1)).c_str()));
}

void SplitStringS(string& s, vector<string>& v, const string& c)
{
	std::string::size_type pos1, pos2;
	pos2 = s.find(c);
	pos1 = 0;
	while (std::string::npos != pos2)
	{
		v.push_back(s.substr(pos1, pos2 - pos1).c_str());
		pos1 = pos2 + c.size();
		pos2 = s.find(c, pos1);
	}
	if (pos1 != s.length())
	{
		v.push_back(s.substr(pos1).c_str());
	}
}

void SplitStringInt(string& s, vector<int>& v, const string& c)
{
	std::string::size_type pos1, pos2;
	pos2 = s.find(c);
	pos1 = 0;
	while (std::string::npos != pos2)
	{
		v.push_back(atoi(s.substr(pos1, pos2 - pos1).c_str()));
		pos1 = pos2 + c.size();
		pos2 = s.find(c, pos1);
	}
	if (pos1 != s.length())
	{
		v.push_back(atoi(s.substr(pos1).c_str()));
	}
}

void read_file(string worker_filename, string task_filename, vector<Worker>* all_workers, vector<Task>* tasks, float self_radius)
{
	cout << "reading files .... " << endl;
	ifstream fin_w(worker_filename);
	
	string id;
	int type, arrive, endtime;
	float x, y, radius;
	float value;
	float end_x, end_y;
	float totdis, totvalue;
	string history;
	float d;
	string avgp, avgpn, value_history;

	while (fin_w >> id >> type >> radius >> arrive >> x >> y >> value_history >> avgp >> avgpn)
	{

		vector<string> avgPrice;
		vector<int> avgPNum;

		map<string, int> avP;
		Worker worker;
		worker.type = type;
		worker.arrive = arrive;
		worker.x = x;
		worker.y = y;
		worker.radius = self_radius;
		SplitString(value_history, worker.history_value, ",");
		SplitStringS(avgp, avgPrice, ",");
		SplitStringInt(avgpn, avgPNum, ",");
		for (int i = 0; i < avgPrice.size(); i++)
		{
			avP.insert(make_pair(avgPrice[i], avgPNum[i]));
		}
		worker.num_avgPrice = avP;
		SplitString(value_history, worker.history_value, ",");
		all_workers->push_back(worker);

	}

	ifstream fin_t(task_filename);
	while (fin_t >> id >> type >> arrive >> endtime >> x >> y >> end_x >> end_y >> d >> value)
	{
		Task task;
		task.type = type;
		task.arrive = arrive;
		task.endtime = endtime;
		task.x = x;
		task.y = y;
		task.distance = d;
		task.end_x = end_x;
		task.end_y = end_y;
		task.value = value;
		tasks->push_back(task);
	}

	sort(all_workers->begin(), all_workers->end(), cmp_worker_arrive);
	sort(tasks->begin(), tasks->end(), cmp_task_arrive);
	maxtime = tasks->rbegin()->arrive;
	for (int i = 0; i < tasks->size(); i++)
	{
		(*tasks)[i].id = i;
	}

	for (int i = 0; i < all_workers->size(); i++)
	{
		(*all_workers)[i].id = i;
	}
	cout << "reading complete ..." << endl;
}




void get_workers(vector<Worker>* all_workers, map<int, multiset<Worker>>* workers)
{
	for (int i = 0; i < (*all_workers).size(); i++)
	{
		Worker worker = (*all_workers)[i];
		int type = worker.type;
		map<int, multiset<Worker>>::iterator it = workers->find(type);
		if (it != workers->end())
		{
			it->second.insert(worker);
		}
		else {
			multiset<Worker> nWorker;
			nWorker.insert(worker);
			workers->insert(make_pair(type, nWorker));
		}
	}

}


void get_otherworkers(vector<Worker>* all_workers, map<int, multiset<Worker>>* myWorkers, map<int, multiset<Worker>>* otherWorkers, int type)
{
	for (int i = 0; i < (*all_workers).size(); i++)
	{
		Worker worker = (*all_workers)[i];
		int type = worker.type;
		map<int, multiset<Worker>>::iterator it = myWorkers->find(type);
		if (it != myWorkers->end())
		{
			it->second.insert(worker);
		}
		else {
			multiset<Worker> nWorker;
			nWorker.insert(worker);
			myWorkers->insert(make_pair(type, nWorker));
		}
	}
	for (int i = 1; i <= type; i++)
	{
		multiset<Worker> otherWorker;

		for (int j = 1; j <= type; j++)
		{
			if (i != j)
			{
				multiset<Worker> subworker;
				subworker = myWorkers->find(j)->second;
				otherWorker.insert(subworker.begin(), subworker.end());
			}
		}
		otherWorkers->insert(make_pair(i, otherWorker));
	}
}
void getWorkers(vector<Worker>* all_workers, map<int, multiset<Worker>*>* workers)
{
	vector<Worker>::iterator iw = all_workers->begin();
	for (; iw != all_workers->end(); iw++)
	{
		int type = iw->type;
		map<int, multiset<Worker>*>::iterator im = workers->find(type);
		if (im != workers->end())
		{
			im->second->insert(*iw);
		}
			
	}
	cout << "";
}



void printresult(map<int, vector<int>> matchNum, map<int, vector<double>> rev, string outputfile)
{
	ofstream fout(outputfile, ios::app);
	map<int, vector<int>>::iterator matchit = matchNum.begin();
	for (; matchit != matchNum.end(); matchit++)
	{
		float ac = (float)matchit->second[1] / (float)matchit->second[2];
		cout << "Type " << matchit->first << " Total Matching " << matchit->second[0] << " Out Matching " << matchit->second[2] << " Real Out " << matchit->second[1] << endl;
		cout << "Type " << matchit->first << " AC: " << ac << endl;
		fout << "Type " << matchit->first << " Total Matching " << matchit->second[0] << " Out Matching " << matchit->second[2] << " Real Out " << matchit->second[1] << endl;
		fout << "Type " << matchit->first << " AC: " << ac << endl;
		fout <<  matchit->second[0] << "\t" << matchit->second[2] << "\t" << matchit->second[1] << endl;
	}

	fout << "Total Matching Num: ";
	matchit = matchNum.begin();
	for (; matchit != matchNum.end(); matchit++)
	{
		fout << matchit->second[0] << "\t";
	}
	fout << endl;
	fout << "Total Cooperative Matching Num: ";
	matchit = matchNum.begin();
	for (; matchit != matchNum.end(); matchit++)
	{
		fout << matchit->second[1] << "\t";
	}


	fout << endl;
	fout << "AC: ";
	matchit = matchNum.begin();
	for (; matchit != matchNum.end(); matchit++)
	{
		float ac = (float)matchit->second[1] / (float)matchit->second[2];
		fout << ac << "\t";
	}
	fout << endl;
	map<int, vector<double>>::iterator matchrev = rev.begin();
	for (; matchrev != rev.end(); matchrev++)
	{
		string val = to_string(matchrev->second[0]);
		cout << "Type " << matchrev->first << " Total Rev " << val << " Give to other Out Rev: " << matchrev->second[1] <<" Finish other r num:"<< matchrev->second[2] <<" Finish other r rev: "<< matchrev->second[3] << endl;
		fout << "Type " << matchrev->first << " Total Rev " << val << "Give to other Out Rev" << matchrev->second[1] << " Finish other r num:" << matchrev->second[2] << " Finish other r rev: " << matchrev->second[3] << endl;
		
	}
	matchrev = rev.begin();
	for (; matchrev != rev.end(); matchrev++)
	{
		string val = to_string(matchrev->second[0]);
		fout << val <<"\t";

	}
	fout << endl;
}

void get_timerange(vector<Task>* tasks, vector<int>* timerange) {
	vector<Task>::iterator ite = tasks->begin();
	int min_time = ite->arrive, max_time = ite->arrive;
	for (; ite != tasks->end(); ite++) {
		if (ite->arrive < min_time) {
			min_time = ite->arrive;
		}
		else if (ite->arrive > max_time) {
			max_time = ite->arrive;
		}
	}
	int timescale = max_time - min_time;
	int pertime = timescale / 24;
	int timetemp = min_time + 4 * pertime;
	timerange->push_back(timetemp);
	timetemp = timetemp + 3 * pertime;
	timerange->push_back(timetemp);
	timetemp = timetemp + 3 * pertime;
	timerange->push_back(timetemp);
	timetemp = timetemp + 7 * pertime;
	timerange->push_back(timetemp);
	timetemp = timetemp + 3 * pertime;
	timerange->push_back(timetemp);
	timetemp = timetemp + 2 * pertime;
	timerange->push_back(timetemp);
}

void cal_rev(vector<Worker>* all_workers, vector<Task>* tasks, vector<tuple<int, int, float> >* M, float* rev1, float* rev2, int* task_1, int* task_2)
{
	*rev1 = 0.0, * rev2 = 0.0;
	*task_1 = 0, * task_2 = 0;
	for (int i = 0; i < M->size(); i++)
	{
		int worker_id = std::get<0>((*M)[i]);
		int task_id = std::get<1>((*M)[i]);
		float value = std::get<2>((*M)[i]);
		if ((*tasks)[task_id].type == 1)
		{
			if ((*all_workers)[worker_id].type == 1)
			{
				*rev1 += value;
			}
			else
			{
				*rev2 += (*tasks)[task_id].value - value;
				*rev1 += value;
			}
			*task_1 += 1;
		}
		else
		{
			if ((*all_workers)[worker_id].type == 2)
			{
				*rev2 += value;
			}
			else
			{
				*rev1 += (*tasks)[task_id].value - value;
				*rev2 += value;
			}
			*task_2 += 1;
		}
	}

}

void get_tasksrange(vector<Task>* tasks, vector<int> timerange, map<int, vector<Task>>* tasksrange) { // 根据时间段将所有的任务化为3部分
	vector<Task>::iterator ite = tasks->begin();
	for (; ite != tasks->end(); ite++) {
		Task tmp = *ite;
		if (ite->arrive < timerange[0] || ite->arrive > timerange[5]) {
			map<int, vector<Task>>::iterator it = tasksrange->find(0);
			if (it != tasksrange->end()) {
				it->second.push_back(tmp);
			}
			else {
				vector<Task> ntask;
				ntask.push_back(tmp);
				tasksrange->insert(make_pair(0, ntask));
			}
		}
		else if ((ite->arrive > timerange[1] && ite->arrive < timerange[2]) || (ite->arrive > timerange[3] && ite->arrive < timerange[4])) {
			map<int, vector<Task>>::iterator it = tasksrange->find(1);
			if (it != tasksrange->end()) {
				it->second.push_back(tmp);
			}
			else {
				vector<Task> ntask;
				ntask.push_back(tmp);
				tasksrange->insert(make_pair(1, ntask));
			}
		}
		else {
			map<int, vector<Task>>::iterator it = tasksrange->find(2);
			if (it != tasksrange->end()) {
				it->second.push_back(tmp);
			}
			else {
				vector<Task> ntask;
				ntask.push_back(tmp);
				tasksrange->insert(make_pair(2, ntask));
			}
		}
	}
}

int main(int argc, char* argv[])
{
	string workerfile = "new_syn/worker-500";
	vector<float> radius;
	radius.push_back(0.5);
	radius.push_back(1);
	radius.push_back(1.5);
	radius.push_back(2);
	vector<int> tasknum;
	tasknum.push_back(500);
	tasknum.push_back(1000);
	tasknum.push_back(2500);
	tasknum.push_back(5000);
	tasknum.push_back(10000);
	tasknum.push_back(20000);
	tasknum.push_back(50000);
	vector<float> eset, geset;
	eset.push_back(0.1);
	eset.push_back(0.4);
	eset.push_back(0.7);
	eset.push_back(1);
	//geset.push_back(0);
	geset.push_back(0.1);
	geset.push_back(0.4);
	geset.push_back(0.7);
	geset.push_back(1);
	vector<int> workernum;
	workernum.push_back(100);
	workernum.push_back(200);
	workernum.push_back(500);
	workernum.push_back(1000);
	workernum.push_back(2500);
	workernum.push_back(5000);
	string requestfile = "new_syn/request-2500_new";

	float geoe = 0.7;
	float este = 0.7;
	string outputfile = "new_syn/new_R/testmemo";

	vector<Worker> all_workers; 
	vector<Task> tasks; 
	vector<tuple<int, int, float> > M;
	int typenum = 3; // 公司数量
	map<int, Worker> workerSet;
	string output;
	map<int, vector<int>> matchNum;
	map<int, vector<double>> rev;

	for (int i = 0; i < workernum.size(); i++) {
		
		all_workers.clear();
		tasks.clear();
		clock_t start, end;
		string worker_file, task_file;
		workerfile = "new_syn/worker-";
		workerfile = workerfile + to_string(workernum[i]);
		worker_file = workerfile;
		task_file = requestfile;
		ofstream fout(outputfile, ios::app);
		read_file(worker_file, task_file, &all_workers, &tasks, 1);
		map<int, multiset<Worker>*> workers;
		multiset<Worker> w1, w2, w3;
		workers.insert(make_pair(1, &w1));
		workers.insert(make_pair(2, &w2));
		workers.insert(make_pair(3, &w3));
		getWorkers(&all_workers, &workers);
		cout << "Geo epsilon: " << geoe << "Pricing epsilon: " << este << endl;
		fout << " Workers: " << 500;
		fout << " Tasks: " << 2500;
		fout << " Radius: " << 1;
		fout << " Geo epsilon: " << geoe << "Pricing epsilon: " << este << endl;
		start = clock();
		basedMatch(&workers, &tasks, &M, &matchNum, &rev, este, geoe);
		end = clock();
		double duration;
		duration = ((double)(end - start) / CLK_TCK);
		duration /= 60;
		printresult(matchNum, rev, outputfile);
		fout << "Time：" << duration << "min" << endl << endl;
		workers.clear();
		M.clear();
		matchNum.clear();
		rev.clear();
	}


}
