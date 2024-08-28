#define _USE_MATH_DEFINES
#include "stdafx.h";
#include <math.h>
static int rejectByPrice = 0;
static int rejectByNoAv = 0;
static int rejectFinal = 0;
static int priceFailed = 0;
static int probZero = 0;
static int rejectBydis = 0;
static int rejectByProb = 0;
static double maxduration=0;
static double maxCalAvaw = 0;
map<int, vector<int>> finishedLocalWorkerTasks;
map<int, vector<int>> finishedOuterWorkerTasks;

void find_avai_worker(multiset<Worker>* workers, Task* task, map<int, multiset<Worker>> *avaiWorker)
{
	float task_x = task->x, task_y = task->y;
	for (multiset<Worker>::iterator it = workers->begin(); it != workers->end(); it++)
	{
		if (it->arrive > task->arrive)
		{
			break;
		}
		float worker_x = it->x, worker_y = it->y;
		float dis = get_distance(task_x, task_y, worker_x, worker_y);
		if (dis <= it->radius)
		{
			map<int, multiset<Worker>>::iterator itw = avaiWorker->find(it->type);
			if (itw != avaiWorker->end())
			{
				itw->second.insert(*it);
			}
			else
			{
				multiset<Worker> w;
				w.insert(*it);
				avaiWorker->insert(make_pair(it->type, w));
			}

		}

	}	
}

void find_extended_avai_worker(multiset<Worker>* workers, Task* task, float epsilon, int angle, map<int, multiset<Worker>> *avaiWorker, map<int, multiset<Worker>> *perturbedWorker, map<int,multiset<Worker>> *usedCalPrice)
{
	float task_x = task->x, task_y = task->y;
	double duration = 0;
	for (multiset<Worker>::iterator it = workers->begin(); it != workers->end(); it++)
	{
		Worker perturbedw = *it;
		double perturbed_lat = 0, perturbed_lon = 0;
		avaWLocationPerturb(perturbedw.x, perturbedw.y, angle, &perturbed_lat, &perturbed_lon);
		if (maxCalAvaw < duration)
		{
			maxCalAvaw = duration;
		}
		perturbedw.x = perturbed_lat;
		perturbedw.y = perturbed_lon;
		if (it->arrive > task->arrive)
		{
			break;
		}
		float worker_x = it->x, worker_y = it->y;
		float dis = get_distance(task_x, task_y, worker_x, worker_y);
		float extend_r = 0;
		extendRadius(it->radius, epsilon, &extend_r);
		if (maxCalAvaw < duration)
		{
			maxCalAvaw = duration;
		}
	
		if (dis <= extend_r)
		{
			map<int, multiset<Worker>>::iterator itw = avaiWorker->find(it->type);
			map<int, multiset<Worker>>::iterator ipw = perturbedWorker->find(it->type);
			if (itw != avaiWorker->end())
			{
				itw->second.insert(*it);
				ipw->second.insert(perturbedw);
			}
			else
			{
				multiset<Worker> w1,w2;
				w1.insert(*it);
				w2.insert(perturbedw);
				avaiWorker->insert(make_pair(it->type, w1));
				perturbedWorker->insert(make_pair(it->type, w2));
			}
			if (dis <= it->radius) {
				map<int, multiset<Worker>>::iterator icp = usedCalPrice->find(it->type);
				if (icp != usedCalPrice->end())
					{
						icp->second.insert(*it);
					}
					else
					{
						multiset<Worker> w;
						w.insert(*it);
						usedCalPrice->insert(make_pair(it->type, w));
					}
			}
		}
	}
}


multiset<Worker>::iterator find_my_worker(multiset<Worker>* workers, Task* task)
{
	
	multiset<Worker>::iterator ans_it = workers->end();
	float ans_dis = 1000000000.0;
	float task_x = task->x, task_y = task->y;
	for (multiset<Worker>::iterator it = workers->begin(); it != workers->end(); it++)
	{
		Worker worker = *it;
		if (worker.arrive > task->arrive) 
		{
			break;
		}
		float worker_x = worker.x, worker_y = worker.y;
		float dis = get_distance(task_x, task_y, worker_x, worker_y);
		if (dis <= worker.radius && dis < ans_dis)
		{
			ans_it = it;
			ans_dis = dis;
			break;
		}
	}
	return ans_it;
}


Worker find_other_worker(multiset<Worker> workers, float value, Task* task) {
	multiset<Worker>::iterator ans_it = workers.end();
	float ans_dis = 1000000000.0;
	float ans_p = 0;
	float task_x = task->x, task_y = task->y;
	int cum = 0;
	
	for (multiset<Worker>::iterator it = workers.begin(); it != workers.end(); it++)
	{
		Worker worker = *it;
		if (worker.arrive > task->arrive) 
		{
			break;
		}
		float worker_x = worker.x, worker_y = worker.y;
		float dis = get_distance(task_x, task_y, worker_x, worker_y);
	
		if (dis <= worker.radius)
		{
			cum++;
			vector<float> history_value = worker.history_value;
			int sum = 0;
			for (int i = 0; i < history_value.size(); i++)
			{
				if (history_value[i] <= value)
				{
					sum++;
				}
			}
			float p = rand() / double(RAND_MAX);
			float sp = 1.0 * sum / history_value.size();
			
			if (p > sp)
			{
				rejectByProb++;
			}
			if (p <= sp && dis < ans_dis)
			{
				ans_it = it;
				ans_dis = dis;
				
			}
		}
	}
	if (ans_it != workers.end())
	{
		return *ans_it;
	}
	else
	{
		Worker w;
		w.id = -1;
		return w;
	}
}

void calculatePricePlus(map<string, int>* totalnum, float epsilon, float* outprice)
{
	*outprice = 0;
	map<string, int> basePro;
	map<string, float> meanMap;
	map<string, int>::iterator itot = totalnum->begin();
	int totNum = 0;
	for (; itot != totalnum->end(); itot++)
	{
		totNum += itot->second;
		basePro.insert(make_pair(itot->first, totNum));
	}
	itot = basePro.begin();
	for (; itot != basePro.end(); itot++)
	{
		float avp = atof(itot->first.c_str());
		float mean = avp * (1 - 1.0 * itot->second / totNum);
		meanMap.insert(make_pair(itot->first, mean));
	}
	map<string, float> expFunc;
	float totalFunc = 0;
	*outprice = 0;
	map<string, float>::iterator ti = meanMap.begin();
	float maxAvg = atof(meanMap.rbegin()->first.c_str());
	for (; ti != meanMap.end(); ti++)
	{
		float num = ti->second;
		float expNum = 0;
		if (epsilon != 0)
		{
			expNum = exp(epsilon * num / (2 * 0.05 * maxAvg));
		}
		else {
			expNum = num;
		}
		totalFunc += expNum;
		expFunc.insert(make_pair(ti->first, expNum));
	}
	map<string, float> probExp;
	map<string, float>::iterator expi = expFunc.begin();
	float tempPro = 0;
	for (; expi != expFunc.end(); expi++)
	{
		float prob = expi->second / totalFunc;
		tempPro += prob;
		probExp.insert(make_pair(expi->first, tempPro));
	}
	map<string, float>::iterator proi = probExp.begin();
	float randPro = rand() % (NN + 1) / (float)(NN + 1);
	float minP = 0, maxP;
	for (; proi != probExp.end(); proi++)
	{
		maxP = proi->second;
		if (randPro >= minP && randPro < maxP)
		{
			string value = proi->first;
			*outprice = atof(value.c_str());
			return;
		}
		minP = maxP;
	}
}

void calEveryPlus(map<int, multiset<Worker>>* ava_worker, float dis, float epsilon, map<int, float>* priceSet)
{
	int mintype = -1;
	float minPrice = 10000;
	map<int, multiset<Worker>>::iterator iw = ava_worker->begin();
	if (dis <= 1)
	{
		minPrice = 9;
		priceSet->insert(make_pair(ava_worker->begin()->first, minPrice));
		priceSet->insert(make_pair(ava_worker->rbegin()->first, minPrice));

	}
	else
	{
		for (; iw != ava_worker->end(); iw++)
		{
			minPrice = 10000;
			multiset<Worker>::iterator w = iw->second.begin();
			for (; w != iw->second.end(); w++)
			{
				float outprice;
				map<string, int> per_w = w->num_avgPrice;
				calculatePricePlus(&per_w, epsilon, &outprice);
				if (outprice != 0 && outprice < minPrice)
				{
					minPrice = outprice;
				}
			}
			if (minPrice != 10000)
			{
				priceSet->insert(make_pair(iw->first, (minPrice+0.3) * dis + 9));
			}
			
		}

	}


}



void calMinForEveryPlus(map<int, multiset<Worker>>* ava_worker, float dis, int* mintype, float* minPrice, float epsilon)
{
	*minPrice = 10000;
	*mintype = -1;
	map<int, multiset<Worker>>::iterator iw = ava_worker->begin();
	if (dis<=1)
	{
		*minPrice = 9;
		int randtype = rand() % 2;
		if (randtype == 0)
		{
			*mintype = ava_worker->begin()->first;
		}
		else
		{
			*mintype = ava_worker->rbegin()->first;
		}

	}
	else {
		for (; iw != ava_worker->end(); iw++)
		{
			multiset<Worker>::iterator w = iw->second.begin();
			for (; w != iw->second.end(); w++)
			{
				float outprice;
				map<string, int> per_w = w->num_avgPrice;
				calculatePricePlus(&per_w, epsilon, &outprice);
				if (outprice != 0 && outprice < *minPrice)
				{
					*minPrice = outprice;
					*mintype = w->type;
				}
			}
		}
		*minPrice = *minPrice * dis + 9;
	}
	
}

void calProb(float r_lat, float r_lng, float radius, int angle, float perw_lat, float perw_lng, float *prob)
{
	float a, b, c;
	float ax, bx;
	double r_x = 0, r_y = 0, w_x = 0,w_y=0;
	ZtGeographyCoordinateTransform ztGCT;
	ztGCT.BL2XY(r_lat, r_lng, r_x, r_y);
	r_x = r_x / 1000;
	r_y = r_y / 1000;
	ZtGeographyCoordinateTransform ztGCT1;
	ztGCT1.BL2XY(perw_lat, perw_lng, w_x, w_y);
	w_x = w_x / 1000;
	w_y = w_y / 1000;
	ax = tan((float)angle * M_PI / 180.0);
	bx = w_y - tan((float)angle * M_PI / 180.0)*w_x;
	a = 1 + pow(ax, 2);
	b = -2 * r_x + 2 * ax*(bx - r_y);
	c = pow(r_x, 2) + pow(((bx - r_y)), 2) - pow(radius,2);
	float delta = pow(b, 2) - 4 * a*c;
	if (delta > 0)
	{
		float x1 = (-b + sqrt(delta)) / (2 * a);
		float x2 = (-b - sqrt(delta)) / (2 * a);
		float d1 = cos((float)angle * M_PI / 180.0)*(x1 - w_x);
		float d2 = cos((float)angle * M_PI / 180.0)*(x2 - w_x);
		float c1, c2;
		if (d1 > 0)
		{
			c1 = 1 - exp(-d1 * 0.7);
		}
		else
		{
			c1 = exp(0.7*d1) - 1;
		}
		if (d2 > 0)
		{
			c2 = 1 - exp(-d2 * 0.7);
		}
		else
		{
			c2 = exp(0.7*d2) - 1;
		}
		*prob = 0.5*abs((c1 - c2));

	}
		
	
	if (delta <= 0)
	{
		*prob = 0;
	}
	
}


void Tripsolve(Task* task, multiset<Worker>* myWorkers, map<int, multiset<Worker>*> *otherWorkers, vector<tuple<int, int, float> >* M, map<int, vector<int>>* matchNum, map<int, vector<double>>* rev, float epsilon, float geoepsilon)//matchNum: id, 此id下完成的task数量,外派出去的数量,外派接受的数量
{
	clock_t start, end, start1, end1, start2, end2;
	int type = task->type;

	multiset<Worker>::iterator my_aviWorker;	
	if (task->value > exp(3))
	{
		my_aviWorker = find_my_worker(myWorkers, task);
		if (my_aviWorker != myWorkers->end())
		{
			map<int, vector<int>>::iterator findLocal = finishedLocalWorkerTasks.find(my_aviWorker->id);
			if (findLocal != finishedLocalWorkerTasks.end())
			{

				findLocal->second.push_back(task->id);
			}
			else {
				vector<int> temp;
				finishedLocalWorkerTasks.insert(make_pair(my_aviWorker->id, temp));
			}
			Worker oldworker = *my_aviWorker;
			Worker worker = oldworker;
			worker.arrive = task->endtime + 1;
			worker.history_value.push_back(task->value);
			worker.x = task->end_x;
			worker.y = task->end_y;
			M->push_back(make_tuple(worker.id, task->id, task->value));			
			myWorkers->erase(my_aviWorker);
			myWorkers->insert(worker);
			map<int, vector<int>>::iterator numI = matchNum->find(type);
			if (numI != matchNum->end())
			{
				int num = numI->second[0];
				num++;
				numI->second[0] = num;
			}
			else {
				vector<int> nu;
				nu.push_back(1);
				nu.push_back(0);
				nu.push_back(0);
				matchNum->insert(make_pair(type, nu));
			}

			map<int, vector<double>>::iterator revI = rev->find(type);
			if (revI != rev->end())
			{
				float revn = revI->second[0];
				revn = revn + task->value;
				revI->second[0] = revn;
			}
			else
			{
				vector<double> re;
				re.push_back(task->value);
				re.push_back(0);
				re.push_back(0);
				re.push_back(0);
				rev->insert(make_pair(type, re));
			}

		}
		
	}
	else
	{	
	
		Task pri_task;
		string finalPrice = "FinalPrice";
		string outerRDetails = "outerRDetails";
		string durationFile = "maxDuration";
		start = clock();
		int angle = 0;
		map<int, multiset<Worker>> ava_worker;
		map<int, multiset<Worker>> per_worker;
		map<int, multiset<Worker>> cal_PriceWorker;
		double duration = 0;
		if (geoepsilon != 0) {
			pri_task = *task;
			float x, y, prip, newx, newy;
			float noise1=0, noise2=0;
			int randomA = 0;
			x = task->x;
			y = task->y;
			prip = task->value;
			double lat=0, lon=0, ht;
			start = clock();
			locationPerturb(x,y, geoepsilon, myWorkers->begin()->radius, &lat, &lon, &angle);
			end = clock();
			pri_task.x = lat;
			pri_task.y = lon;
			map<int, multiset<Worker>*>::iterator oti = otherWorkers->begin();
			for (;oti!=otherWorkers->end();oti++)
			{
				map<int, multiset<Worker>> avaiWorker;
				map<int, multiset<Worker>> perWorker;
				map<int, multiset<Worker>> calPriceWorker;
				start = clock();
				find_extended_avai_worker(oti->second, &pri_task, geoepsilon,angle, &avaiWorker,&perWorker, &calPriceWorker);
				end = clock();
				duration = ((double)(end - start) / CLK_TCK);
				if (maxduration < duration)
				{
					maxduration = duration;
				}
				
				if (!avaiWorker.empty())
				{
					ava_worker.insert(*avaiWorker.begin());
					per_worker.insert(*perWorker.begin());
				}	
				
				if (!calPriceWorker.empty())
				{
					cal_PriceWorker.insert(*calPriceWorker.begin());
				}
			}

		}
		else {			
			map<int, multiset<Worker>*>::iterator oti = otherWorkers->begin();
			for (; oti != otherWorkers->end(); oti++)
			{
				map<int, multiset<Worker>> avaiWorker;
				find_avai_worker(oti->second, task, &avaiWorker);
				if (!avaiWorker.empty())
				{
					ava_worker.insert(*avaiWorker.begin());
				}
			}
		}


		map<int, vector<int>>::iterator numI = matchNum->find(type);
		if (numI != matchNum->end())
		{
			int num = numI->second[2];
			num++;
			numI->second[2] = num;
		}
		else {
			vector<int> nu;
			nu.push_back(0);
			nu.push_back(0);
			nu.push_back(1);
			matchNum->insert(make_pair(type, nu));
		}
		if (ava_worker.size() == 0)
		{
			rejectByNoAv++;
			return;
		}


		map<int, float> priceSet;
		int mintype;
		float minpri = 0;
		start = clock();
		if (cal_PriceWorker.size() != 0)
		{
			calEveryPlus(&cal_PriceWorker, task->distance, epsilon, &priceSet);
		}
		else {
			calEveryPlus(&ava_worker, task->distance, epsilon, &priceSet);
		}
		end = clock();
		
		if (priceSet.size() == 0)
		{
			priceFailed++;
			return;
		}
		map<int, float>::iterator ip = priceSet.begin();
		vector<int> avaPlatform;
		for (; ip != priceSet.end(); ip++)
		{
			if (ip->second <= task->value)
			{
				avaPlatform.push_back(ip->first);
			}
		}
		if (avaPlatform.size() == 0)//没有平台价格小于订单价格
		{
			rejectByPrice++;
			return;
		}

		string probfile = "eachprob";
		int selectedPlatform = -1;
		float maxPro = 0; 
		for (int i = 0; i < avaPlatform.size(); i++)
		{
			float tempMax = 0;
			start = clock();
			multiset<Worker> perturbedWorker = per_worker.find(avaPlatform[i])->second;
			if (perturbedWorker.size() != 0)
			{
				multiset<Worker>::iterator ppw = perturbedWorker.begin();
				for (; ppw != perturbedWorker.end(); ppw++)
				{
					float prob;
					calProb(task->x, task->y, ppw->radius, angle,  ppw->x,  ppw->y, &prob);
					if (prob > maxPro)
					{
						maxPro = prob;
						selectedPlatform = ppw->type;
					}
					if (prob > tempMax)
					{
						tempMax = prob;
					}
				}	
			}
			end = clock();
		}
		
		if (maxPro == 0)
		{
			probZero++;
			return;
		}


		multiset<Worker> selectedWorker = ava_worker.find(selectedPlatform)->second; 
		float outerPrice = priceSet.find(selectedPlatform)->second;

		Worker other_aviWorker;
		start = clock();

		other_aviWorker = find_other_worker(selectedWorker, outerPrice, &pri_task);
		
		if (other_aviWorker.id == -1)
		{
			rejectFinal++;
		}
		int outworkertype;
		float dis = get_distance(task->x, task->y, other_aviWorker.x, other_aviWorker.y);
		if (other_aviWorker.id != -1 & dis <= other_aviWorker.radius)
		{
			map<int, vector<int>>::iterator findOuter = finishedOuterWorkerTasks.find(other_aviWorker.id);
			if (findOuter != finishedOuterWorkerTasks.end())
			{

				findOuter->second.push_back(task->id);
			}
			else {
				vector<int> temp;
				finishedOuterWorkerTasks.insert(make_pair(other_aviWorker.id, temp));
			}
			outworkertype = other_aviWorker.type;
			map<int, multiset<Worker>*>::iterator imi = otherWorkers->find(other_aviWorker.type);
			for (multiset<Worker>::iterator it = imi->second->begin(); it != imi->second->end(); it++)
			{
				if (it->id == other_aviWorker.id)
				{
					Worker sltWorker = other_aviWorker;
					sltWorker.arrive = task->endtime + 1;
					sltWorker.history_value.push_back(outerPrice);
					sltWorker.x = task->end_x;
					sltWorker.y = task->end_y;
					M->push_back(make_tuple(sltWorker.id, task->id, task->value - outerPrice));
					imi->second->erase(it);
					imi->second->insert(sltWorker);
					break;
				}			
			}
			

			map<int, vector<int>>::iterator numII = matchNum->find(type);
			if (numII != matchNum->end())
			{
				int num = numII->second[0];
				int othernum = numII->second[1];
				num++;
				othernum++;
				numII->second[0] = num;
				numII->second[1] = othernum;
			}
			else {
				vector<int> nu;
				nu.push_back(1);
				nu.push_back(1);
				nu.push_back(1);
				matchNum->insert(make_pair(type, nu));
			}

			map<int, vector<double>>::iterator revI = rev->find(type);
			if (revI != rev->end())
			{
				float revn = revI->second[0];
				float otherrev = revI->second[1];
				revn = revn + task->value - outerPrice;
				otherrev = otherrev + outerPrice;
				revI->second[0] = revn;
				revI->second[1] = otherrev;
	
			}
			else
			{
				vector<double> re;
				re.push_back(task->value - outerPrice);
				re.push_back(outerPrice);
				re.push_back(0);
				re.push_back(0);
				rev->insert(make_pair(type, re));
			}
			revI = rev->find(outworkertype);
			if (revI != rev->end())
			{
				float revn = revI->second[0];				
				revn = revn + outerPrice;
				revI->second[0] = revn;
				float outnum = revI->second[2];
				revI->second[2] += 1;
				float outrev = revI->second[3];
				revI->second[3] += outerPrice;

			}
			else
			{
				vector<double> re;
				//re.push_back(totPrice);
				re.push_back(outerPrice);
				re.push_back(0);
				re.push_back(1);
				re.push_back(outerPrice);
				rev->insert(make_pair(type, re));
			}

			
		}
		end = clock();
	}
}

void TripNoSelectionsolve(Task* task, multiset<Worker>* myWorkers, map<int, multiset<Worker>*> *otherWorkers, vector<tuple<int, int, float> >* M, map<int, vector<int>>* matchNum, map<int, vector<double>>* rev, float epsilon, float geoepsilon)//matchNum: id, 此id下完成的task数量,外派出去的数量,外派接受的数量
{
	clock_t start, end, start1, end1, start2, end2;
	int type = task->type;

	multiset<Worker>::iterator my_aviWorker;
	if (task->value > exp(3))
	{
		my_aviWorker = find_my_worker(myWorkers, task);
		if (my_aviWorker != myWorkers->end())
		{
			map<int, vector<int>>::iterator findLocal = finishedLocalWorkerTasks.find(my_aviWorker->id);
			if (findLocal != finishedLocalWorkerTasks.end())
			{

				findLocal->second.push_back(task->id);
			}
			else {
				vector<int> temp;
				finishedLocalWorkerTasks.insert(make_pair(my_aviWorker->id, temp));
			}
			Worker oldworker = *my_aviWorker;
			Worker worker = oldworker;
			worker.arrive = task->endtime + 1;
			worker.history_value.push_back(task->value);
			worker.x = task->end_x;
			worker.y = task->end_y;
			M->push_back(make_tuple(worker.id, task->id, task->value));
			myWorkers->erase(my_aviWorker);
			myWorkers->insert(worker);
			map<int, vector<int>>::iterator numI = matchNum->find(type);
			if (numI != matchNum->end())
			{
				int num = numI->second[0];
				num++;
				numI->second[0] = num;
			}
			else {
				vector<int> nu;
				nu.push_back(1);
				nu.push_back(0);
				nu.push_back(0);
				matchNum->insert(make_pair(type, nu));
			}

			map<int, vector<double>>::iterator revI = rev->find(type);
			if (revI != rev->end())
			{
				float revn = revI->second[0];
				revn = revn + task->value;
				revI->second[0] = revn;
			}
			else
			{
				vector<double> re;
				re.push_back(task->value);
				re.push_back(0);
				re.push_back(0);
				re.push_back(0);
				rev->insert(make_pair(type, re));
			}

		}
		
	}
	else
	{

		Task pri_task;
		string finalPrice = "FinalPrice";
		string outerRDetails = "outerRDetails";
		string durationFile = "maxDuration";
		start = clock();
		int angle = 0;
		map<int, multiset<Worker>> ava_worker;
		map<int, multiset<Worker>> per_worker;
		map<int, multiset<Worker>> cal_PriceWorker;
		double duration = 0;
		if (geoepsilon != 0) {
			
			pri_task = *task;
			float x, y, prip, newx, newy;
			float noise1 = 0, noise2 = 0;
			int randomA = 0;
			x = task->x;
			y = task->y;
			prip = task->value;
			double lat = 0, lon = 0, ht;
			start = clock();
			locationPerturb(x, y, geoepsilon, myWorkers->begin()->radius, &lat, &lon, &angle);
			end = clock();
			
			pri_task.x = lat;
			pri_task.y = lon;
			map<int, multiset<Worker>*>::iterator oti = otherWorkers->begin();
			for (; oti != otherWorkers->end(); oti++)
			{
				map<int, multiset<Worker>> avaiWorker;
				map<int, multiset<Worker>> perWorker;
				map<int, multiset<Worker>> calPriceWorker;
				start = clock();
				find_extended_avai_worker(oti->second, &pri_task, geoepsilon, angle, &avaiWorker, &perWorker, &calPriceWorker);
				end = clock();
				duration = ((double)(end - start) / CLK_TCK);
				if (maxduration < duration)
				{
					maxduration = duration;
				}

				if (!avaiWorker.empty())
				{
					ava_worker.insert(*avaiWorker.begin());
					per_worker.insert(*perWorker.begin());
				}

				if (!calPriceWorker.empty())
				{
					cal_PriceWorker.insert(*calPriceWorker.begin());
				}
			}

		}
		else {
			map<int, multiset<Worker>*>::iterator oti = otherWorkers->begin();
			for (; oti != otherWorkers->end(); oti++)
			{
				map<int, multiset<Worker>> avaiWorker;
				find_avai_worker(oti->second, task, &avaiWorker);
				if (!avaiWorker.empty())
				{
					ava_worker.insert(*avaiWorker.begin());
				}
				
			}

		}


		map<int, vector<int>>::iterator numI = matchNum->find(type);
		if (numI != matchNum->end())
		{
			int num = numI->second[2];
			num++;
			numI->second[2] = num;
		}
		else {
			vector<int> nu;
			nu.push_back(0);
			nu.push_back(0);
			nu.push_back(1);
			matchNum->insert(make_pair(type, nu));
		}

		if (ava_worker.size() == 0)
		{
			return;
		}


		map<int, float> priceSet;
		int mintype;
		float minpri = 0;
		start = clock();
		if (cal_PriceWorker.size() != 0)
		{
			calEveryPlus(&cal_PriceWorker, task->distance, epsilon, &priceSet);
		}
		else {
			calEveryPlus(&ava_worker, task->distance, epsilon, &priceSet);
		}
		end = clock();
		if (priceSet.size() == 0)
		{
			priceFailed++;
			return;
		}

		map<int, float>::iterator ip = priceSet.begin();
		vector<int> avaPlatform;
		for (; ip != priceSet.end(); ip++)
		{
			if (ip->second <= task->value)
			{
				avaPlatform.push_back(ip->first);
			}
		}
		if (avaPlatform.size() == 0)//没有平台价格小于订单价格
		{
			rejectByPrice++;
			return;
		}

		string probfile = "eachprob";
		float minPrice = 1000000;
		int selectedPlatform = -1;
		for (map<int, float>::iterator im = priceSet.begin(); im!=priceSet.end(); im++)
		{
			if (minPrice > im->second) {
				minPrice = im->second;
				selectedPlatform = im->first;
			}
		}
		multiset<Worker> selectedWorker = ava_worker.find(selectedPlatform)->second;
		float outerPrice = priceSet.find(selectedPlatform)->second;

		Worker other_aviWorker;
		start = clock();
		//other_aviWorker = find_other_worker(selectedWorker, outerPrice, task);
		other_aviWorker = find_other_worker(selectedWorker, outerPrice, &pri_task);
		float dis = get_distance(task->x, task->y, other_aviWorker.x, other_aviWorker.y);

		if (other_aviWorker.id == -1)
		{
			rejectFinal++;
		}
		int outworkertype;
		if (other_aviWorker.id != -1 & dis <= other_aviWorker.radius)
		{
			map<int, vector<int>>::iterator findOuter = finishedOuterWorkerTasks.find(other_aviWorker.id);
			if (findOuter != finishedOuterWorkerTasks.end())
			{

				findOuter->second.push_back(task->id);
			}
			else {
				vector<int> temp;
				finishedOuterWorkerTasks.insert(make_pair(other_aviWorker.id, temp));
			}
			outworkertype = other_aviWorker.type;
			map<int, multiset<Worker>*>::iterator imi = otherWorkers->find(other_aviWorker.type);
			for (multiset<Worker>::iterator it = imi->second->begin(); it != imi->second->end(); it++)
			{
				if (it->id == other_aviWorker.id)
				{
					Worker sltWorker = other_aviWorker;
					sltWorker.arrive = task->endtime+1;
					sltWorker.history_value.push_back(outerPrice);
					sltWorker.x = task->end_x;
					sltWorker.y = task->end_y;
					M->push_back(make_tuple(sltWorker.id, task->id, task->value - outerPrice));
					imi->second->erase(it);
					imi->second->insert(sltWorker);
					break;
				}
			}


			map<int, vector<int>>::iterator numII = matchNum->find(type);
			if (numII != matchNum->end())
			{
				int num = numII->second[0];
				int othernum = numII->second[1];
				num++;
				othernum++;
				numII->second[0] = num;
				numII->second[1] = othernum;
			}
			else {
				vector<int> nu;
				nu.push_back(1);
				nu.push_back(1);
				nu.push_back(1);
				matchNum->insert(make_pair(type, nu));
			}

			map<int, vector<double>>::iterator revI = rev->find(type);
			if (revI != rev->end())
			{
				float revn = revI->second[0];
				float otherrev = revI->second[1];
				revn = revn + task->value - outerPrice;
				otherrev = otherrev + outerPrice;
				revI->second[0] = revn;
				revI->second[1] = otherrev;

			}
			else
			{
				vector<double> re;
				//re.push_back(totPrice);
				re.push_back(task->value - outerPrice);
				re.push_back(outerPrice);
				re.push_back(0);
				re.push_back(0);
				rev->insert(make_pair(type, re));
			}
			revI = rev->find(outworkertype);
			if (revI != rev->end())
			{
				float revn = revI->second[0];
				revn = revn + outerPrice;
				revI->second[0] = revn;
				float outnum = revI->second[2];
				revI->second[2] += 1;
				float outrev = revI->second[3];
				revI->second[3] += outerPrice;

			}
			else
			{
				vector<double> re;
				//re.push_back(totPrice);
				re.push_back(outerPrice);
				re.push_back(0);
				re.push_back(1);
				re.push_back(outerPrice);
				rev->insert(make_pair(type, re));
			}
		}
		end = clock();
	}
}

void get_othertasks(vector<Task>* tasks, map<int,multiset<vector<float>>>* othertasks, int typenum) { 
	cout << "get_othertasks start .... " << endl;
	map<int, multiset<vector<float>>> mytasks;
	for (int i = 0; i < (*tasks).size(); i++) {
		Task task = (*tasks)[i];
		int type = task.type;
		map<int, multiset<vector<float>>>::iterator it = mytasks.find(type);
		if (it != mytasks.end()) {
			vector<float> vi;
			vi.push_back(task.distance);
			vi.push_back(task.value);
			it->second.insert(vi);
		}
		else {
			multiset<vector<float>> ntask;
			vector<float> vi;
			vi.push_back(task.distance);
			vi.push_back(task.value);
			ntask.insert(vi);
			mytasks.insert(make_pair(type, ntask));
		}
	}
	map<int, multiset<vector<float>>>::iterator ite = mytasks.begin();
	map<int, multiset<vector<float>>> mytasks_val;
	for (; ite != mytasks.end(); ite++) {
		multiset<vector<float>>::iterator it = (ite->second).begin();
		int type = ite->first;
		float min_dis = floor((*it)[0]);
		float dis_l = min_dis;
		float dis_h = dis_l + 1;
		multiset<float> history_value;
		for (; it != (ite->second).end(); ) {
			float distance = (*it)[0];
			if (distance == 0) {
				it++;
				continue;
			}
			float value = (*it)[1] - 9;
			if (distance < 1) {
				value = 0;
			}
			float avg_value = value / distance;
			if (distance >= dis_l && distance < dis_h) {
				history_value.insert(avg_value);
				it++;
			}
			else {
				multiset<float>::iterator it2 = history_value.begin();
				float min_val = (*it2), max_val = (*it2);
				float valscale = *(--history_value.end()) - *(history_value.begin());
				float perval = valscale / 10;
				vector<float> val_nums, val_num;
				float val = min_val + perval;
				if (perval == 0) {
					val += 0.1;
				}
				int num = 0;
				for (; it2 != history_value.end(); it2++) {
					if ((*it2) < val) {
						num++;
					}
					else {
						val_nums.push_back(float(num) / history_value.size());
						num = 1;
						val += perval;
						while ((*it2) >= val) {
							val_nums.push_back(0);
							val += perval;
						}
					}
				}
				if (num != 0) {
					val_nums.push_back(float(num) / history_value.size());
					num = 0;
				}
				map<int, multiset<vector<float>>>::iterator it3 = mytasks_val.find(type);
				vector<float> vi;
				vi.push_back(type);
				vi.push_back(dis_l);
				vi.push_back(min_val);
				vi.push_back(perval);
				for (int i = 0; i < val_nums.size(); i++) {
					vi.push_back(val_nums[i]);
				}
				if (it3 == mytasks_val.end()) {
					multiset<vector<float>> ntask;
					ntask.insert(vi);
					mytasks_val.insert(make_pair(type, ntask));
				}
				else {
					it3->second.insert(vi);
				}
				dis_l = floor(distance);
				dis_h = dis_l + 1;
				history_value.clear();
			}
		}
	}
	for (int i = 1; i <= typenum; i++) {
		multiset<vector<float>> othertask;
		for (int j = 1; j <= typenum; j++) {
			if (i != j) {
				multiset<vector<float>> subtask;
				subtask = mytasks_val.find(j)->second;
				othertask.insert(subtask.begin(), subtask.end());
			}
			
		}
		othertasks->insert(make_pair(i, othertask));
	}
	cout << "get_othertasks end .... " << endl;
}



void printresult1(map<int, vector<int>> matchNum, map<int, vector<float>> rev, string outputfile)
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
	}
	map<int, vector<float>>::iterator matchrev = rev.begin();
	for (; matchrev != rev.end(); matchrev++)
	{
		cout << "Type " << matchrev->first << " Total Rev " << matchrev->second[0] << " Out Rev" << matchrev->second[1] << endl;
		fout << "Type " << matchrev->first << " Total Rev " << matchrev->second[0] << " Out Rev" << matchrev->second[1] << endl;
	}
}

void basedMatch(map<int, multiset<Worker>*>* workers, vector<Task>* tasks, vector<tuple<int, int, float> >* M, map<int, vector<int>>* matchNum, map<int, vector<double>>* rev,float epsilon, float geoepsilon)
{
	int test=0;
	map<int, multiset<Worker>*> otherworker1, otherworker2,otherworker3;
	map<int, multiset<Worker>*>::iterator mulot = workers->begin();
	for (;mulot!=workers->end();mulot++)
	{
		if (mulot->first == 1)
		{
			otherworker2.insert(make_pair(1, mulot->second));
			otherworker3.insert(make_pair(1, mulot->second));
		}
		else if (mulot->first == 2)
		{
			otherworker1.insert(make_pair(2, mulot->second));
			otherworker3.insert(make_pair(2, mulot->second));
		}
		else
		{
			otherworker1.insert(make_pair(3, mulot->second));
			otherworker2.insert(make_pair(3, mulot->second));
		}

	}
	

	for (int i = 0; i < tasks->size(); i++)
	{
		Task temp = (*tasks)[i];	
		cout << "\r" << i << "/" << tasks->size();
		multiset<Worker>* myworker = workers->find(temp.type)->second;
		if (temp.type == 1)
		{
			Tripsolve(&temp, myworker, &otherworker1,M, matchNum, rev, epsilon,geoepsilon);
			//TripNoSelectionsolve(&temp, myworker, &otherworker1, M, matchNum, rev, epsilon, geoepsilon);
		}
		else if (temp.type == 2)
		{
			Tripsolve(&temp, myworker, &otherworker2, M, matchNum, rev, epsilon, geoepsilon);
			//TripNoSelectionsolve(&temp, myworker, &otherworker2, M, matchNum, rev, epsilon, geoepsilon);;
		}
		else
		{
			Tripsolve(&temp, myworker, &otherworker3, M, matchNum, rev, epsilon,geoepsilon);
			//TripNoSelectionsolve(&temp, myworker, &otherworker3, M, matchNum, rev, epsilon, geoepsilon);
		}
		
	}
	string calFile = "calFile";
	ofstream ofc(calFile,ios::app);
	cout  << " rejectByNoAv: " << rejectByNoAv << " RejectFinal:" << rejectFinal <<" rejectBydis: "<< rejectBydis << " probZero" << probZero;
	ofc << " rejectByNoAv: " << rejectByNoAv << " RejectFinal:" << rejectFinal << " rejectBydis" << rejectBydis << " probZero" << probZero << "rejectByProb" << rejectByProb<<endl;
	


}
