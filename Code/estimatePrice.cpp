#include "stdafx.h"






void caltotMin(map<int, float> *minAvg, float distance, int *mintype, float *minpri)
{
	map<int, float>::iterator ot = minAvg->begin();
	*minpri = ot->second;
	*mintype = ot->first;
	for (; ot != minAvg->end(); ot++)
	{
		if (*minpri > ot->second)
		{
			*minpri = ot->second;
			*mintype = ot->first;
		}
	}
	*minpri = *minpri * distance;
}




/*
void (map<int, multiset<Worker>> *avaWorker, map<int, float> *avgForEP)
{
	map<int, multiset<Worker>>::iterator it = avaWorker->begin();
	for (; it != avaWorker->end(); it++)
	{
		float mean=0;
		calMean(&it->second, &mean);
		avgForEP->insert(make_pair(it->first, mean));
	}
}
*/





void calAvgMin(map<int, float> *avgForEP, int *mintype, float *minpri)
{
	if (avgForEP->size() == 0) {//外面公司没有符合条件的报价
		*mintype = -1;
		return;
	}
	map<int, float>::iterator it = avgForEP->begin();
	*mintype = it->first;
	*minpri = it->second;
	for (;it!=avgForEP->end();it++)
	{
		if (it->second < *minpri)
		{
			*mintype = it->first;
			*minpri = it->second;
		}
	}
}