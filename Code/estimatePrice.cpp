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


void calAvgMin(map<int, float> *avgForEP, int *mintype, float *minpri)
{
	if (avgForEP->size() == 0) {
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
