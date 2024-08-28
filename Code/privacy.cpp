#define _USE_MATH_DEFINES
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include<iomanip>
#include <math.h>
#include <random>
#include "stdafx.h"

using namespace std;



double fasterlambertw(double x)
{



	double logterm = log(-x);
	double loglogterm = log(-logterm);


	double w = logterm - loglogterm + loglogterm / logterm;
	double expw = exp(-w);

	return (w * w + expw * x) / (1.0f + w);
}

double unidistribute(double a, double b)
{
	double rnd;
	for (int i = 0; i < 10; i++)
	{
		rnd = (double)rand() / RAND_MAX * (b - a) + a;
	}
	rnd = (double)rand() / RAND_MAX * (b - a) + a;
	return rnd;

}

double cdf(float epsilon)
{
	double c;
	double privacy, cdfpara;
	double ramdomp = unidistribute(0, 1);
	privacy = -1 / epsilon;
	cdfpara = (ramdomp - 1) * exp(-1);
	c = privacy * (fasterlambertw(cdfpara) + 1);
	return c;
}

int Sgn(float d)
{
	if (d < 0) return -1;
	else if (d == 0) return 0;
	else return 1;
}

float absfloat(float b)
{
	if (b > 0)
	{
		return b;
	}
	else
	{
		return -b;
	}
}

void set_privacy(float x, float y, float epsilon, float* newx, float* newy)
{

	double theta;
	double p;
	theta = unidistribute(0, M_PI);
	double si, co;
	si = sin(theta);
	co = cos(theta);
	double r;
	r = cdf(epsilon);
	*newx = x + r * co;
	*newy = y + r * si;
}


float noisyCount(float sensitivity, float epsilon)
{
	//assert(epsilon > 0);
	double d = unidistribute(0, 1);

	float uniform = (float)d - 0.5;
	while (abs(uniform - 0.5) < 0.00009)
	{
		d = unidistribute(0, 1);
		uniform = (float)d - 0.5;
	}
	float result = sensitivity / epsilon * Sgn(uniform) * log(1 - 2.0 * absfloat(uniform));
	return result;
}

int randomAngle(int* random)
{
	random_device rd;
	mt19937 eng(rd());
	uniform_int_distribution<> distr(0, 360);
	int random_angle = distr(eng);
	*random = random_angle;
	return 0;
}


float randomNum(float *random)
{
	random_device rd;
	mt19937 eng(rd());
	uniform_real_distribution<float> distr(0.0f, 1.0f);

	float random_num = distr(eng);
	*random = random_num;
	return 0;
}



//[-r,r] 概率分布满足ppt2
float noiseSettingForLocal(float epsilon, float radius, float *noise)
{
	float random_u = 0;
	randomNum(&random_u);
	float beta = 1 / epsilon;
	float multi = 1.0 / (1 - exp(-radius * epsilon));
	float temp = exp(-radius * epsilon);
	if (random_u <= 0.5)
	{
		*noise = beta * log(temp + 2 * random_u / multi);
	}
	else {
		*noise = -beta * log(temp + 2 * (1 - random_u) / multi);
	}

	return *noise;
}


void locationPerturb(double real_lat, double real_lon, float epsilon, float radius, double *perturbed_lat, double *perturbed_lon, int *angle)
{
	
	int random = 0;
	randomAngle(&random);
	*angle = random;
	float cosValue = cos((float)random * M_PI / 180.0);
	float sinValue = sin((float)random * M_PI / 180.0);
	float rrand = 0;
	noiseSettingForLocal(epsilon, radius, &rrand);//noise的单位是km,要在km上做扰动
	ZtGeographyCoordinateTransform ztGCT;
	double lx = 0, ly = 0;
	float newx = 0, newy = 0;
	double lat, lon, ht;
	ztGCT.BL2XY(real_lat, real_lon, lx, ly);
	//lx,ly都是以米为单位的
	//cout << lx << " " << ly << endl;
	lx = lx / 1000;
	ly = ly / 1000;
	//cout << rrand << endl;
	//location privacy alogrithm need to be changed. Input lat,lng, out put: perturbed lat lng, theta
	//set_privacy(lx, ly, geoepsilon, &newx, &newy);
	newx = lx + cosValue *rrand;
	newy = ly + sinValue *rrand;
	newx = newx * 1000;
	newy = newy * 1000;
	//cout << newx << " " << newy << endl;
	
	ztGCT.XY2BL(newx, newy, lat, lon);
	*perturbed_lat = lat;
	*perturbed_lon = lon;
	//cout << lat << " " << lon << endl;

}


void extendRadius(float rad, float epsilon, float *exrad)
{

	float beta = 1 / epsilon;
	float multi = (1 - exp(-1 *rad* epsilon));
	float random_u = 0;
	randomNum(&random_u);
	float noise = -beta * log(1 - random_u * multi);
	*exrad = rad + noise;

}

void avaWLocationPerturb(double real_lat, double real_lon, int angle, double *perturbed_lat, double *perturbed_lon)
{
	float cosValue = cos((float)angle * M_PI / 180.0);
	float sinValue = sin((float)angle * M_PI / 180.0);
	float randz = 0;
	float epsilon = 0.7;
	float beta = 1 / epsilon;
	float random_u = 0;
	randomNum(&random_u);
	if (random_u <= 0.5)
	{
		randz =  beta * log(2 * random_u);
	}
	else {
		randz = -beta * log(2 * (1 - random_u));
	}
	ZtGeographyCoordinateTransform ztGCT;
	double lx = 0, ly = 0;
	float newx = 0, newy = 0;
	double lat, lon, ht;
	ztGCT.BL2XY(real_lat, real_lon, lx, ly);
	//lx,ly都是以米为单位的
	//cout << lx << " " << ly << endl;
	lx = lx / 1000;
	ly = ly / 1000;
	//cout << rrand << endl;
	//location privacy alogrithm need to be changed. Input lat,lng, out put: perturbed lat lng, theta
	//set_privacy(lx, ly, geoepsilon, &newx, &newy);
	newx = lx + cosValue * randz;
	newy = ly + sinValue * randz;
	newx = newx * 1000;
	newy = newy * 1000;
	//cout << newx << " " << newy << endl;
	ztGCT.XY2BL(newx, newy, *perturbed_lat, *perturbed_lon);

}
