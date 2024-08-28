#include "stdafx.h"

#define PI                      3.141592654
#define EARTH_RADIUS            6378.137        //地球近似半径

//两点的距离（纬度，经度）
double get_distance(double lat1, double lng1, double lat2, double lng2)
{
	double radLat1 = lat1 * PI / 180.0;   //角度1? = π / 180
	double radLat2 = lat2 * PI / 180.0;   //角度1? = π / 180
	double a = radLat1 - radLat2;//纬度之差
	double b = lng1 * PI / 180.0 - lng2 * PI / 180.0;  //经度之差
	double dst = 2 * asin((sqrt(pow(sin(a / 2), 2) + cos(radLat1) * cos(radLat2) * pow(sin(b / 2), 2))));
	dst = dst * EARTH_RADIUS;
	dst = round(dst * 10000) / 10000;
	return dst;
}

//计算角度
int get_angle(double lat1, double lng1, double lat2, double lng2)
{
	double x = lat1 - lat2;//t d
	double y = lng1 - lng2;//z y
	int angle = -1;
	if (y == 0 && x > 0) angle = 0;
	if (y == 0 && x < 0) angle = 180;
	if (x == 0 && y > 0) angle = 90;
	if (x == 0 && y < 0) angle = 270;
	if (angle == -1)
	{
		double dislat = get_distance(lat1, lng2, lat2, lng2);
		double dislng = get_distance(lat2, lng1, lat2, lng2);
		if (x > 0 && y > 0) angle = atan2(dislng, dislat) / PI * 180;
		if (x < 0 && y > 0) angle = atan2(dislat, dislng) / PI * 180 + 90;
		if (x < 0 && y < 0) angle = atan2(dislng, dislat) / PI * 180 + 180;
		if (x > 0 && y < 0) angle = atan2(dislat, dislng) / PI * 180 + 270;
	}
	return angle;
}



ZtGeographyCoordinateTransform::ZtGeographyCoordinateTransform()
	: meridianLine(-360), projType('g')
{

}

ZtGeographyCoordinateTransform::~ZtGeographyCoordinateTransform()
{

}

bool ZtGeographyCoordinateTransform::XY2BL(double x, double y, double& lat, double& lon)
{
	if (projType == 'u')
	{
		y = y / 0.9996;
	}

	double bf0 = y / ellipPmt.a0, bf;
	double threshould = 1.0;
	while (threshould > 0.00000001)
	{
		double y0 = -ellipPmt.a2 * sin(2 * bf0) / 2 + ellipPmt.a4 * sin(4 * bf0) / 4 - ellipPmt.a6 * sin(6 * bf0) / 6;
		bf = (y - y0) / ellipPmt.a0;
		threshould = bf - bf0;
		bf0 = bf;
	}

	double t, j2;
	t = tan(bf);
	j2 = ellipPmt.ep2 * pow(cos(bf), 2);

	double v, n, m;
	v = sqrt(1 - ellipPmt.e2 * sin(bf) * sin(bf));
	n = ellipPmt.a / v;
	m = ellipPmt.a * (1 - ellipPmt.e2) / pow(v, 3);

	x = x - 500000;
	if (projType == 'u')
	{
		x = x / 0.9996;
	}

	double temp0, temp1, temp2;
	temp0 = t * x * x / (2 * m * n);
	temp1 = t * (5 + 3 * t * t + j2 - 9 * j2 * t * t) * pow(x, 4) / (24 * m * pow(n, 3));
	temp2 = t * (61 + 90 * t * t + 45 * pow(t, 4)) * pow(x, 6) / (720 * pow(n, 5) * m);
	lat = (bf - temp0 + temp1 - temp2) * 57.29577951308232;

	temp0 = x / (n * cos(bf));
	temp1 = (1 + 2 * t * t + j2) * pow(x, 3) / (6 * pow(n, 3) * cos(bf));
	temp2 = (5 + 28 * t * t + 6 * j2 + 24 * pow(t, 4) + 8 * t * t * j2) * pow(x, 5) / (120 * pow(n, 5) * cos(bf));
	lon = (temp0 - temp1 + temp2) * 57.29577951308232 + meridianLine;

	return true;
}

bool ZtGeographyCoordinateTransform::BL2XY(double lat, double lon, double& x, double& y)
{
	if (meridianLine < -180)
	{
		meridianLine = int((lon + 1.5) / 3) * 3;
	}

	lat = lat * 0.0174532925199432957692;
	double dL = (lon - meridianLine) * 0.0174532925199432957692;

	double X = ellipPmt.a0 * lat - ellipPmt.a2 * sin(2 * lat) / 2 + ellipPmt.a4 * sin(4 * lat) / 4 - ellipPmt.a6 * sin(6 * lat) / 6;
	double tn = tan(lat);
	double tn2 = tn * tn;
	double tn4 = tn2 * tn2;

	double j2 = (1 / pow(1 - ellipPmt.f, 2) - 1) * pow(cos(lat), 2);
	double n = ellipPmt.a / sqrt(1.0 - ellipPmt.e2 * sin(lat) * sin(lat));

	double temp[6] = { 0 };
	temp[0] = n * sin(lat) * cos(lat) * dL * dL / 2;
	temp[1] = n * sin(lat) * pow(cos(lat), 3) * (5 - tn2 + 9 * j2 + 4 * j2 * j2) * pow(dL, 4) / 24;
	temp[2] = n * sin(lat) * pow(cos(lat), 5) * (61 - 58 * tn2 + tn4) * pow(dL, 6) / 720;
	temp[3] = n * cos(lat) * dL;
	temp[4] = n * pow(cos(lat), 3) * (1 - tn2 + j2) * pow(dL, 3) / 6;
	temp[5] = n * pow(cos(lat), 5) * (5 - 18 * tn2 + tn4 + 14 * j2 - 58 * tn2 * j2) * pow(dL, 5) / 120;

	y = X + temp[0] + temp[1] + temp[2];
	x = temp[3] + temp[4] + temp[5];

	if (projType == 'g')
	{
		x = x + 500000;
	}
	else if (projType == 'u')
	{
		x = x * 0.9996 + 500000;
		y = y * 0.9996;
	}

	return true;
}

bool ZtGeographyCoordinateTransform::XYZ2BLH(double x, double y, double z, double& lat, double& lon, double& ht)
{
	double preB, preN;
	double nowB = 0, nowN = 0;
	double threshould = 1.0;

	preB = atan(z / sqrt(x * x + y * y));
	preN = ellipPmt.a / sqrt(1 - ellipPmt.e2 * sin(preB) * sin(preB));
	while (threshould > 0.0000000001)
	{
		nowN = ellipPmt.a / sqrt(1 - ellipPmt.e2 * sin(preB) * sin(preB));
		nowB = atan((z + preN * ellipPmt.e2 * sin(preB)) / sqrt(x * x + y * y));

		threshould = fabs(nowB - preB);
		preB = nowB;
		preN = nowN;
	}
	ht = sqrt(x * x + y * y) / cos(nowB) - nowN;
	lon = atan2(y, x) * 57.29577951308232;    // 180 / pi
	lat = nowB * 57.29577951308232;

	return true;
}

bool ZtGeographyCoordinateTransform::BLH2XYZ(double lat, double lon, double ht, double& x, double& y, double& z)
{
	double sinB = sin(lat / 57.29577951308232);
	double cosB = cos(lat / 57.29577951308232);
	double sinL = sin(lon / 57.29577951308232);
	double cosL = cos(lon / 57.29577951308232);

	double N = ellipPmt.a / sqrt(1.0 - ellipPmt.e2 * sinB * sinB);
	x = (N + ht) * cosB * cosL;
	y = (N + ht) * cosB * sinL;
	z = (N * ellipPmt.b * ellipPmt.b / (ellipPmt.a * ellipPmt.a) + ht) * sinB;

	return true;
}