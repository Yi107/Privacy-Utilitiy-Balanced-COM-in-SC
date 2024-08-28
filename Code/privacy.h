#ifndef PRIVACY_H
#define PRIVACY_H
void set_privacy(float x, float y, float epsilon, float* newx, float* newy);
float noisyCount(float sensitivity, float epsilon);
float noiseSettingForLocal(float epsilon, float radius, float *noise);
int randomAngle(int* random);
void locationPerturb(double real_lat, double real_lon, float epsilon, float radius, double *perturbed_lat, double *perturbed_lon, int *angle); 
void extendRadius(float rad, float epsilon, float *exrad);
void avaWLocationPerturb(double real_lat, double real_lon, int angle, double *perturbed_lat, double *perturbed_lon);
#endif // !PRIVACY_H
