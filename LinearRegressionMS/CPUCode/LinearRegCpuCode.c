#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Maxfiles.h"
#include "MaxSLiCInterface.h"

#define noise_percentage 25.0

void generisi_ulaz(float * data_x1, float * data_x2, float * data_y,
		int numValues);
void ispisi_ulaz(float * data_x1, float * data_x2, float * data_y,
		int numValues);
double AWGN_generator();
double vreme(void);

float slope1;
float slope2;
float intercept;

int main() {

	const int dataPoints = 16000; //total number of data points, should be a multiples of 4

	int i; // loop counter
	double startTime = 0.0, cpuDuration = 0.0, dfeDuration = 0.0; // timer values


	//allocate memory for input data and coefficients
	float *x1 = calloc(dataPoints, sizeof(float)); // input data x1
	float *x2 = calloc(dataPoints, sizeof(float)); // input data x2
	float *y = calloc(dataPoints, sizeof(float)); // input data y

	float *a_cpu = calloc(dataPoints, sizeof(float)); // coefficients a (intercept), from CPU
	float *b1_cpu = calloc(dataPoints, sizeof(float)); // coefficients b (slope), from CPU
	float *b2_cpu = calloc(dataPoints, sizeof(float)); // coefficients b (slope), from CPU
	float *a_dfe = calloc(dataPoints, sizeof(float)); // coefficients a (intercept), from DFE
	float *b1_dfe = calloc(dataPoints, sizeof(float)); // coefficients b (slope), from DFE
	float *b2_dfe = calloc(dataPoints, sizeof(float)); // coefficients b (slope), from CPU


	//generating input data
	printf("Generating input data...\n");
	generisi_ulaz(x1, x2, y, dataPoints);

	//CPU computing
	printf("Computing on CPU...\n");

	float cov_X1_X1, cov_X1_X2, cov_X2_X2;
	float cov_Y_X1, cov_Y_X2;

	cov_X1_X1 = cov_X1_X2 = cov_X2_X2 = cov_Y_X1 = cov_Y_X2 = 0;

	float meanY = 0, meanX1 = 0, meanX2 = 0;

	// racunanje srednjih vrednosti je ulaz za kovarijanse

	startTime = vreme();

	for (i = 0; i < dataPoints; i++) {
		meanY += y[i];
		meanX1 += x1[i];
		meanX2 += x2[i];
	}
	meanY /= dataPoints;
	meanX1 /= dataPoints;
	meanX2 /= dataPoints;

	// -- racunanje srednjih vrednosti

	for (i = 0; i < dataPoints; i++) {

		//TODO: deljenje ukloniti
		//TODO: ovako prilagodjen kod je za DFE;
		// optimalniji gold standard u duhu control flow?

		/*A*/
		cov_Y_X1 += (x1[i] - meanX1) * (y[i] - meanY) / dataPoints;
		/*E*/
		cov_Y_X2 += (x2[i] - meanX2) * (y[i] - meanY) / dataPoints;
		/*D*/
		cov_X1_X2 += (x2[i] - meanX2) * (x1[i] - meanX1) / dataPoints;

		/*C*/
		cov_X1_X1 += (x1[i] - meanX1) * (x1[i] - meanX1) / dataPoints;
		/*F*/
		cov_X2_X2 += (x2[i] - meanX2) * (x2[i] - meanX2) / dataPoints;

		if (i == 0) {
			//cannot compute slope for first point, set to 0 instead
			b1_cpu[i] = 0;
			b2_cpu[i] = 0;
		} else {

			float k = (cov_X1_X1 * cov_X2_X2 - cov_X1_X2 * cov_X1_X2);
			b1_cpu[i] = (cov_Y_X1 * cov_X2_X2 - cov_X1_X2 * cov_Y_X2) / k;
			b2_cpu[i] = (cov_X1_X1 * cov_Y_X2 - cov_Y_X1 * cov_X1_X2) / k;
		}

		a_cpu[i] = meanY - b1_cpu[i] * meanX1 - b2_cpu[i] * meanX2;

	}

	cpuDuration = vreme() - startTime;

	//DFE computing
	printf("Computing on DFE...\n");

	startTime = vreme();

	LinearReg(dataPoints, x1, x2, y, a_dfe, b1_dfe, b2_dfe);

	dfeDuration = vreme() - startTime;

	printf("DFE Done\n");

	//ispisi_ulaz(x1, x2, y, dataPoints);

	//ispis izlaznih rezultata
	for (i = 0; i < dataPoints; i++) {
		printf("Intel:   y = %f + x1 * %f + x2 * %f\n", a_cpu[i], b1_cpu[i],
				b2_cpu[i]);
		printf("DFE:     y = %f + x1 * %f + x2 * %f\n", a_cpu[i], b1_cpu[i],
				b2_cpu[i]);
	}

	printf("delta:     a = %f \t b = %f \t c = %f\n", intercept - a_cpu[i - 1],
			slope1 - b1_cpu[i - 1], slope2 - b2_cpu[i - 1]);
	printf("greska:     a = %f% \t b = %f% \t c = %f%\n",
			(intercept - a_cpu[i - 1])/intercept*100,
				(slope1 - b1_cpu[i - 1])/slope1*100,
				(slope2 - b2_cpu[i - 1])/slope2*100);

	//	printf("Intel:   y = %f + x1 * %f + x2 * %f\n", a_cpu[dataPoints - 1],
	//			b1_cpu[dataPoints - 1], b2_cpu[dataPoints - 1]);

	printf("CPU compute time: %f ms\n", cpuDuration);
	printf("DFE compute time: %f ms\n", dfeDuration);

	//dealokacija memorije
	free(x1);
	free(x2);
	free(y);
	free(a_cpu);
	free(b1_cpu);
	free(b2_cpu);
	free(a_dfe);
	free(b1_dfe);
	free(b2_dfe);

	return 0;
}

//merenje vremena
double vreme() {
	struct timeval vreme;
	gettimeofday(&vreme, NULL);
	return vreme.tv_usec / 1000.0; // ms

	//return 1000000 * vreme.tv_usec+ vreme.;
}

// generisanje ulaznih podataka za racunanje regresije
void generisi_ulaz(float * data_x1, float * data_x2, float * data_y,
		int numValues) {
	srand(time(NULL));

	slope1 = 1 - (float) rand() / RAND_MAX * 2;
	slope2 = 1 - (float) rand() / RAND_MAX * 2;
	// nagib: -1 .. 1
	intercept = 10 - (float) rand() / RAND_MAX * 20;
	// tacka preseka: -10 .. 10

	printf(
			"Ulazni podaci odgovaraju linearnoj fji: (%2f) + (%2f) * x1 + (%2f) * x2 bez dodate buke.\n",
			intercept, slope1, slope2);
	int i;
	// access array elements like array1[ii + dataPoints *jj]

	for (i = 0; i < numValues; i++) {

		float noise = AWGN_generator() * noise_percentage / 100; // 0..100%

		data_x1[i] = 10 * ((float) rand() / RAND_MAX); // 0..10
		data_x2[i] = 10 * ((float) rand() / RAND_MAX); // 0..10

		data_y[i] = intercept + slope1 * data_x1[i] + slope2 * data_x2[i];
		data_y[i] = data_y[i] * (1 - 2 * noise);

		// linear function with gaussian noise limited to a fixed percentage
	}

}

void ispisi_ulaz(float * data_x1, float * data_x2, float * data_y,
		int numValues) {
	printf("*** GENERISANE VREDNOSTI ***\n");
	for (int i = 0; i < numValues; i++) {

		printf("[%d] :\tx1 = %f\tx2 = %f\ty = %f\n", i, data_x1[i], data_x2[i],
				data_y[i]);
	}
	printf("*** KRAJ ISPISA ***\n");

}

#define PI 3.1415926536

double AWGN_generator() {

	/* Generates additive white Gaussian Noise samples with zero mean and a standard deviation of 1. */

	double temp1;
	double temp2;
	double result;
	int p;

	p = 1;

	while (p > 0) {
		temp2 = (rand() / ((double) RAND_MAX));
		if (temp2 == 0) {// temp2 is >= (RAND_MAX / 2)
			p = 1;
		} else {// temp2 is < (RAND_MAX / 2)
			p = -1;
		}
	}

	temp1 = cos((2.0 * (double) PI) * rand() / ((double) RAND_MAX));
	result = sqrt(-2.0 * log(temp2)) * temp1;

	return result;

}
// credit: Additive White Gaussian Noise by Dr Cagri Tanriover @ EmbeddedRelated.com

