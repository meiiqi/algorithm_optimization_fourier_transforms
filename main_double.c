#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#define nsamples 360 // # samples per iteration
#define PI 3.14159265

double complex complexExp(double x)
{
    double arg = 2 * PI * x;
    return cos(arg) + sin(arg) * I;
}

double complex * DFT (double * time_signal)
{
    int s = -1;

    double complex * freq_signal = calloc(nsamples, sizeof(double complex));

    for (int i = 0; i < nsamples; i++) // change to input_size
    {
        for (int j = 0; j < nsamples; j++)
        {
            freq_signal[i] = freq_signal[i] + time_signal[j] * complexExp( (s * i * j) / (double) nsamples);
        }
    }

    return freq_signal;
}

double complex * iDFT (double * time_signal)
{
    int s = 1;
    double complex * freq_signal = malloc(sizeof(double complex) * nsamples);

    for (int i = 0; i < nsamples; i++)
    {
        for (int j = 0; j < nsamples; j++)
        {
            freq_signal[i] = freq_signal[i] + time_signal[j] / nsamples * complexExp( (s * i * j) / (double) nsamples);
        }
    }

    return freq_signal;
}

double complex * FFT (double * time_signal)
{
    return NULL;
}

void generate_signal (double * signal)
{
    int nperiod = 2;
    double sampling_time_rad = nperiod * 2 * PI;

    double t_resolution_rad = sampling_time_rad / nsamples;

    for (int i = 0; i < nsamples; i++)
    {
        double t_rad = i * t_resolution_rad;

        signal[i] = cos(50*t_rad);
        // printf("%d: %f\n", i, signal[i]);
    }
}

void print_complex_to_csv(double complex * array, int array_size)
{ 
    FILE *fp;

    fp = fopen("./complex.csv", "w");

    fprintf(fp, "%f", creal(array[0]));

    for (int i = 1; i < array_size; i++)
    {
        fprintf(fp,",%f", creal(array[i]));
    }

    fprintf(fp, "\n%f", cimag(array[0]));

    for (int i = 1; i < array_size; i++)
    {
        fprintf(fp,",%f", cimag(array[i]));
    }
        
    fclose(fp);
}

void print_real_to_csv(double * array, int array_size)
{
    FILE *fp;
    fp = fopen("./real.csv", "w");

    fprintf(fp, "%f", array[0]);

    for (int i = 1; i < array_size; i++)
    {
        fprintf(fp,",%f", array[i]);
    }

    fclose(fp);    
}

int main(void)
{
    double time_signal[nsamples];
    generate_signal(time_signal);

    // time
    double complex * new_signal = DFT(time_signal);
    // time

    print_complex_to_csv(new_signal, nsamples);
    print_real_to_csv(time_signal, nsamples); 

    // time
    // FFT(NULL);
    // time

    // compare performance
    return 0;
}
