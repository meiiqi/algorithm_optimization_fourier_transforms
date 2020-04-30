#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#define nsamples 360 // # samples per iteration
#define PI 3.14159265

float complex complexExp(float x)
{
    float arg = 2 * PI * x;
    return cos(arg) + sin(arg) * I;
}

float complex * DFT (float * time_signal)
{
    int s = -1;

    float complex * freq_signal = malloc(sizeof(float complex) * nsamples);

    for (int i = 0; i < nsamples; i++) // change to input_size
    {
        for (int j = 0; j < nsamples; j++)
        {
            freq_signal[i] = freq_signal[i] + time_signal[j] * complexExp( (s * i * j) / nsamples);
            // printf("Freq signal = %.3f + %.3fi\n", creal(freq_signal[i]), cimag(freq_signal[i]));
        }
    }

    return freq_signal;
}

float complex * iDFT (float * time_signal)
{
    int s = 1;
    float complex * freq_signal = malloc(sizeof(float complex) * nsamples);

    for (int i = 0; i < nsamples; i++)
    {
        for (int j = 0; j < nsamples; j++)
        {
            freq_signal[i] = freq_signal[i] + time_signal[j] / nsamples * complexExp( (s * i * j) / nsamples);
            // printf("Freq signal = %.3f + %.3fi\n", creal(freq_signal[i]), cimag(freq_signal[i]));
        }
    }

    return freq_signal;
}

float complex * FFT (float * time_signal)
{
    return NULL;
}

void generate_signal (float * signal)
{
    int nperiod = 2;
    float sampling_time_rad = nperiod * 2 * PI;
    float sampling_time_deg = nperiod * 360;

    float t_resolution_rad = sampling_time_rad / nsamples;
    float t_resolution_deg = sampling_time_deg / nsamples;

    for (int i = 0; i < nsamples; i++)
    {
        float t_rad = i * t_resolution_rad;
        float t_deg = i * t_resolution_deg;

        signal[i] = cos(t_rad);
        // signal[i] = 5 + 4 * cos(2 * t_rad) + 2 * cos(8 * t_rad - PI / 2) + 3 * cos(32 * t_rad + PI/2);
        // printf("generating: %f with t = %f\n", sin(t_rad), t_deg);
    }
}

void print_complex_to_csv(char * filename, float complex * array, int array_size)
{ 
    FILE *fp;

    fp = fopen("./complex.csv", "w");

    fprintf(fp, "%f", crealf(array[0]));

    for (int i = 1; i < array_size; i++)
    {
        fprintf(fp,",%f", crealf(array[i]));
    }

    fprintf(fp, "\n%f", cimagf(array[0]));

    for (int i = 1; i < array_size; i++)
    {
        fprintf(fp,",%f", cimagf(array[i]));
    }
        
    fclose(fp);
}

void print_real_to_csv(char * filename, float * array, int array_size)
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
    float time_signal[nsamples];
    generate_signal(time_signal);

    // time
    float complex * new_signal = DFT(time_signal);

    print_complex_to_csv("", new_signal, nsamples);
    print_real_to_csv("", time_signal, nsamples); 

    // time
    // FFT(NULL);
    // time

    // compare performance

    return 0;
}
