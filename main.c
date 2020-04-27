#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define sample_rate 180 // # samples per iteration
#define PI 3.14159265

float complex complexExp(float x)
{
    float arg = 2 * PI * x;
    return cos(arg) + sin(arg) * I;
}

float complex * DFT (float * time_signal)
{
    int s = -1;
    float complex * freq_signal = malloc(sizeof(float complex) * sample_rate);

    for (int i = 0; i < sample_rate; i++)
    {
        for (int j = 0; j < sample_rate; j++)
        {
            freq_signal[i] = freq_signal[i] + time_signal[j] * complexExp( (s * i * j) / sample_rate);
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
    float resolution_rad = 2 * PI / sample_rate;
    float resolution_deg = 360 / sample_rate;

    for (int i = 0; i < sample_rate; i++)
    {
        float t_rad = i * resolution_rad;
        float t_deg = i * resolution_deg;

        signal[i] = sin(t_rad);
        // printf("generating: %f with t = %f\n", sin(t_rad), t_deg);
    }
}

int main(void)
{
    // printf("hello!\n");

    float time_signal[sample_rate];
    generate_signal(time_signal);

    // time
    float complex * new_signal = DFT(time_signal);

    // test
    // for (int i = 0; i < sample_rate; i++)
    // {
    //     printf("[%i]: %f + 1 ==? %f \n", i, time_signal[i], new_signal[i]);
    // }

    for (int i = 0; i < sample_rate; i++)
    {
        printf("[%d] %.2f + %.2fi\n", i, creal(new_signal[i]), cimag(new_signal[i]));
    }

    // float complex z1 = 1.0 + 3.0 * I;
    // float complex z2 = 1.0 - 4.0 * I;

    // printf("Starting values: Z1 = %.2f + %.2fi\tZ2 = %.2f %+.2fi\n", 
    //        creal(z1), 
    //        cimag(z1), 
    //        creal(z2), 
    //        cimag(z2));

    // float complex z3 = ccos(45.0 * PI / 180) + csin(45.0 * PI / 180) * I;
    // float complex z3 = cos(45.0 * PI / 180) + sin(45.0 * PI / 180) * I;

    // printf("Z3 = %.3f + %.3fi", creal(z3),cimag(z3));

 

    // time
    // FFT(NULL);
    // time

    // compare performance

    return 0;
}
