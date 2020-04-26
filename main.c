#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define sample_rate 360 // # samples per iteration
#define PI 3.14159265

float * DFT (float * time_signal)
{
    float * freq_signal = malloc(sizeof(float) * sample_rate);
    for (int i = 0; i < sample_rate; i++)
    {
        freq_signal[i] = time_signal[i] + 1;
    }
    return freq_signal;
}

void * FFT (void * time_signal)
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
        printf("generating: %f with t = %f\n", sin(t_rad), t_deg);
    }
}

int main(void)
{
    printf("hello!");

    float time_signal[sample_rate];
    generate_signal(time_signal);

    // time
    float * new_signal = DFT(time_signal);

    // test
    for (int i = 0; i < sample_rate; i++)
    {
        printf("[%i]: %f + 1 ==? %f \n", i, time_signal[i], new_signal[i]);
    }


    // time
    FFT(NULL);
    // time

    // compare performance

    return 0;
}
