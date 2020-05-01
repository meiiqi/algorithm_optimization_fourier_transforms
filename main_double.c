#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define PI 3.14159265

double complex * DFT (double complex * time_signal, int N)
{
    double complex * freq_signal = calloc(N, sizeof(double complex));

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            freq_signal[i] = freq_signal[i] + creal(time_signal[j]) * cexp( (-2 * I * PI * i * j) / N);
        }
    }

    return freq_signal;
}

double complex * iDFT (double * time_signal, int N)
{
    int s = 1;
    double complex * freq_signal = malloc(sizeof(double complex) * N);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            freq_signal[i] = freq_signal[i] + time_signal[j] / N * cexp( (s * i * j) / (double) N);
        }
    }

    return freq_signal;
}

// https://rosettacode.org/wiki/Fast_Fourier_transform#C
void fft_recursion(double complex * copy, double complex * output, int N, int step)
{
	if (step < N) {
		fft_recursion(output, copy, N, step * 2);
		fft_recursion(output + step, copy + step, N, step * 2);
 
		for (int i = 0; i < N; i += 2 * step) {
			double complex t = output[i + step] * cexp(-I * PI * i / N);
			copy[i / 2] = output[i] + t;
			copy[(i + N)/2] = output[i] - t;
		}
	}
}

double complex * FFT(double complex * input, int N)
{
    double complex * output = calloc(N, sizeof(double complex));
    double complex copy[N];
	for (int i = 0; i < N; i++)
    {
        output[i] = input[i];
        copy[i] = input[i];
    }
    
	fft_recursion(copy, output, N, 1);
    return output;
}

void generate_signal (double complex * signal, int N)
{
    int nperiod = 2;
    double sampling_time_rad = nperiod * 2 * PI;

    double t_resolution_rad = sampling_time_rad / N;

    for (int i = 0; i < N; i++)
    {
        double t_rad = i * t_resolution_rad;

        signal[i] = cos(50 * t_rad);
    }
}

void print_complex_to_csv(double complex * array, int N)
{ 
    FILE *fp;

    fp = fopen("./complex.csv", "w");

    // print on first row
    fprintf(fp, "%f", creal(array[0]));

    for (int i = 1; i < N; i++)
    {
        fprintf(fp,",%f", creal(array[i]));
    }

    // print on second row
    fprintf(fp, "\n%f", cimag(array[0]));

    for (int i = 1; i < N; i++)
    {
        fprintf(fp,",%f", cimag(array[i]));
    }
        
    fclose(fp);
}

void print_real_to_csv(double complex * array, int N)
{
    FILE *fp;
    fp = fopen("./real.csv", "w");

    fprintf(fp, "%f", creal(array[0]));

    for (int i = 1; i < N; i++)
    {
        fprintf(fp,",%f", creal(array[i]));
    }

    fclose(fp);    
}

void print_performance_to_csv(int * row_1, double * row_2, double * row_3, int N)
{
    FILE *fp;
    fp = fopen("./performance.csv", "w");

    fprintf(fp, "%d", row_1[0]);

    for (int i = 1; i < N; i++)
    {
        fprintf(fp,",%d", row_1[i]);
    }

    fprintf(fp, "\n%f", row_2[0]);

    for (int i = 1; i < N; i++)
    {
        fprintf(fp,",%f", row_2[i]);
    }

    fprintf(fp, "\n%f", row_3[0]);

    for (int i = 1; i < N; i++)
    {
        fprintf(fp,",%f", row_3[i]);
    }

    fclose(fp);    
}

int main(void)
{
    int nsamples_array[] = {64, 128, 256, 512, 1024, 2048}; // # samples per iteration (need to be powers of two)
    int n_iterations = 100;
    int array_size = sizeof(nsamples_array)/sizeof(nsamples_array[0]);
    double DFT_avg_time_array[array_size];
    double FFT_avg_time_array[array_size];

    for (int i = 0; i < array_size; i++)
    {
        int N = nsamples_array[i];
        printf("%d: Number of samples: %d\n", i, N);

        double complex time_signal[N];
        generate_signal(time_signal, N);
        double DFT_avg_time = 0;
        double FFT_avg_time = 0;

        /////////////////////////////////////////////////
        // NAIVE APPROACH
        for (int j = 0; j < n_iterations; j++)
        {
            clock_t begin = clock();
            double complex * new_signal = DFT(time_signal, N);
            clock_t end = clock();
            DFT_avg_time += (double)(end - begin);
            free(new_signal);
        }
        DFT_avg_time_array[i] = DFT_avg_time / n_iterations;

        printf("Average Time DFT (ms): %f\n", DFT_avg_time_array[i]);
        
        // // --- Exporting the signals to CSV for plotting in Python --- 
        // // print_real_to_csv(time_signal, N); 
        // // print_complex_to_csv(new_signal, N);
        // /////////////////////////////////////////////////

        // /////////////////////////////////////////////////
        // // OPTIMIZED APPROACH

        for (int j = 0; j < n_iterations; j++)
        {
            clock_t begin = clock();
            double complex * new_signal = FFT(time_signal, N);
            clock_t end = clock();
            FFT_avg_time += (double)(end - begin);
            free(new_signal);
        }
        FFT_avg_time_array[i] = FFT_avg_time / n_iterations;

        printf("Average Time FFT (ms): %f\n", FFT_avg_time_array[i]);

        // --- Exporting the signals to CSV for plotting in Python --- 
        // print_real_to_csv(time_signal, N); 
        // print_complex_to_csv(new_signal, N);
        /////////////////////////////////////////////////
    }

    print_performance_to_csv(nsamples_array, DFT_avg_time_array, FFT_avg_time_array, array_size);

    return 0;
}
