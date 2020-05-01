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

void fft_recursion_inplace(double complex * input_output, double complex * copy, int N, int step)
{
	if (step < N) {
		fft_recursion_inplace(copy, input_output, N, step * 2);
		fft_recursion_inplace(copy + step, input_output + step, N, step * 2);
 
		for (int i = 0; i < N; i += 2 * step) {
			double complex t = copy[i + step] * cexp(-I * PI * i / N);
			input_output[i / 2] = copy[i] + t;
			input_output[(i + N)/2] = copy[i] - t;
		}
	}
}

void FFT_inplace(double complex * input_output, int N)
{
    double complex tmp[N];
	for (int i = 0; i < N; i++)
    {
        tmp[i] = input_output[i];
    }
    
	fft_recursion_inplace(input_output, tmp, N, 1);
}

double complex * FFT_outofplace(double complex * input, int N)
{
    int s = -1;

    // create output
    double complex* output = calloc(N, sizeof(double complex));
    if (N == 1)
    {
        output[0] = input[0];
        return output;
    }

    // create f_e_bar and f_o_bar
    double complex* input_even = calloc(N / 2, sizeof(double complex));
    double complex* input_odd = calloc(N / 2, sizeof(double complex));
    for (int i = 0; i < N / 2; i++)
    {
        input_even[i] = input[2 * i + 0];
        input_odd[i] = input[2 * i + 1];
    }

    double complex* f_e_bar = FFT_outofplace(input_even, N/2);
    double complex* f_o_bar = FFT_outofplace(input_odd, N/2);

    // create w
    double complex* w = calloc(N, sizeof(double complex));
    for (int i = 0; i < N; i++)
    {
      w[i] = cexp((2*I* PI*(double)s * (double)i)/ (double)N);
    }

    // concatenate into output
    for (int i = 0; i < N/2; i++)
    {
      output[   i   ] = f_e_bar[i] + w[   i   ] * f_o_bar[i];
      output[i + N/2] = f_e_bar[i] + w[i + N/2] * f_o_bar[i];
    }

    // free memory
    free(w);
    free(input_even);
    free(input_odd);
    free(f_e_bar);
    free(f_o_bar);
    return output;
}

void generate_signal (float ang_freq, double complex * signal, int N)
{
    int nperiod = 2;
    double sampling_time_rad = nperiod * 2 * PI;

    double t_resolution_rad = sampling_time_rad / N;

    for (int i = 0; i < N; i++)
    {
        double t_rad = i * t_resolution_rad;

        signal[i] = cos(ang_freq * t_rad);
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

void print_performance_to_csv(char * filename, int * row_1, double * row_2, double * row_3, double * row_4, int N)
{
    FILE *fp;
    fp = fopen(filename, "w");

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

    fprintf(fp, "\n%f", row_4[0]);

    for (int i = 1; i < N; i++)
    {
        fprintf(fp,",%f", row_4[i]);
    }

    fclose(fp);    
}

int main(void)
{
    // Test parameters
    int nsamples_array[] = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384}; // # samples per iteration (need to be powers of two)
    int n_iterations = 50;
    float ang_freq = 50;

    int array_size = sizeof(nsamples_array)/sizeof(nsamples_array[0]);
    double DFT_avg_time_array[array_size];
    double FFT1_avg_time_array[array_size];
    double FFT2_avg_time_array[array_size];

    for (int i = 0; i < array_size; i++)
    {
        int N = nsamples_array[i];
        printf("[Iteration #%d] Number of samples: %d\n", i, N);

        double complex time_signal[N];
        generate_signal(ang_freq, time_signal, N);
        double DFT_avg_time = 0;
        double FFT1_avg_time = 0;
        double FFT2_avg_time = 0;

        /////////////////////////////////////////////////
        // NAIVE APPROACH - DFT
        for (int j = 0; j < n_iterations; j++)
        {
            clock_t begin = clock();
            double complex * new_signal = DFT(time_signal, N);
            clock_t end = clock();
            DFT_avg_time += (double)(end - begin);
            free(new_signal);
        }
        DFT_avg_time_array[i] = DFT_avg_time / n_iterations;

        printf("Average Time DFT (ms): %f\n", DFT_avg_time_array[i]/CLOCKS_PER_SEC*1000);        
        // /////////////////////////////////////////////////


        // /////////////////////////////////////////////////
        // // OPTIMIZED APPROACH 1 - FFT out-of-place

        for (int j = 0; j < n_iterations; j++)
        {
            clock_t begin = clock();
            double complex * new_signal = FFT_outofplace(time_signal, N);
            clock_t end = clock();
            FFT1_avg_time += (double)(end - begin);
            free(new_signal);
        }
        FFT1_avg_time_array[i] = FFT1_avg_time / n_iterations;

        printf("Average Time FFT outofplace (good) (ms): %f\n", FFT1_avg_time_array[i]/CLOCKS_PER_SEC*1000);
        /////////////////////////////////////////////////

        // /////////////////////////////////////////////////
        // // OPTIMIZED APPROACH 2 - FFT in-place

        for (int j = 0; j < n_iterations; j++)
        {
            clock_t begin = clock();
            FFT_inplace(time_signal, N);
            clock_t end = clock();
            FFT2_avg_time += (double)(end - begin);
        }
        FFT2_avg_time_array[i] = FFT2_avg_time / n_iterations;

        printf("Average Time FFT inplace (ms): %f\n", FFT2_avg_time_array[i]/CLOCKS_PER_SEC*1000);
        /////////////////////////////////////////////////
    }

    print_performance_to_csv("./Performance/performance.csv", nsamples_array, DFT_avg_time_array, FFT1_avg_time_array, FFT2_avg_time_array, array_size);

    return 0;
}
