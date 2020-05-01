# ECSE444_Project

### Algorithms
1. Naive Implementation: DFT
2. Optimization #1: Cooley-Tukey out-of-place FFT
3. Optimization #2: Cooley-Tukey in-place FFT

### Environment
* Linux

### Running the Performance Test
_Note: Code should run without changing the default values. Feel free to change any of these parameters in [the code](./main_run_multiple.c)._
* number of iterations = 50 (default)
* number of samples to test = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384} (default)
* angular frequency of the waveform = 50 (default)
* the waveform = cos ( ang_freq * t_rad) (default)
* output csv filename to store performance results = `"./Performance/performance.csv"` (default)

```
gcc -O0 main_run_multiple.c -o performance_test.out -std=c99 -lm

./performance_test.out
```

### Running a Single Iteration to Save the Signals to csv

_Note: Code should run without changing the default values. Feel free to change any of these parameters in [the code](./main_run_single.c)._
* number of iterations = 1 (default)
* number of samples to test = {256} (default)
* angular frequency of the waveform = 50 (default)
* the waveform = cos ( ang_freq * t_rad) (default)
* output csv filename to store DFT time signal results = `"./DFT/real.csv"` (default)
* output csv filename to store DFT frequency signal results = `"./DFT/complex.csv"` (default)
* output csv filename to store out-of-place FFT time signal results = `"./FFT_out/real.csv"` (default)
* output csv filename to store out-of-place FFT  frequency signal results = `"./FFT_out/complex.csv"` (default)
* output csv filename to store in-place FFT time signal results = `"./FFT_in/real.csv"` (default)
* output csv filename to store in-place FFT frequency signal results = `"./FFT_in/complex.csv"` (default)

```
gcc -O0 main_run_single.c -o run_single_iteration.out -std=c99 -lm

./run_single_iteration.out
```

### Run Jupyter Notebooks to Plot or Test the Algorithm Integrity
* ./plotting.ipynb
*./test_main.ipynb

### csv files structure

* real.csv: only contains real amplitudes
```
y0, y1, y2,...,yN
```
---
* complex.csv: first row are the real values and second row are the imaginary values of complex numbers

If z = a + b * I, where I is the imaginary number, then the first row is: 
```
a0, a1, a2,...,aN
```
and the second row is:
```
b0, b1, b2,...,bN
```
---
* performance.csv: first row is list of nsamples used for testing, second row is DFT's performance (ms), third row is out-of-place FFT's performance (ms), and fourth row is in-place FFT's performance

First row: nsamples
```
64, 128. 256. 512, 1024, 2048, 4096, 8192, 16384
```
Second row: average computation times of DFT per nsamples
```
t_DFT_0, t_DFT_1,..., t_DFT_8
```
etc.
