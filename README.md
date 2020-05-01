# ECSE444_Project

### Algorithms
1. Naive Implementation: DFT
2. Optimization #1: Cooley-Tukey out-of-place FFT
3. Optimization #2: Cooley-Tukey in-place FFT

### Environment
* Linux

### Running the Performance Test
_Note: Feel free to change any of these parameters in [the code](./main_run_multiple.c)._
* number of iterations = 50 (default)
* number of samples to test = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384} (default)
* output csv filename to store performance results = `"./Performance/performance.csv"` (default)

```
gcc -O0 main_run_multiple.c -o performance_test.out -std=c99 -lm

./performance_test.out
```

### Running a Single Iteration to Save the Signals to csv

_Note: Feel free to change any of these parameters in [the code](./main_run_single.c)._
* number of iterations = 50 (default)
* number of samples to test = 256 (default)
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

