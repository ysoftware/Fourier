#ifndef FFT_H_
#define FFT_H_

#include <assert.h>
#include <math.h>
#include <stdbool.h>

typedef struct {
    double real;
    double imaginary;
} FFTValue;

// Using Cooley-Tukey fast Fourier transform algorithm, process given samples and write the data into results buffer.
//
// Buffer must be at least of size number_of_samples*sizeof(*results).
void fft_process(FFTValue *results, int number_of_samples, int channels, signed short *interleaved_samples);

/// Supports logarithmic scaling of the frequencies (represent lower frequencies in mode definition than higher).
///
/// Frequency of each bucket depends on the number of samples analyzed,
/// so for 1 second of audio, each bucket is 1hz, so we would need 20k buckets.
/// This basically makes number_of_samples the same as number_of_buckets.
double fft_get_bar_value(FFTValue *results, int number_of_samples, int bar, int number_of_bars, bool is_logarithmic);

#endif // FFT_H_

#define FFT_IMPLEMENTATION
#ifdef FFT_IMPLEMENTATION

#ifndef FFT_IMPLEMENTATION_
#define FFT_IMPLEMENTATION_

#define FFT_PI 3.14159265358979323846

void fft_process(FFTValue *results, int number_of_samples, int channels, signed short *interleaved_samples) {
    assert((number_of_samples & (number_of_samples - 1)) == 0 && "number_of_samples must be a power of 2.");

    // copy samples into results buffer as real values, imaginary = 0
    for (int i = 0; i < number_of_samples; i++) {
        results[i].real = (double)interleaved_samples[i * channels + 0];
        results[i].imaginary = 0.0;
    }

    // compute log2(number_of_samples)
    int log2n = 0;
    for (int i = number_of_samples; i > 1; i >>= 1) log2n++;

    // bit-reversal permutation: reorder samples so the butterfly stages
    // can work in-place. for example with 8 samples: index 1 (001) swaps
    // with index 4 (100), index 3 (011) swaps with index 6 (110), etc.
    // we only swap when reversed > i to avoid swapping back
    for (int i = 0; i < number_of_samples; i++) {
        int reversed = 0;
        int index = i;
        for (int bit = 0; bit < log2n; bit++) {
            reversed = (reversed << 1) | (index & 1);
            index >>= 1;
        }
        if (reversed > i) {
            FFTValue temp    = results[i];
            results[i]        = results[reversed];
            results[reversed] = temp;
        }
    }

    // butterfly stages: step doubles each time (2, 4, 8, ..., number_of_samples)
    // each stage merges pairs of sub-transforms into larger ones
    for (int step = 2; step <= number_of_samples; step *= 2) {

        // base twiddle factor for this stage: e^(-2*PI*i / step)
        // this is the rotation we apply to the odd element of each butterfly
        double angle             = -2.0 * FFT_PI / step;
        double twiddle_real      = cos(angle);
        double twiddle_imaginary = sin(angle);

        // walk through each sub-transform of size step in the array
        for (int k = 0; k < number_of_samples; k += step) {

            // current is the running twiddle, starts at e^0 = (1, 0)
            // and gets multiplied by the base twiddle each inner iteration
            double current_real      = 1.0;
            double current_imaginary = 0.0;

            // process each butterfly pair within this sub-transform
            for (int j = 0; j < step / 2; j++) {
                FFTValue even = results[k + j];
                FFTValue odd  = results[k + j + step / 2];

                // rotate the odd element by the current twiddle factor
                double odd_rotated_real      = odd.real * current_real      - odd.imaginary * current_imaginary;
                double odd_rotated_imaginary = odd.real * current_imaginary + odd.imaginary * current_real;

                // butterfly: even+odd goes to top slot, even-odd to bottom slot
                // this is the core reuse: one multiply, two resultss
                results[k + j].real               = even.real + odd_rotated_real;
                results[k + j].imaginary           = even.imaginary + odd_rotated_imaginary;
                results[k + j + step / 2].real      = even.real - odd_rotated_real;
                results[k + j + step / 2].imaginary = even.imaginary - odd_rotated_imaginary;

                // advance the running twiddle by multiplying by base twiddle
                // equivalent to e^(-2*PI*i*(j+1)/step) without calling cos/sin
                double next_real      = current_real * twiddle_real      - current_imaginary * twiddle_imaginary;
                double next_imaginary = current_real * twiddle_imaginary + current_imaginary * twiddle_real;

                current_real      = next_real;
                current_imaginary = next_imaginary;
            }
        }
    }
}

/// Supports logarithmic scaling of the frequencies (represent lower frequencies in mode definition than higher).
double fft_get_bar_value(FFTValue *results, int number_of_samples, int bar, int number_of_bars, bool is_logarithmic) {
    // second half of buckets is mirror image (not to be used)
    int valid_buckets = number_of_samples/2;

    assert(number_of_bars < valid_buckets/2);
    assert(bar < number_of_bars);

    int accumulate_bars = valid_buckets / number_of_bars;
    double value = 0;

    int start, end;
    if (is_logarithmic) {
        start = (int)(valid_buckets * pow((double)bar / number_of_bars, 2.0));
        end = (int)(valid_buckets * pow((double)(bar + 1) / number_of_bars, 2.0));
    } else {
        start = accumulate_bars * bar;
        end = start + accumulate_bars;
    }

    for (int i = start; i < end; ++i) {
        FFTValue result = results[i];
        float amplitude = 2.0 * sqrt(result.real * result.real + result.imaginary * result.imaginary) / number_of_samples;
        value += amplitude;
    }

    return value / accumulate_bars;
}

#endif // FFT_IMPLEMENTATION_
#endif // FFT_IMPLEMENTATION
