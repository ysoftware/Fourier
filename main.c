#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <raylib.h>

typedef signed short s16;

#define CHANNELS 2
#define SAMPLE_RATE 44000

typedef struct {
    double real;
    double imaginary;
} DFTResult;

void dft(s16 *interleaved_samples, int number_of_buckets, int number_of_samples, DFTResult *output) {
    for (int i = 0; i < number_of_buckets; i++) {
        output[i].real = 0;
        output[i].imaginary = 0;
        double frequency = (double)i / (double)number_of_samples;
        for (int n = 0; n < number_of_samples; n++) {
            double value = (double)interleaved_samples[n*CHANNELS+0];
            double angle = -2.0 * PI * n * frequency;
            output[i].real += value * cos(angle);
            output[i].imaginary += value * sin(angle);
        }
    }
}

// Cooley-Tukey
void fft(s16 *interleaved_samples, int number_of_samples, DFTResult *output) {
    assert((number_of_samples & (number_of_samples - 1)) == 0 && "number_of_samples must be a power of 2.");

    // copy samples into output buffer as real values, imaginary = 0
    for (int i = 0; i < number_of_samples; i++) {
        output[i].real = (double)interleaved_samples[i * CHANNELS + 0];
        output[i].imaginary = 0.0;
    }

    // compute log2(number_of_samples) — safe because we asserted power-of-two
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
            DFTResult temp   = output[i];
            output[i]        = output[reversed];
            output[reversed] = temp;
        }
    }

    // butterfly stages: step doubles each time (2, 4, 8, ..., number_of_samples)
    // each stage merges pairs of sub-transforms into larger ones
    for (int step = 2; step <= number_of_samples; step *= 2) {

        // base twiddle factor for this stage: e^(-2*PI*i / step)
        // this is the rotation we apply to the odd element of each butterfly
        double angle             = -2.0 * PI / step;
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
                DFTResult even = output[k + j];
                DFTResult odd  = output[k + j + step / 2];

                // rotate the odd element by the current twiddle factor
                double odd_rotated_real      = odd.real * current_real      - odd.imaginary * current_imaginary;
                double odd_rotated_imaginary = odd.real * current_imaginary + odd.imaginary * current_real;

                // butterfly: even+odd goes to top slot, even-odd to bottom slot
                // this is the core reuse: one multiply, two outputs
                output[k + j].real               = even.real + odd_rotated_real;
                output[k + j].imaginary           = even.imaginary + odd_rotated_imaginary;
                output[k + j + step / 2].real      = even.real - odd_rotated_real;
                output[k + j + step / 2].imaginary = even.imaginary - odd_rotated_imaginary;

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

//
// frequency of each bucket depends on the number of samples analyzed
// for 1 second of audio, each bucket is 1hz, so we would need 20k buckets
//

#define NUMBER_OF_BUCKETS 200 // how many buckets of data we collect from the given samples
#define SAMPLES 16384         // how many samples we analyze

int main(void) {
    InitAudioDevice();
    InitWindow(800, 1200, "Fourier");

    s16 *interleaved_samples = malloc(CHANNELS*SAMPLES*sizeof(*interleaved_samples));
    memset(interleaved_samples, 0, CHANNELS*SAMPLES*sizeof(*interleaved_samples));

    double time = GetTime();

#define ADD_FREQUENCY(f, p, a) { phase = p; for (int i = 0; i < SAMPLES; i++) { phase += 2.0f * PI * f / SAMPLE_RATE; interleaved_samples[i*CHANNELS+0] += (s16)(sinf(phase) * a); }}
    float phase = 0;
    ADD_FREQUENCY( 40, 0, 1000);
    ADD_FREQUENCY( 80, 0, 2000);
    ADD_FREQUENCY(120, 0, 3000);
    ADD_FREQUENCY(160, 0, 4000);
    ADD_FREQUENCY(200, 0, 5000);
    ADD_FREQUENCY(240, 0, 6000);
    ADD_FREQUENCY(280, 0, 1000);
    ADD_FREQUENCY(320, 0, 2000);
    ADD_FREQUENCY(360, 0, 3000);
    ADD_FREQUENCY(400, 0, 4000);
    ADD_FREQUENCY(440, 0, 5000);
    ADD_FREQUENCY(480, 0, 6000);

    printf("sample generation done in %f msec\n", (GetTime() - time)*1000.f);

    int number_of_buckets = SAMPLES;

    time = GetTime();
    DFTResult *results = malloc(number_of_buckets*sizeof(*results));
    dft(interleaved_samples, number_of_buckets, SAMPLES, results);
    printf("dft done in %f usec\n", (GetTime() - time)*1000000.f);

    time = GetTime();
    DFTResult *results2 = malloc(number_of_buckets*sizeof(*results2));
    fft(interleaved_samples, SAMPLES, results2);
    printf("fft done in %f usec\n", (GetTime() - time)*1000000.f);

    double epsilon = 1e-3 * (6000.0 * SAMPLES);
    printf("Number of buckets: %d, epsilon: %f\n", number_of_buckets, epsilon);
    for (int i = 0; i < number_of_buckets; i++) {
        double diff_real      = fabs(results[i].real - results2[i].real);
        double diff_imaginary = fabs(results[i].imaginary - results2[i].imaginary);

        if (!(diff_real < epsilon && diff_imaginary < epsilon)) {
            printf("failed to match [%d]: real: %f, imaginary: %f\n", i, diff_real, diff_imaginary);
            printf("values1 = real: %f, imaginary: %f\n", results[i].real, results[i].imaginary);
            printf("values2 = real: %f, imaginary: %f\n", results2[i].real, results2[i].imaginary);
            exit(1);
        }
    }

    printf("Number of buckets: %d\n", NUMBER_OF_BUCKETS);

    Wave wave = {
        .frameCount = SAMPLES,
        .sampleRate = SAMPLE_RATE,
        .sampleSize = 16,
        .channels = CHANNELS,
        .data = interleaved_samples,
    };
    char text[512];

    Sound sound = LoadSoundFromWave(wave);
    PlaySound(sound);

    int scroll_offset = 0;

    while (!WindowShouldClose()) {
        scroll_offset += GetMouseWheelMoveV().y * 2000;

        BeginDrawing(); {
            ClearBackground(BLACK);
            DrawFPS(700, 10);

            for (int i = 0; i < number_of_buckets/2; i++) {
                DFTResult result = results[i];
                float amplitude = 2.0 * sqrt(result.real * result.real + result.imaginary * result.imaginary) / SAMPLES;

                int y = i * 18 + scroll_offset;

                sprintf(text, "%4d Hz", i*SAMPLE_RATE/SAMPLES);
                DrawText(text, 10, y, 18, WHITE);

                DrawRectangle(80, y, amplitude / 10, 18, RED);
            }
        } EndDrawing();
    }

    CloseAudioDevice();

    return 0;
}
