#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <raylib.h>

typedef signed short s16;

#define CHANNELS 2
#define SAMPLE_RATE 44100

typedef struct {
    double real;
    double imaginary;
} DFTResult;

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

int main(int argc, char **argv) {
    InitAudioDevice();
    InitWindow(800, 1200, "Fourier");
    SetWindowState(FLAG_VSYNC_HINT | FLAG_WINDOW_RESIZABLE);

    if (argc == 1) {
        printf("Specify path to the file to visualize.\n");
        return 1;
    }

    Wave wave = LoadWave(argv[1]);
    if (wave.data == NULL) {
        printf("Unable to load wave from file: '%s'.\n", argv[1]);
        return 1;
    }

    WaveFormat(&wave, SAMPLE_RATE, 16, CHANNELS);
    Sound sound = LoadSoundFromWave(wave);
    s16 *samples = (s16*)wave.data;

    int number_of_samples = 8192; // also number of buckets
    DFTResult *results = malloc(number_of_samples*sizeof(*results));

    double start_time = GetTime();
    PlaySound(sound);

    char text[512];
    while (!WindowShouldClose()) {
        double elapsed = GetTime() - start_time;
        if (elapsed >= 0) {
            int current_sample = (int)(elapsed * SAMPLE_RATE);
            int offset = current_sample * CHANNELS;
            fft(samples + offset, number_of_samples, results);
        }

        BeginDrawing(); {
            ClearBackground(BLACK);
            DrawFPS(700, 10);

            sprintf(text, "%f seconds", elapsed);
            DrawText(text, 10, 10, 18, WHITE);

            // second half of buckets is mirror image (not to be used)
            for (int i = 0; i < number_of_samples/2; i++) {
                DFTResult result = results[i];
                float amplitude = 2.0 * sqrt(result.real * result.real + result.imaginary * result.imaginary) / number_of_samples;

                float value = amplitude / 10;
                DrawRectangle(i * 1, 1000 - value, 1, value, RED);
            }
        } EndDrawing();
    }

    CloseAudioDevice();

    return 0;
}
