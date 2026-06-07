#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <raylib.h>

typedef signed short s16;
#define PRINT_LINE 20

#define CHANNELS 2
#define NUMBER_OF_BUCKETS (SAMPLES/2 + 1)

#define SAMPLE_RATE 8000

typedef struct {
    double real;
    double imaginary;
} DFTResult;

void dft(s16 *interleaved_samples, int number_of_buckets, int samples, DFTResult *output) {
    for (int i = 0; i < number_of_buckets; i++) {
        for (int n = 0; n < samples; n++) {
            float value = (float)interleaved_samples[n*CHANNELS+0];
            float angle = -2.0 * PI * i * n / samples;
            output[i].real += value * cos(angle);
            output[i].imaginary += value * sin(angle);
        }
    }
}


#define SAMPLES 1000

int main(void) {
    s16 *interleaved_samples = malloc(CHANNELS*SAMPLES*sizeof(*interleaved_samples));
    memset(interleaved_samples, 0, CHANNELS*SAMPLES*sizeof(*interleaved_samples));

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

    DFTResult *results = malloc(NUMBER_OF_BUCKETS*sizeof(*results));
    dft(interleaved_samples, NUMBER_OF_BUCKETS, SAMPLES, results);

    InitWindow(800, 1200, "Fourier");

    char text[512];

    while (!WindowShouldClose()) {
        BeginDrawing(); {
            for (int i = 0; i*SAMPLE_RATE/SAMPLES < 500; i++) {
                DFTResult result = results[i];
                float amplitude = 2.0 * sqrt(result.real * result.real + result.imaginary * result.imaginary) / SAMPLES;

                sprintf(text, "%4d Hz", i*SAMPLE_RATE/SAMPLES);
                DrawText(text, 10, i * 18, 18, WHITE);
                DrawRectangle(80, i * 18, amplitude / 10, 18, RED);
            }
        } EndDrawing();
    }

    return 0;
}
