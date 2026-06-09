#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <raylib.h>

#define FFT_IMPLEMENTATION
#include "fft.h"

#define CHANNELS 2
#define SAMPLE_RATE 44100

#define WINDOW_SIZE 800

int main(int argc, char **argv) {
    InitAudioDevice();
    InitWindow(WINDOW_SIZE, WINDOW_SIZE, "Fourier");
    SetWindowState(FLAG_VSYNC_HINT);

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
    signed short *samples = (signed short*)wave.data;

    double start_time = GetTime();
    PlaySound(sound);

    // fft settings

    int number_of_samples = 2048;  // also number of buckets
    float frequencies_scale = 1.8; // increase resolution of lower frequencies
    float window_alpha = 0.5;      // smooth out samples at edges of each batch to reduce noise

    FFTValue *results = malloc(number_of_samples*sizeof(*results));

    char text[512];
    while (!WindowShouldClose()) {
        double elapsed = GetTime() - start_time;
        int current_sample = (int)(elapsed * SAMPLE_RATE);
        int offset = current_sample * CHANNELS;

        if (current_sample + number_of_samples > wave.frameCount) break;

        fft_process(results, number_of_samples, CHANNELS, &samples[offset], window_alpha);

        BeginDrawing(); {
            ClearBackground(BLACK);

            int number_of_bars = 30;
            int padding = 10;
            int spacing = 5;

            int bar_width = (WINDOW_SIZE - spacing*(number_of_bars-1) - padding*2) / number_of_bars;
            int max_bar_height = WINDOW_SIZE - padding*2;

            for (int i = 0; i < number_of_bars; i++) {
                float value = fft_get_bar_value(results, number_of_samples, i, number_of_bars, frequencies_scale);
                value = fmin(max_bar_height, fmax(1, value));

                Rectangle frame = {
                    .x = padding + bar_width*i + fmax(0, spacing*(i-1)),
                    .y = padding + max_bar_height/2 - value/2,
                    .width = bar_width,
                    .height = value + 2,
                };
                DrawRectangleRec(frame, RED);
            }
        } EndDrawing();
    }

    CloseAudioDevice();
    CloseWindow();

    return 0;
}
