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

int main(int argc, char **argv) {
    InitAudioDevice();
    InitWindow(800, 800, "Fourier");
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

    int number_of_samples = 8192; // also number of buckets
    FFTValue *results = malloc(number_of_samples*sizeof(*results));

    double start_time = GetTime();
    PlaySound(sound);

    char text[512];
    while (!WindowShouldClose()) {
        double elapsed = GetTime() - start_time;
        int current_sample = (int)(elapsed * SAMPLE_RATE);
        int offset = current_sample * CHANNELS;

        if (current_sample + number_of_samples > wave.frameCount) break;
        fft_process(results, number_of_samples, CHANNELS, &samples[offset]);

        int number_of_bars = 100;
        int padding = 10;
        int size = 780;

        BeginDrawing(); {
            ClearBackground(BLACK);

            for (int i = 0; i < number_of_bars; i++) {
                float value = fft_get_bar_value(results, number_of_samples, i, number_of_bars, true);
                Rectangle frame = {
                    .x = padding + i * ((float)size / number_of_bars),
                    .y = padding + size - fmin(size, fmax(1, value)),
                    .width = ((float)size / number_of_bars),
                    .height = fmin(size, fmax(1, value))
                };
                DrawRectangleRec(frame, RED);
            }
        } EndDrawing();
    }

    CloseAudioDevice();
    CloseWindow();

    return 0;
}
