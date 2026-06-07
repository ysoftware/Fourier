### STB-style single header fast fourier transform C library primarily for visualizing frequencies

```c
#define FFT_IMPLEMENTATION
#include "fft.h"

int main(void) {
    signed short *samples; // allocate and fill up

    int number_of_samples = 4096; // must be power of 2
    FFTValue *results = malloc(number_of_samples*sizeof(*results));

    fft_process(results, number_of_samples, 1, samples);

    for (int i = 0; i < number_of_samples; i++) {
        float value = fft_get_bar_value(results, number_of_samples, i, number_of_samples, 1);
        // DrawFrequencyBar(i, value);
    }
}
```

## Please refer to main.c for a complete example using Raylib

```bash
$ make
$ ./main Naraka.mp3
```
