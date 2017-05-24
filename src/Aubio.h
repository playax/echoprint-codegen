#ifndef ECHOPRINT_CODEGEN_AUBIO_H
#define ECHOPRINT_CODEGEN_AUBIO_H

#include "/Users/devstation2/dev/aubio/aubio/src/aubio.h"
#include "/Users/devstation2/dev/aubio/aubio/src/onset/onset.h"
#include "/Users/devstation2/dev/aubio/aubio/src/fvec.h"
#include "../../aubio/aubio/src/fvec.h"

// FFT params
//uint_t buf_size = 512;
//uint_t hop_size = 256;
//uint_t sample_rate = 44100;

uint_t buf_size = 256;
uint_t hop_size = 128;
uint_t sample_rate = 5512;


aubio_onset_t *echo_o;

fvec_t *echo_onset;
fvec_t *echo_ibuf;

smpl_t is_onset;

void echo_aubio_prepare() {
    echo_o = new_aubio_onset("default", buf_size, hop_size, sample_rate);
    echo_onset = new_fvec(1);
    echo_ibuf = new_fvec(hop_size);
}

uint_t echo_aubio_onset_do(float *aubio_buffer) {

    for (int i = 0; i < hop_size; ++i) {
        echo_ibuf->data[i] = aubio_buffer[i];
    }

    aubio_onset_do(echo_o, echo_ibuf, echo_onset);

    is_onset = fvec_get_sample(echo_onset, 0);

    if (is_onset > 0) {
        return aubio_onset_get_last(echo_o);
    }
    else {
        return 0;
    }

}

#endif //ECHOPRINT_CODEGEN_AUBIO_H
