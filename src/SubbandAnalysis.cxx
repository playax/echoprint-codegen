//
//  echoprint-codegen
//  Copyright 2011 The Echo Nest Corporation. All rights reserved.
//


#include "SubbandAnalysis.h"
#include "AudioStreamInput.h"

#ifdef _WIN32
#include "win_funcs.h"
#endif

SubbandAnalysis::SubbandAnalysis(AudioStreamInput* pAudio) {
    _pSamples = pAudio->getSamples();
    _NumSamples = pAudio->getNumSamples();
    Init();
}

SubbandAnalysis::SubbandAnalysis(const float* pSamples, uint numSamples) :
    _pSamples(pSamples), _NumSamples(numSamples) {
    Init();
}

SubbandAnalysis::~SubbandAnalysis() {
}

void SubbandAnalysis::Init() {
    // Calculate the analysis filter bank coefficients
    _Mr = matrix_f(M_ROWS, M_COLS);
    _Mi = matrix_f(M_ROWS, M_COLS);
    for (uint i = 0; i < M_ROWS; ++i) {
        for (uint k = 0; k < M_COLS; ++k) {
            _Mr(i,k) = cos((2*i + 1)*(k-4)*(M_PI/16.0));
            _Mi(i,k) = sin((2*i + 1)*(k-4)*(M_PI/16.0));
        }
    }
}

void SubbandAnalysis::Compute() {
    uint frameCounter, i, j;

    // Z buffer: contents of the X buffer (128 samples) * C window function.
    float Z[C_LEN];

    // Sum of 16 samples.
    float Y[M_COLS];

    _NumFrames = (_NumSamples - C_LEN + 1)/SUBBANDS;
    assert(_NumFrames > 0);

    _Data = matrix_f(SUBBANDS, _NumFrames);

    std::cout << "_NumSamples: " << _NumSamples << std::endl; // Exemplo: 2988012
    std::cout << "_NumFrames (for subband analysis): " << _NumFrames << std::endl; // Exemplo: 373485

    // Hop size: 8
    for (frameCounter = 0; frameCounter < _NumFrames; ++frameCounter) {
        // std::cout << "------------ In frameCounter: " << frameCounter << std::endl;
        // Multiply the contents of the 128-sample buffer by the C window function (tabulated in the standard).
        for (i = 0; i < C_LEN; ++i) {
            // std::cout << "Reading sample #" << (frameCounter * SUBBANDS + i) << std::endl;
            Z[i] = _pSamples[frameCounter * SUBBANDS + i] * SubbandFilterBank::C[i];
        }

        // Next: the contents of the Z buffer are divided into 16 8-element
        // vectors which are summed to form a Y vector.

        // Copy the first 16 samples from Z to Y as-is.
        for (i = 0; i < M_COLS; ++i) {
            Y[i] = Z[i];
        }

        // Continue with the copying and summing (128 samples in this frame divided by 16)
        // Result: Y[0]  = Z[0]  + Z[16] + Z[32] + ... + Z[112] (8 samples)
        //         Y[1]  = Z[1]  + Z[17] + Z[33] + ... + Z[113] (8 samples)
        //         ...
        //         Y[15] = Z[15] + Z[31] + Z[47] + ... + Z[127] (8 samples)
        for (i = 0; i < M_COLS; ++i) {
            for (j = 1; j < M_ROWS; ++j) {
                Y[i] += Z[i + M_COLS*j];
            }
        }

        // The Y vector is transformed using a variant of the Modified Discrete
        // Cosine Transform (MDCT) to yield the desired 8 subband values.
        for (i = 0; i < M_ROWS; ++i) {
            float Dr = 0, Di = 0;
            for (j = 0; j < M_COLS; ++j) {
                Dr += _Mr(i,j) * Y[j];
                Di -= _Mi(i,j) * Y[j];
            }
            _Data(i, frameCounter) = Dr*Dr + Di*Di;
        }
    }
}

