//
//  echoprint-codegen
//  Copyright 2011 The Echo Nest Corporation. All rights reserved.
//


#include "Fingerprint.h"
#include "Params.h"
#include <string.h>

#ifdef _WIN32
#include "win_funcs.h"
#endif

#define SATURATE(var, val) if ((var) > (val)) (var) = (val);
#define FRAMES_TO_PLOT 100000

unsigned int MurmurHash2 ( const void * key, int len, unsigned int seed ) {
    // MurmurHash2, by Austin Appleby http://sites.google.com/site/murmurhash/
    // m and r are constants set by austin
    const unsigned int m = 0x5bd1e995;
    const int r = 24;
    // Initialize the hash to a 'random' value
    unsigned int h = seed ^ len;
    // Mix 4 bytes at a time into the hash
    const unsigned char * data = (const unsigned char *)key;
    while(len >= 4)    {
        unsigned int k = *(unsigned int *)data;
        k *= m;
        k ^= k >> r;
        k *= m;
        h *= m;
        h ^= k;
        data += 4;
        len -= 4;
    }

    // Handle the last few bytes of the input array
    switch(len)    {
        case 3: h ^= data[2] << 16;
        case 2: h ^= data[1] << 8;
        case 1: h ^= data[0];
                h *= m;
    };

    // Do a few final mixes of the hash to ensure the last few
    // bytes are well-incorporated.
    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;
    return h;
}

Fingerprint::Fingerprint(SubbandAnalysis* pSubbandAnalysis, int offset)
    : _pSubbandAnalysis(pSubbandAnalysis), _Offset(offset) { }

/**
 *
 * @param ttarg T_targ: target spacing between onsets reported in any subband. The default value of T_targ is 345 to
 *              aim for about one onset per subband per second. You can reduce it to increase the density of detected
 *              subband onsets (for finer, but more voluminous, fingerprints).
 * @param out
 * @param onset_counter_for_band
 * @return
 */
uint Fingerprint::adaptiveOnsets(int ttarg, matrix_u&out, uint*&onset_counter_for_band) {
    //  E is a sgram-like matrix of energies.

    // pE: an 8-element array representing a single frame.
    const float *pE;
    int bands, frames, i, j, k;

    // The minimum number of frames that have to exist between two onsets for
    // the first onset to be considered "permanent".
    int deadtime = 128;

    double H[SUBBANDS],taus[SUBBANDS];
    int contact[SUBBANDS], lcontact[SUBBANDS], tsince[SUBBANDS];
    double overfact = 1.1;  /* threshold rel. to actual peak */
    uint onset_counter = 0;

    matrix_f E = _pSubbandAnalysis->getMatrix();

    // Take successive stretches of 8 subband samples and sum their energy under a hann window, then hop by 4 samples (50% window overlap).
    int hop = 4;
    int nsm = 8;
    float ham[8];

    // Calculate the Hann window coefficients for 8 samples.
    for(int i = 0 ; i != nsm ; i++)
        ham[i] = .5 - .5*cos( (2.*M_PI/(nsm-1))*i);

    // Calculate the number of frame sets that will be computed?
    int nc =  floor((float)E.size2()/(float)hop)-(floor((float)nsm/(float)hop)-1);

    // Eb is a (nc x 8) matrix.
    matrix_f Eb = matrix_f(nc, 8);

    // Initialize Eb with 0.
    for(uint r=0;r<Eb.size1();r++) for(uint c=0;c<Eb.size2();c++) Eb(r,c) = 0.0;

    // For each set of 8 frames from the subband analysis
    for(i=0;i<nc;i++) {
        // For each subband in the analysis
        for(j=0;j<SUBBANDS;j++) {
            // For each point in the Hann window
            for(k=0;k<nsm;k++)  Eb(i,j) = Eb(i,j) + ( E(j,(i*hop)+k) * ham[k]);
            Eb(i,j) = sqrtf(Eb(i,j));
        }
    }

    frames = Eb.size1();
    bands = Eb.size2();
    pE = &Eb.data()[0];

    out = matrix_u(SUBBANDS, frames);
    onset_counter_for_band = new uint[SUBBANDS];

    double bn[] = {0.1883, 0.4230, 0.3392}; /* preemph filter */   // new
    int nbn = 3;
    double a1 = 0.98;

    // The last filtered level for each subband.
    double Y0[SUBBANDS];

    for (j = 0; j < bands; ++j) {
        onset_counter_for_band[j] = 0;

        taus[j] = 1.0;

        // H: the current threshold
        H[j] = pE[j];

        // If the value exceeded the threshold.
        contact[j] = 0;
        lcontact[j] = 0;
        tsince[j] = 0;
        Y0[j] = 0;
    }

    std::cout << "Detecting onsets for " << frames << " frames." << std::endl;

//    for (j = 0; j < SUBBANDS; ++j) {
//        std::cout << "Values in band: " << j << std::endl;
//        std::cout.precision(std::numeric_limits<float>::max_digits10);
//
//        for (i = 0; i < 1000 /* frames */; ++i) {
//            std::cout << std::fixed << Eb(i, j) << ",";
//        }
//        std::cout << std::endl;
//    }

    for (i = 0; i < frames; ++i) {
        for (j = 0; j < SUBBANDS; ++j) {

            double xn = 0;

            /* Calculate the filter - FIR part */
            if (i >= 2 * nbn) {
                for (int k = 0; k < nbn; ++k) {
                    xn += bn[k] * (pE[j - SUBBANDS * k] - pE[j - SUBBANDS * (2 * nbn - k)]);
                }
            }

            /* IIR part */
            xn = xn + a1 * Y0[j];

            /* Remember the last filtered level */
            Y0[j] = xn;

            contact[j] = (xn > H[j])? 1 : 0;

            if (i < FRAMES_TO_PLOT) {
                std::cout.precision(std::numeric_limits<float>::max_digits10);
                std::cout << "          pE[frame=" << i << "][band=" << j << "] = " << std::fixed << pE[j] << std::endl;
                std::cout << "          xn[frame=" << i << "][band=" << j << "] = " << std::fixed << xn << std::endl;
                std::cout << "           H[frame=" << i << "][band=" << j << "] = " << std::fixed << H[j] << std::endl;
                std::cout << "     contact[frame=" << i << "][band=" << j << "] = " << std::fixed << contact[j] << std::endl;
            }


            if (contact[j] == 1) {
                /* update with new threshold */
                H[j] = xn * overfact;
            } else {
                /* apply decays */
                H[j] = H[j] * exp(-1.0/(double)taus[j]);
            }

            // std::cout << "i = " << i << " - j = " << j <<  " - contact[j] = " << contact[j] << " - lcontact[j] = " << lcontact[j] << std::endl;

            if (contact[j] == 0 && lcontact[j] == 1) {
                /* detach */

                // If the recently detected onset is too close to the previous one,
                // remove the previous one and keep this new instead.
                if (onset_counter_for_band[j] > 0 && (int)out(j, onset_counter_for_band[j] - 1) > i - deadtime) {
                    // overwrite last-written time
                    --onset_counter_for_band[j];
                    --onset_counter;
                }

                out(j, onset_counter_for_band[j]++) = i;
                ++onset_counter;
                tsince[j] = 0;
            }
            ++tsince[j];

            if (tsince[j] > ttarg) {
                taus[j] = taus[j] - 1;
                if (taus[j] < 1) taus[j] = 1;
            } else {
                taus[j] = taus[j] + 1;
            }

            lcontact[j] = contact[j];
        }

        // Increment pE to the next frame (increment the pointer by the number
        // of bands, i.e., the num of columns)
        pE += bands;
    }

    return onset_counter;
}


// dan is going to beat me if i call this "decimated_time_for_frame" like i want to
uint Fingerprint::quantized_time_for_frame_delta(uint frame_delta) {
    double time_for_frame_delta = (double)frame_delta / ((double)Params::AudioStreamInput::SamplingRate / 32.0);
    return ((int)floor((time_for_frame_delta * 1000.0) / (float)QUANTIZE_DT_S) * QUANTIZE_DT_S) / floor(QUANTIZE_DT_S*1000.0);
}

/**
 * Quantize frames with a 8-sized window.
 * Frames 0-7  => 0
 * Frames 8-15 => 1
 * Frames 8-15 => 1
 * @param frame
 * @return
 */
uint Fingerprint::quantized_time_for_frame_absolute(uint frame) {
    // std::cout << "Computing quantized_time_for_frame_absolute for frame: " << frame << std::endl;
    double time_for_frame = _Offset + (double)frame / ((double)Params::AudioStreamInput::SamplingRate / 32.0);

    // std::cout.precision(std::numeric_limits<float>::max_digits10);
    // std::cout << "time_for_frame = " << std::fixed << time_for_frame << std::endl;
    // std::cout << "not_int = " << std::fixed << ((rint((time_for_frame * 1000.0) /  (float)QUANTIZE_A_S) * QUANTIZE_A_S) / floor(QUANTIZE_A_S*1000.0)) << std::endl;

    return ((int)rint((time_for_frame * 1000.0) /  (float)QUANTIZE_A_S) * QUANTIZE_A_S) / floor(QUANTIZE_A_S*1000.0);
}


void Fingerprint::Compute() {
    uint actual_codes = 0;
#if !defined(UNHASHED_CODES)
    unsigned char hash_material[5];
    for(uint i=0;i<5;i++) hash_material[i] = 0;
#endif
    uint * onset_counter_for_band;
    matrix_u out;

    uint onset_count = adaptiveOnsets(345, out, onset_counter_for_band);
    _Codes.resize(onset_count*6);

    std::cout << "Total onsets count: " << onset_count << std::endl;

    for (int b = 0; b < 8; ++b) {
        std::cout << "Onsets on band " << b << ": " << onset_counter_for_band[b] << std::endl;
    }

   for(int band = 0; band < SUBBANDS; band++) {
        std::cout << "Band: " << band << std::endl;
        int line_counter = 0;
        for(int onset = 0; onset < onset_counter_for_band[band]; onset++) {
            while (line_counter < out(band,onset) && line_counter < FRAMES_TO_PLOT) {
                std::cout << "real_contact[frame=" << line_counter << "][band=" << band << "] = " << 0 << std::endl;
                ++line_counter;
            }
            if (out(band,onset) >= FRAMES_TO_PLOT) {
                break;
            }
            std::cout << "real_contact[frame=" << out(band,onset) << "][band=" << band << "] = " << 1 << std::endl;
            ++line_counter;

        }
   }

//    std::cout << "Onsets per band (band, onset, frame, quantized_frame):" << std::endl;
//
//    for(int band = 0; band < SUBBANDS; band++) {
//        std::cout << "Band: " << band << std::endl;
//        for(int onset = 0; onset < onset_counter_for_band[band]; onset++) {
//            std::cout << "(" << band << ", " << onset << ", " << out(band,onset) << ", " << quantized_time_for_frame_absolute(out(band,onset)) << ")" << std::endl;
//        }
//    }


#if defined(UNHASHED_CODES)
    assert(SUBBANDS <= 8);
#endif
    for(unsigned char band=0;band<SUBBANDS;band++) {
        if (onset_counter_for_band[band]>2) {
            for(uint onset=0;onset<onset_counter_for_band[band]-2;onset++) {
                // What time was this onset at?
                uint time_for_onset_ms_quantized = quantized_time_for_frame_absolute(out(band,onset));

                uint p[2][6];
                for (int i = 0; i < 6; i++) {
                    p[0][i] = 0;
                    p[1][i] = 0;
                }
                uint nhashes = 6;

                if ((int)onset == (int)onset_counter_for_band[band]-4)  { nhashes = 3; }
                if ((int)onset == (int)onset_counter_for_band[band]-3)  { nhashes = 1; }

                p[0][0] = (out(band,onset+1) - out(band,onset));
                p[1][0] = (out(band,onset+2) - out(band,onset+1));

                if(nhashes > 1) {
                    p[0][1] = (out(band,onset+1) - out(band,onset));
                    p[1][1] = (out(band,onset+3) - out(band,onset+1));
                    p[0][2] = (out(band,onset+2) - out(band,onset));
                    p[1][2] = (out(band,onset+3) - out(band,onset+2));

                    if(nhashes > 3) {
                        p[0][3] = (out(band,onset+1) - out(band,onset));
                        p[1][3] = (out(band,onset+4) - out(band,onset+1));
                        p[0][4] = (out(band,onset+2) - out(band,onset));
                        p[1][4] = (out(band,onset+4) - out(band,onset+2));
                        p[0][5] = (out(band,onset+3) - out(band,onset));
                        p[1][5] = (out(band,onset+4) - out(band,onset+3));
                    }
                }

                // For each pair emit a code
                for(uint k=0; k < nhashes; k++) {
                    // Quantize the time deltas to 23ms
                    short time_delta0 = (short)quantized_time_for_frame_delta(p[0][k]);
                    short time_delta1 = (short)quantized_time_for_frame_delta(p[1][k]);
                    if (k == 0 && time_delta0 == 0 && time_delta1 == 0) {
                        continue;
                    }
                    uint hashed_code;
#if defined(UNHASHED_CODES)
                    assert(time_delta0 <= 1023);
                    assert(time_delta1 <= 1023);
#if defined(NDEBUG)
                    SATURATE(time_delta0, 1023);
                    SATURATE(time_delta1, 1023);
#endif
                    hashed_code = ((band & 7) << 20) | ((time_delta0 & 1023) << 10) | (time_delta1 & 1023);
#else
                    // Create a key from the time deltas and the band index
                    memcpy(hash_material+0, (const void*)&time_delta0, 2);
                    memcpy(hash_material+2, (const void*)&time_delta1, 2);
                    memcpy(hash_material+4, (const void*)&band, 1);
                    hashed_code = MurmurHash2(&hash_material, 5, HASH_SEED) & HASH_BITMASK;
#endif
                    // Set the code alongside the time of onset
                    _Codes[actual_codes++] = FPCode(time_for_onset_ms_quantized, hashed_code);
                    //fprintf(stderr, "whee %d,%d: [%d, %d] (%d, %d), %d = %u at %d\n", actual_codes, k, time_delta0, time_delta1, p[0][k], p[1][k], band, hashed_code, time_for_onset_ms_quantized);
                }
            }
        }
    }

    _Codes.resize(actual_codes);
    delete [] onset_counter_for_band;
}


