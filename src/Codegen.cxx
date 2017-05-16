//
//  echoprint-codegen
//  Copyright 2011 The Echo Nest Corporation. All rights reserved.
//

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include "Codegen.h"
#include "AudioBufferInput.h"
#include "Fingerprint.h"
#include "Whitening.h"
#include "SubbandAnalysis.h"
#include "Fingerprint.h"
#include "Common.h"

#include "Base64.h"
#include <zlib.h>

using std::string;
using std::vector;

Codegen::Codegen(const float* pcm, unsigned int numSamples, int start_offset) {
    for (int i = 0; i < 2; ++i) {
        is_code_string_cached[i] = false;
    }

    if (Params::AudioStreamInput::MaxSamples < (uint)numSamples)
        throw std::runtime_error("File was too big\n");

    Whitening *pWhitening = new Whitening(pcm, numSamples);
    pWhitening->Compute();

    AudioBufferInput *pAudio = new AudioBufferInput();
    pAudio->SetBuffer(pWhitening->getWhitenedSamples(), pWhitening->getNumSamples());

    SubbandAnalysis *pSubbandAnalysis = new SubbandAnalysis(pAudio);
    pSubbandAnalysis->Compute();

    // ========================================================================

    std::ofstream subbandSpectrogram;
    matrix_f E = pSubbandAnalysis->getMatrix();
    subbandSpectrogram.open("subband_spectrogram.txt");
    subbandSpectrogram.precision(std::numeric_limits<float>::max_digits10);

    for (int i = 0; i < pSubbandAnalysis->getNumFrames(); ++i) {
        for (int j = 0; j < pSubbandAnalysis->getNumBands(); ++j) {
            subbandSpectrogram << i << "," << j << "," << std::fixed << E(j,i) << "\n";
        }
    }

    subbandSpectrogram.close();

    // ========================================================================
    Fingerprint *pFingerprint = new Fingerprint(pSubbandAnalysis, start_offset);
    pFingerprint->Compute();

#if defined(UNHASHED_CODES)
    _CodeString = createCodeStringJSON(pFingerprint->getCodes());
#else
    _CodeString = createCodeString(pFingerprint->getCodes());
#endif
    _NumCodes = pFingerprint->getCodes().size();

    delete pFingerprint;
    delete pSubbandAnalysis;
    delete pWhitening;
    delete pAudio;
}

string Codegen::createCodeString(vector<FPCode> vCodes) {
    if (vCodes.size() < 3) {
        return "";
    }
    std::ostringstream codestream;
    codestream << std::setfill('0') << std::hex;
    for (uint i = 0; i < vCodes.size(); i++)
        codestream << std::setw(5) << vCodes[i].frame;

    for (uint i = 0; i < vCodes.size(); i++) {
        int hash = vCodes[i].code;
        codestream << std::setw(5) << hash;
    }
    return compress(codestream.str());
}

string Codegen::createCodeStringJSON(vector<FPCode> vCodes) {
    std::ostringstream codestream;
    codestream << "[";
    for (uint i = 0; i < vCodes.size(); i++) {
        int hash = vCodes[i].code;
        // codestream << std::setw(5) << hash;
        codestream << "[" <<  vCodes[i].frame << ", "
            << ((hash >> 20) & 7) << ", "
            << ((hash >> 10) & 1023) << ", "
            << ((hash >>  0) & 1023)
            << "]";
        if (i < vCodes.size()-1) {
            codestream << ", ";
        }
    }
    codestream << "]";
    return codestream.str();
}

string Codegen::compress(const string& s) {
    long max_compressed_length = s.size()*2;
    unsigned char *compressed = new unsigned char[max_compressed_length];

    // zlib the code string
    z_stream stream;
    stream.next_in = (Bytef*)(unsigned char*)s.c_str();
    stream.avail_in = (uInt)s.size();
    stream.zalloc = (alloc_func)0;
    stream.zfree = (free_func)0;
    stream.opaque = (voidpf)0;
    deflateInit(&stream, Z_DEFAULT_COMPRESSION);
    do {
        stream.next_out = compressed;
        stream.avail_out = max_compressed_length;
        if(deflate(&stream, Z_FINISH) == Z_STREAM_END) break;
    } while (stream.avail_out == 0);
    uint compressed_length = stream.total_out;
    deflateEnd(&stream);

    // base64 the zlib'd code string
    string encoded = base64_encode(compressed, compressed_length, true);
    delete [] compressed;
    return encoded;
}

std::string Codegen::getCodeString(bool human_readable) {
    const uint n = human_readable;
    if (!is_code_string_cached[n]) {
        is_code_string_cached[n] = true;
        if (human_readable) {
            if (_CodeString.size() > 0) {
                code_string_cache[n] = _CodeString;
            } else {
                code_string_cache[n] = "[]";
            }
        } else {
            if (_CodeString.size() > 0) {
                code_string_cache[n] = '"' + compress(_CodeString) + '"';
            } else {
                code_string_cache[n] = "\"\"";
            }
        }
    }
    return code_string_cache[n];
}
