#ifndef NOISE_H
#define NOISE_H

namespace CurlNoise
{

void curlNoiseSum( 
		float *pos,  	// IN: position (4D)
		float freq,   // IN: frequency
		float freq_tm,
		float alpha,
		float beta,
		float seed,
		float octaves,     // IN: ocatves 
		int oct_leadin,
		float *ampl,
		float *ampl_df,  // divergence free amplitude factor
		float ampl_hf_gain, // high frequency boost
		float ampl_hf_min_oct,
		float *out	// OUT: result vector (3D)
    	);
};

#endif /* NOISE_H */

