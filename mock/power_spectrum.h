#ifndef POWER_SPECTRUM_H
#define POWER_SPECTRUM_H 1

#include <fftw3.h>

void compute_power_spectrum(const char filename[],
		    fftw_complex* fk,
		    const int nc, const double boxsize,
		    const double k_min, const double k_max, const double dk,
			    const int mas_correction);



#endif
