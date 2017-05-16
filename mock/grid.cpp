#include "grid.h"

FFTmesh::FFTmesh(const int nc_) :
  nc(nc_), ncz(2*(nc_/2+1))
{
  const size_t ngrid= (size_t) nc*nc*ncz;
  fx= (double*) fftwf_malloc(sizeof(double)*ngrid);
  assert(fx);
  
  plan= fftwf_plan_dft_r2c_3d(nc, nc, nc, fx, (fftwf_complex*) fx, 
			      FFTW_ESTIMATE);
    
}

void FFTmesh::fft()
{
  fftwf_execute(plan);
}

FFTmesh::~FFTmesh()
{
  fftwf_free(fx);
  fftwf_destroy_plan(plan);
}

void FFTmesh::clear()
{
  memset(fx, 0, sizeof(float)*(nc*nc*ncz));
}
