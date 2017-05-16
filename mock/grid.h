#ifndef GRID_H
#define GRID_H 1

#include <fftw3.h>

class Grid {
 public:
  FFTmesh(const int nc);
  ~FFTmesh();

  void fft();
  int get_nc() const { return nc; }
  float* data(){ return mesh; }
  void clear();

 private:
  const int nc, ncz;
  double* fx;
  fftwf_plan plan;
};

#endif

