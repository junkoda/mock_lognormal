#ifndef GRID_H
#define GRID_H 1

#include <fftw3.h>

enum FFTMode {fft_mode_unknown, fft_mode_x, fft_mode_k};

class Grid {
 public:
  Grid(const int nc_, const double boxsize_);
  ~Grid();

  void fft_forward();
  void fft_inverse();
  void clear();
  void print();

  double* fx;
  fftw_complex* fk;
  const int nc;
  const double boxsize;
  FFTMode mode;


 private:
  const int ncz;
  fftw_plan plan_forward;
  fftw_plan plan_inverse;

  void print_k();
  void print_x();
};

#endif

