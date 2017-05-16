#include <iostream>
#include <cassert>
#include <cstring>
#include <cmath>
#include "grid.h"

using namespace std;

Grid::Grid(const int nc_, const double boxsize_) :
  nc(nc_), ncz(2*(nc_/2+1)), boxsize(boxsize_)
{
  const size_t ngrid= (size_t) nc*nc*ncz;
  fx= (double*) fftw_malloc(sizeof(double)*ngrid);
  assert(fx);
  fk= reinterpret_cast<fftw_complex*>(fx);
  
  plan_forward= fftw_plan_dft_r2c_3d(nc, nc, nc, fx, fk,
				     FFTW_ESTIMATE);

  plan_inverse= fftw_plan_dft_c2r_3d(nc, nc, nc, fk, fx,
				     FFTW_ESTIMATE);

  cerr << "allocated " << ngrid << " grids\n";
}

Grid::~Grid()
{
  fftw_free(fx);
  fftw_destroy_plan(plan_forward);
  fftw_destroy_plan(plan_inverse);
}

void Grid::fft_forward()
{
  // fk= V/nc^3 sum_x f(x) exp(-ik x)
  
  assert(mode == fft_mode_x);
  fftw_execute(plan_forward);
  mode = fft_mode_k;

  const double fac= pow(boxsize/nc, 3.0);
  const size_t ngrid= (size_t) nc*nc*2*(nc/2 + 1);
  for(int i=0; i<ngrid; ++i) {
    fx[i] *= fac;
  }
}

void Grid::fft_inverse()
{
  // fx = 1/V sum_k f(k) exp(ik x)
  assert(mode == fft_mode_k);
  fftw_execute(plan_inverse);
  mode= fft_mode_x;

  const double fac= 1.0/(boxsize*boxsize*boxsize);
  const size_t ngrid= (size_t) nc*nc*2*(nc/2 + 1);
  for(int i=0; i<ngrid; ++i) {
    fx[i] *= fac;
  }
}

void Grid::clear()
{
  memset(fx, 0, sizeof(double)*(nc*nc*ncz));
}

void Grid::print()
{
  if(mode == fft_mode_x)
    print_x();
  else if(mode == fft_mode_k)
    print_k();
  else
    assert(false);
}

void Grid::print_k()
{
  const size_t nckz= nc/2 + 1;

  for(int ix=0; ix<nc; ++ix) {
    double kx= ix < nc/2 ? ix : ix - nc;
    for(int iy=0; iy<nc; ++iy) {
      double ky= iy < nc/2 ? iy : iy - nc;
      for(int iz=0; iz<nckz; ++iz) {
	double kz= iz;
	double k= sqrt(kx*kx + ky*ky + kz*kz);
	size_t index= (ix*nc + iy)*nckz + iz;
	double delta_re= fk[index][0];
	double delta_im= fk[index][1];
	double Pk= delta_re*delta_re + delta_im*delta_im;
	printf("%e %e %e %e\n", k, Pk, fk[index][0], fk[index][1]);
	//printf("%d %d %d %e\n", ix, iy, iz, fk[index][0]);
      }
    }
  }
}

void Grid::print_x()
{
  const size_t nck= 2*(nc/2 + 1);

  for(int ix=0; ix<nc; ++ix) {
    double rx= ix < nc/2 ? ix : ix - nc;
    for(int iy=0; iy<nc; ++iy) {
      double ry= iy < nc/2 ? iy : iy - nc;
      for(int iz=0; iz<nc; ++iz) {
	double rz= iz < nc/2 ? iz : iz - nc;
	double r= sqrt(rx*rx + ry*ry + rz*rz);
	size_t index= (ix*nc + iy)*ncz + iz;
	printf("%d %d %d %e %e\n", ix, iy, iz, r, fx[index]);
	//if(isnan(fx[index]))
	//  printf("%d %d %d %e\n", ix, iy, iz, fx[index]);
      }
    }
  }
}
