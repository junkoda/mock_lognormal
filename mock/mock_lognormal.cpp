//
// Using Boost program options
//   style:   ./options [options] <required arg>
//   example: ./options --x=3 filename
//

#include <iostream>
#include <string>
#include <valarray>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <boost/program_options.hpp>


#include "power.h"
#include "cosmology.h"
#include "grid.h"


using namespace std;
using namespace boost::program_options;

void set_power_spectrum(PowerSpectrum const * const ps, const double boxsize, Grid* const grid);

void set_lognormal_xi(Grid* const grid);

void generate_delta_k(const unsigned long seed,
		      const bool fix_amplitude,
		      const int n_mas,
		      Grid* const grid);

double generate_lognormal_density(Grid* const grid);

void print_mock(const unsigned long seed, const size_t np,
		const double boxsize, const double n_max,
		Grid const * const grid);

static void print_pk(Grid const * const grid);
static void print_xi(Grid const * const grid);

int main(int argc, char* argv[])
{
  //
  // command-line options (Boost program_options)
  //
  options_description opt("options [options] filename");
  opt.add_options()
    ("help,h", "display this help")
    ("filename,f", value<string>(), "filename")
    ("nc", value<int>()->default_value(64), "number of grids per dimension")
    ("boxsize", value<double>()->default_value(0.0), "boxsize")
    ("np", value<size_t>()->default_value(1000), "number of output particles")
    //("z", value<double>()->default_value(0.0), "redshift")
    //("omega-m", value<double>()->default_value(0.308), "Omega matter")
    ("fix-amplitude", "fix delta_k amplitude")
    ("seed", value<unsigned long>()->default_value(1), "random seed")
    ("corrct-mas", "fix mass assignment window function")
    ("print-Pk", "print intput P(k)")
    ("print-xi", "print xi(r), FFT of intput P(k)")
    ("print-Pkg", "print P_gaussian(k), FFT of intput P(k)")
    ("test-fft", "print P(k) after forward and inverse FFT")

    ;
  
  positional_options_description p;
  p.add("filename", -1);
  
  variables_map vm;
  store(command_line_parser(argc, argv).options(opt).positional(p).run(), vm);
  notify(vm);

  if(vm.count("help") || ! vm.count("filename")) {
    cout << opt; 
    return 0;
  }

  const int nc= vm["nc"].as<int>(); assert(nc > 0);
  //const double omega_m= vm["omega-m"].as<double>();
  //const double z= vm["z"].as<double>();
  const double boxsize= vm["boxsize"].as<double>();
  const bool fix_amplitude= vm.count("fix-amplitude");
  const unsigned long seed= vm["seed"].as<unsigned long>();
  
  const string filename= vm["filename"].as<string>();

  const int n_mas= 0.0; // DEBUG!!!

  cerr << "fix-amplitude= " << fix_amplitude << endl; 

  // Load power spectrum
  PowerSpectrum* const ps= power_alloc(filename.c_str());

  // Set power spectrum amplitude
  //cosmology_init(omega_m);
  //const double a= 1.0/(1.0 + z); 
  //const double growth= cosmology_D_growth(a);
  //power_rescale(ps, growth);

  Grid* grid= new Grid(nc, boxsize);

  //  P(k)
  cerr << "generate target P(k)\n";
  set_power_spectrum(ps, boxsize, grid); // grid => P(k)

  if(vm.count("print-Pk")) {
    cout << "# P(k) grid\n";
    print_pk(grid);
    return 0;
  }

  // xi(r)
  cerr << "FFT to xi(r)\n";
  
  grid->fft_inverse(); // grid => xi(r)

  if(vm.count("print-xi")) {
    print_xi(grid);
    cerr << "Print xi\n";
    return 0;
  }

  if(vm.count("test-fft")) {
    cout << "# FFT test; k P(k)\n";
    grid->fft_forward();
    grid->print();
    return 0;
  }

  cerr << "convert to gaussian xi = log(1 + xi_target)\n";
  //set_lognormal_xi(grid); // grid => xi_gaussian(r)
  //DEBUG!!!

  // grid of P_gaussian(k)
  cerr << "FFT to P_gaussian(k)\n";
  grid->fft_forward(); // grid => P_gaussian(r)

  if(vm.count("print-Pkg")) {
    cout << "# P_gaussian(k)\n";
    print_pk(grid);
    return 0;
  }

  cerr << "Generate random delta_k\n";
  generate_delta_k(seed, fix_amplitude, n_mas, grid);
  // grid => delta_gaussian(k)

  cerr << "FFT to lognormal density 1 + delta(x)\n";
  grid->fft_inverse(); // grid = 1 + delta_lognormal(x)

  double n_max= generate_lognormal_density(grid);

  cerr << "Print lognormal mock\n";
  cerr << "n_max= " << n_max << endl;
  const size_t np= vm["np"].as<size_t>();

  print_mock(seed + 100, np, boxsize, n_max, grid);

  
  return 0;
}

void set_power_spectrum(PowerSpectrum const * const ps, const double boxsize, Grid* const grid)
{
  // Output: grid of P(k)
  grid->clear();
  int nc= grid->nc;
  const int nckz= nc/2 + 1;

  fftw_complex* const pk= (fftw_complex*) grid->fx; 

  double fac= 2.0*M_PI/boxsize;
  for(int ix=0; ix<nc; ++ix) {
   double kx= ix < nc/2 ? fac*ix : fac*(ix - nc);
   for(int iy=0; iy<nc; ++iy) {
    double ky= iy < nc/2 ? fac*iy : fac*(iy - nc);
    int iz0= (ix == 0 && iy == 0);

    for(int iz=0; iz<nc/2; ++iz) {
      double kz= fac*iz;
      double k= sqrt(kx*kx + ky*ky + kz*kz);
      size_t index= (ix*nc + iy)*nckz + iz;
	  
      pk[index][0]= power(ps, k);
      pk[index][1]= 0.0;

      //printf("%e %e\n", k, pk[index][0]);
    }
   }
  }

  pk[0][0]= pk[0][1]= 0.0;

  grid->mode= fft_mode_k;
}

void set_lognormal_xi(Grid* const grid)
{
  // Convert non-linear xi(r) to gaussian xi_g(r)
  // xi_g = log(1 + xi(r))
  double* xi= grid->fx;
  size_t nc= grid->nc;
  size_t ncz= 2*(nc/2 + 1);
  
  for(int ix=0; ix<nc; ++ix) {
    for(int iy=0; iy<nc; ++iy) {
      size_t index= (ix*nc + iy)*ncz;
      for(int iz=0; iz<nc; ++iz) {
	assert(1.0 + xi[index] > 0.0);
	xi[index]= log(1.0 + xi[index]);
	index++;
      }
    }
  }
}


void compute_window_array(valarray<double>& v, const int n_mas)
{
  const int nc= v.size();
  const int knq = nc/2;
  const double fac= M_PI/nc;

  for(int i=1; i<nc; ++i) {
    int k= i <= knq ? i : i - nc;
    double sinc = sin(fac*k)/(fac*k);
    v[i] = 1.0/pow(sinc, 2*n_mas);
  }
}

  
void generate_delta_k(const unsigned long seed,
		      const bool fix_amplitude,
		      const int n_mas,
		      Grid* const grid) 
{
  // input grid: P_gauusian(k)
  // output grid: delta_gaussian_k(k)
  const int nc= grid->nc;
  const size_t nckz= nc/2 + 1;
  valarray<double> v_corr(1.0, nc);

  gsl_rng* rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed);

  size_t negative= 0;
  double P_min= 0.0;


  //const double fft_fac= 1.0/(nc*nc*nc);

  fftw_complex* const fk= (fftw_complex*) grid->fx;
  
  if(n_mas > 0)
    compute_window_array(v_corr, n_mas);

  //double fac= 2.0*M_PI/boxsize;
  for(int ix=0; ix<nc; ++ix) {
    double kx= ix < nc/2 ? ix : ix - nc;
    //double kx= ix < nc/2 ? fac*ix : fac*(ix - nc);
   double corr_x= v_corr[ix];
   for(int iy=0; iy<nc; ++iy) {
     double ky= iy < nc/2 ? iy : iy - nc;
     //double ky= iy < nc/2 ? fac*iy : fac*(iy - nc);
    double corr_xy= corr_x*v_corr[iy];
    size_t index= (ix*nc + iy)*nckz;
    for(int iz=0; iz<nc/2; ++iz) {
      double kz= iz;
      // double kz= fac*iz;
      double k= sqrt(kx*kx + ky*ky + kz*kz);
      double corr_xyz= corr_xy*v_corr[iz];

      double phase= gsl_rng_uniform(rng)*M_PI;

      double delta2= fk[index][0]/corr_xyz;
      

      if(!fix_amplitude) {
	double ampl;
	do
	  ampl = gsl_rng_uniform(rng);
	while(ampl == 0.0);
	
	delta2 *= -log(ampl);
      }

      double delta_k_mag= 0.0;
      if(fk[index][0] < P_min)
	P_min= fk[index][0];
	  
      if( fk[index][0] > 0.0) {
	//fprintf(stderr, "%e %e\n", k, delta2);
	delta_k_mag= sqrt(delta2);
      }
      else {
	negative++;
      }

      fk[index][0]= delta_k_mag*cos(phase);
      fk[index][1]= delta_k_mag*sin(phase);
	  
      index++;

    }
   }
  }

  gsl_rng_free(rng);

  fprintf(stderr, "P_min= %e\n", P_min);
  fprintf(stderr, "negative P(k): %zu\n", negative);
}

double generate_lognormal_density(Grid* const grid)
{
  // input: delta_gaussian(x)
  // output: 1 + delta_lognormal \prop exp(delta_gaussian)

  assert(grid->mode == fft_mode_x);
  size_t nc= grid->nc;
  size_t ncz= 2*(nc/2 + 1);
  double* const delta= grid->fx;
  long double nsum= 0.0;
  double nmax= 0.0;
  
  for(int ix=0; ix<nc; ++ix) {
   for(int iy=0; iy<nc; ++iy) {
    size_t index= (ix*nc + iy)*ncz;
    for(int iz=0; iz<nc; ++iz) {
      double d= delta[index]= exp(delta[index]);
      nsum += d;
      if(d > nmax) nmax= d;

      //cerr << delta[index] << " " << d << endl;
      
      index++;
    }
   }
  }

  // normalization
  const double fac= nsum/(nc*nc*nc);
  for(size_t ix=0; ix<nc; ++ix) {
   for(size_t iy=0; iy<nc; ++iy) {
    size_t index= (ix*nc + iy)*ncz;
    for(int iz=0; iz<nc; ++iz) {
      delta[index] /= fac;
      
      index++;
    }
   }
  }

  return nmax/fac;
}
  
void print_mock(const unsigned long seed, const size_t np,
		const double boxsize, const double n_max,
		Grid const * const grid)
{
  // input: grid of 1 + delta

  gsl_rng* rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed);

  const size_t nc= grid->nc;
  const size_t ncz= nc/2 + 1;
  const double nmean= np/(nc*nc*nc);
  double const * const n= grid->fx;

  double x[3];
  int ix[3];

  const size_t niter= static_cast<size_t>(n_max*np);
  cerr << "n_iteration= " << niter << endl;
  //cerr << niter << endl;
  size_t count= 0;
				      
  for(size_t i=0; i<niter; ++i) {
    for(int k=0; k<3; ++k) {
      x[k]= boxsize*gsl_rng_uniform(rng);
      ix[k]= (int) floor(x[k]/boxsize*nc);
      ix[k]= ix[k] % nc;
    }
    size_t index= (nc*ix[0] + ix[1])*ncz + ix[2];
    double prob= n[index] / n_max;

    if(gsl_rng_uniform(rng) < prob) {
      printf("%e %e %e\n", x[0], x[1], x[2]);
      count++;
    }

  }
    
  gsl_rng_free(rng);

  fprintf(stderr, "np= %zu, count= %zu\n", np, count);
}

void print_xi(Grid const * const grid)
{
  const double boxsize= grid->boxsize;
  const int nc= grid->nc;
  const size_t ncz= 2*(nc/2 + 1);
  const double fac= boxsize/nc;

  for(int ix=0; ix<nc; ++ix) {
   double rx= ix <= nc/2 ? fac*ix : fac*(ix - nc);
   for(int iy=0; iy<nc; ++iy) {
    double ry= iy <= nc/2 ? fac*iy : fac*(iy - nc);
    for(int iz=0; iz<nc; ++iz) {
      double rz= iz <= nc/2 ? fac*iz : fac*(iz - nc);
      double r= sqrt(rx*rx + ry*ry + rz*rz);

      size_t index= (ix*nc + iy)*ncz + iz;
      printf("%e %e\n", r, grid->fx[index]);
    }
   }
  }
}

void print_pk(Grid const * const grid)
{
  const int nc= grid->nc;
  const size_t nckz= nc/2 + 1;
  const double fac= 2.0*M_PI/grid->boxsize;
  
  for(int ix=0; ix<nc; ++ix) {
    double kx= ix < nc/2 ? fac*ix : fac*(ix - nc);
    for(int iy=0; iy<nc; ++iy) {
      double ky= iy < nc/2 ? fac*iy : fac*(iy - nc);
      for(int iz=0; iz<nckz; ++iz) {
	double kz= fac*iz;
	double k= sqrt(kx*kx + ky*ky + kz*kz);
	size_t index= (ix*nc + iy)*nckz + iz;

	printf("%e %e %e\n", k, grid->fk[index][0], grid->fk[index][1]);
      }
    }
  }
}
  
