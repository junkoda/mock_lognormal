//
// Using Boost program options
//   style:   ./options [options] <required arg>
//   example: ./options --x=3 filename
//

#include <iostream>
#include <string>
#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;



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
  
  const string filename= vm["filename"].as<string>();
  cout << "filename= " << filename << endl;


  Grid* grid= Grid(nc);
  
  plan= fftwf_plan_dft_r2c_3d(nc, nc, nc, mesh, (fftwf_complex*) mesh, 
			      FFTW_ESTIMATE);
    
  memset(mesh, 0, sizeof(float)*nmesh);

  
  return 0;
}

void alloc
