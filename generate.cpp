#include "Util.hpp"
#include "Vec.hpp"

int main(int argc, char** argv)
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " PHI_FILE SIGMA_FILE ###" << std::endl;
    //exit(1);
    // XXX: Remove
    std::cerr << "Using default "
              << PHIDATA << " " << SIGMADATA << " " << NUMPOINTS << std::endl;
    argc = 4;
    char** new_argv = new char*[4];
    new_argv[0] = argv[0];
    new_argv[1] = (char*)& PHIDATA;
    new_argv[2] = (char*)& SIGMADATA;
    new_argv[3] = to_char(NUMPOINTS);
    argv = new_argv;
  }

  srand(time(NULL));
  unsigned N = atoi(argv[3]);

  typedef Vec<3,double> source_type;
  std::ofstream data(argv[1]);
  for (unsigned i = 0; i < N; ++i)
    data << source_type(get_random(), get_random(), get_random()) << std::endl;

  typedef double        charge_type;
  std::ofstream sigma(argv[2]);
  for (unsigned i = 0; i < N; ++i)
    sigma << charge_type(get_random()) << std::endl;
}
