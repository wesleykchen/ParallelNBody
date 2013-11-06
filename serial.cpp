#include "Util.hpp"

#include "Vec.hpp"

#include "kernel/Laplace.kern"

// Serial version of n-body algorithm
int main(int argc, char** argv)
{
  // Create a Kernel
  typedef LaplacePotential kernel_type;
  kernel_type K;

  static_assert(std::is_same<
                typename kernel_type::source_type,
                typename kernel_type::target_type>::value,
                "source_type != target_type");

  typedef typename kernel_type::source_type Point;
  std::vector<Point> data;
  typedef typename kernel_type::charge_type charge_type;
  std::vector<charge_type> sigma;

  // TODO: Generalize on source_type, charge_type, and result_type
  static_assert(std::is_same<typename kernel_type::source_type, Point>::value,
                "source_type != Vec<3,double>");
  static_assert(std::is_same<typename kernel_type::charge_type, double>::value,
                "charge_type != double");
  static_assert(std::is_same<typename kernel_type::result_type, double>::value,
                "result_type != double");

  std::vector<std::string> arg(argv, argv+argc);

  if (arg.size() < 3) {
    std::cerr << "Usage: " << arg[0] << " PHI_FILE SIGMA_FILE" << std::endl;
    //exit(1);
    // XXX: Remove
    std::cerr << "Using default " << PHIDATA << " " << SIGMADATA << std::endl;

    arg.resize(1);
    arg.push_back(PHIDATA);
    arg.push_back(SIGMADATA);
  }

  // Read the data from PHI_FILE interpreted as Points
  std::ifstream data_file(arg[1]);
  data_file >> data;

  // Read the data from SIGMA_FILE interpreted as doubles
  std::ifstream sigma_file(arg[2]);
  sigma_file >> sigma;

  // Make sure this makes sense and get metadata
  assert(data.size() == sigma.size());
  unsigned N = sigma.size();
  std::cout << "N = " << N << std::endl;

  // Compute the matvec
  typedef typename kernel_type::result_type result_type;
  std::vector<result_type> phi(N);

  Clock timer;
  timer.start();
  P2P::block_eval(K,
                  data.begin(), data.end(), sigma.begin(),
                  data.begin(), data.end(), phi.begin());
  double time = timer.elapsed();

  std::cout << "Computed in " << time << " seconds" << std::endl;
  double checkSum = std::accumulate(phi.begin(), phi.end(), 0.0);
  std::cout << "Serial - checksum answer is: " << checkSum << std::endl;

  std::ofstream phi_file("data/phi.txt");
  phi_file << phi << std::endl;
}
