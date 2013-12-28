#include "Util.hpp"

#include "kernel/Laplace.kern"
#include "meta/kernel_traits.hpp"

// Serial version of n-body algorithm
int main(int argc, char** argv)
{
  // Create a Kernel
  typedef LaplacePotential kernel_type;
  kernel_type K;

  // Define source_type, target_type, charge_type, result_type
  typedef kernel_type::source_type source_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::result_type result_type;

  // We are testing symmetric kernels
  static_assert(std::is_same<source_type, target_type>::value,
                "Testing symmetric kernels, need source_type == target_type");

  std::vector<source_type> data;
  std::vector<charge_type> sigma;

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
  std::vector<result_type> phi(N);

  Clock timer;
  timer.start();
  p2p(K,
      data.begin(), data.end(),
      sigma.begin(), phi.begin());
  double time = timer.elapsed();

  std::cout << "Computed in " << time << " seconds" << std::endl;
  result_type check = std::accumulate(phi.begin(), phi.end(), result_type());
  std::cout << "Serial - checksum answer is: " << check << std::endl;

  std::ofstream phi_file("data/phi.txt");
  phi_file << phi << std::endl;
}
