#include "Util.hpp"

#include "kernel/NonParaBayesian.kern" 
#include "meta/kernel_traits.hpp"

// Serial version of n-body algorithm
int main(int argc, char** argv)
{
  bool checkErrors = true;

  // Parse optional command line args
  std::vector<std::string> arg(argv, argv + argc);
  for (unsigned i = 1; i < arg.size(); ++i) {
    if (arg[i] == "-nocheck") {
      checkErrors = false;
      arg.erase(arg.begin() + i);  // Erase this arg
      --i;                         // Reset index
    }
  }

  if (arg.size() < 3) {
    std::cerr << "Usage: " << arg[0] << " SOURCE_FILE CHARGE_FILE" << std::endl;
    //exit(1);
    // XXX: Remove
    std::cerr << "Using default " << SOURCE_DATA << " " << CHARGE_DATA << std::endl;

    arg.resize(1);
    arg.push_back(SOURCE_DATA);
    arg.push_back(CHARGE_DATA);
  }

  // Create a Kernel
  typedef NonParaBayesian kernel_type;
  kernel_type K(1,1);

  // Define source_type, target_type, charge_type, result_type
  typedef kernel_type::source_type source_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::result_type result_type;

  // We are testing symmetric kernels
  static_assert(std::is_same<source_type, target_type>::value,
                "Testing symmetric kernels, need source_type == target_type");

  std::vector<source_type> source;
  std::vector<charge_type> charge;

  // Read the data from SOURCE_FILE interpreted as Points
  std::ifstream source_file(arg[1]);
  source_file >> source;

  // Read the data from CHARGE_FILE interpreted as doubles
  std::ifstream charge_file(arg[2]);
  charge_file >> charge;

  // Make sure this makes sense and get metadata
  assert(source.size() == charge.size());
  unsigned N = charge.size();
  std::cout << "N = " << N << std::endl;

  // Compute the matvec
  std::vector<result_type> result(N);

  Clock timer;
  timer.start();
  p2p(K, source.begin(), source.end(), charge.begin(), result.begin());
  double time = timer.elapsed();

  std::cout << "Computed in " << time << " seconds" << std::endl;

  if (checkErrors) {
    std::cout << "Computing direct matvec..." << std::endl;

    std::vector<result_type> exact(N);

    // Compute the result with a direct matrix-vector multiplication
    p2p(K, source.begin(), source.end(), charge.begin(), exact.begin());

    print_error(exact, result);
  }

  std::ofstream result_file("data/result.txt");
  result_file << result << std::endl;
}
