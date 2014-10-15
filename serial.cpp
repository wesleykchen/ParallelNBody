#include "Util.hpp"

#include "meta/kernel_traits.hpp"
#include "meta/random.hpp"

#include "kernel/InvSq.kern"

// Serial version of n-body algorithm
int main(int argc, char** argv)
{
  bool checkErrors = true;

  // Parse optional command line args
  std::vector<std::string> arg(argv, argv + argc);
  for (unsigned i = 1; i < arg.size(); ++i) {
    if (arg[i] == "-nocheck") {
      checkErrors = false;
      arg.erase(arg.begin() + i, arg.begin() + i + 1);  // Erase this arg
      --i;                                              // Reset index
    }

    if (arg.size() < 2) {
      std::cerr << "Usage: " << arg[0] << " NUMPOINTS [-nocheck]" << std::endl;
      exit(1);
    }
  }

  srand(time(NULL));
  unsigned N = string_to_<int>(arg[1]);

  // Create a Kernel
  typedef InvSq kernel_type;
  kernel_type K;

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

  // Set the seed
  const int seed = 1337;
  meta::default_generator.seed(1337);

  // generate source data
  for (unsigned i = 0; i < N; ++i)
    source.push_back(meta::random<source_type>::get());

  // generate charge data
  for (unsigned i = 0; i < N; ++i)
    charge.push_back(meta::random<charge_type>::get());

  // Display metadata
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
    timer.start();
    detail::block_eval(K,
                       source.begin(), source.end(),
                       charge.begin(), exact.begin());
    double directCompTime = timer.elapsed();

    print_error(exact, result);
    std::cout << "DirectCompTime: " << directCompTime << std::endl;
  }

  std::ofstream result_file("data/result128k.txt");

  std::string result_filename = "data/";
  result_filename += std::string("invsq")
      + "_n" + std::to_string(N)
      + "_s" + std::to_string(seed) + ".txt";

  std::ofstream result_file(result_filename);
}
