#include "Util.hpp"

#include "kernel/Laplace.kern"
#include "meta/kernel_traits.hpp"

#include "meta/random.hpp"

int main(int argc, char** argv)
{
  std::vector<std::string> arg(argv, argv+argc);

  if (arg.size() < 3) {
    std::cerr << "Usage: " << arg[0] << " SOURCE_FILE CHARGE_FILE N" << std::endl;
    //exit(1);
    // XXX: Remove
    std::cerr << "Using default "
              << SOURCE_DATA << " " << CHARGE_DATA << " " << NUMPOINTS << std::endl;

    arg.resize(1);
    arg.push_back(SOURCE_DATA);
    arg.push_back(CHARGE_DATA);
    arg.push_back(to_string(NUMPOINTS));
  }

  srand(time(NULL));
  unsigned N = string_to_<int>(arg[3]);

  // Define the Kernel to use
  typedef LaplacePotential kernel_type;

  // Define source_type, target_type, charge_type, result_type
  typedef kernel_type::source_type source_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::result_type result_type;

  // We are testing symmetric kernels
  static_assert(std::is_same<source_type, target_type>::value,
                "Testing symmetric kernels, need source_type == target_type");

  std::ofstream source(arg[1]);
  for (unsigned i = 0; i < N; ++i)
    source << meta::random<source_type>::get() << std::endl;

  std::ofstream charge(arg[2]);
  for (unsigned i = 0; i < N; ++i)
    charge << meta::random<charge_type>::get() << std::endl;
}
