#include "Util.hpp"

#include "kernel/NonParaBayesian.kern"
#include "meta/kernel_traits.hpp"
#include "meta/random.hpp"

#include <type_traits>

// Broadcast version of n-body algorithm

inline unsigned calcStart(unsigned r, unsigned P, unsigned N) {
  return std::min(N, r * idiv_up(N,P));
}

inline unsigned calcEnd(unsigned r, unsigned P, unsigned N) {
  return calcStart(r+1, P, N);
}


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

  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int P;
  MPI_Comm_size(MPI_COMM_WORLD, &P);

  // Define the kernel
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

  if (rank == MASTER) {
    // generate source data
    for (unsigned i = 0; i < N; ++i)
      source.push_back(meta::random<source_type>::get());

    // generate charge data
    for (unsigned i = 0; i < N; ++i)
      charge.push_back(meta::random<charge_type>::get());

    // display metadata
    std::cout << "N = " << N << std::endl;
    std::cout << "P = " << P << std::endl;
  }

  Clock timer;
  Clock commTimer;
  double totalCommTime = 0;

  // Broadcast the size of the problem to all processes
  timer.start();
  commTimer.start();
  MPI_Bcast(&N, sizeof(N), MPI_CHAR, MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  // Allocate memory on all other processes
  if (rank != MASTER) {
    source = std::vector<source_type>(N);
    charge = std::vector<charge_type>(N);
  }

  // Broadcast the data to all processes
  commTimer.start();
  MPI_Bcast(source.data(), sizeof(source_type) * source.size(), MPI_CHAR,
            MASTER, MPI_COMM_WORLD);
  MPI_Bcast(charge.data(), sizeof(charge_type) * charge.size(), MPI_CHAR,
            MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  // All processors have a chunk to hold their temporary answers
  std::vector<result_type> rI(idiv_up(N,P));

  // Evaluate computation
  p2p(K,
      source.begin(), source.end(), charge.begin(),
      source.begin() + calcStart(rank,P,N),
      source.begin() + calcEnd(rank,P,N),
      rI.begin());

  // Collect results and display
  std::vector<result_type> result;
  if (rank == MASTER)
    result = std::vector<result_type>(N);

  commTimer.start();
  MPI_Gather(rI.data(), sizeof(result_type) * rI.size(), MPI_CHAR,
             result.data(), sizeof(result_type) * rI.size(), MPI_CHAR,
             MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  double time = timer.elapsed();
  printf("[%d] Timer: %e\n", rank, time);
  printf("[%d] CommTimer: %e\n", rank, totalCommTime);

  // Check the result
  if (rank == MASTER && checkErrors) {
    std::cout << "Computing direct matvec..." << std::endl;

    std::vector<result_type> exact(N);

    compTimer.starT();
    // Compute the result with a direct matrix-vector multiplication
    p2p(K, source.begin(), source.end(), charge.begin(), exact.begin());
    double directCompTime = compTime.elapsed();

    print_error(exact, result);
    std::cout << "DirectCompTime: " << directCompTime << std::endl;
  }

  if (rank == MASTER) {
    std::ofstream result_file("data/result.txt");
    result_file << result << std::endl;
  }

  MPI_Finalize();
  return 0;
}
