#include "Util.hpp"

#include "kernel/Laplace.kern"
#include "meta/kernel_traits.hpp"

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
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int P;
  MPI_Comm_size(MPI_COMM_WORLD, &P);

  // Define the kernel
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

  std::vector<source_type> source;
  std::vector<charge_type> charge;
  unsigned N;

  if (rank == MASTER) {
    std::vector<std::string> arg(argv, argv+argc);

    if (arg.size() < 3) {
      std::cerr << "Usage: " << arg[0] << " SOURCE_FILE CHARGE_FILE" << std::endl;
      //exit(1);
      // XXX: Remove
      std::cerr << "Using default " << SOURCE_DATA << " " << CHARGE_DATA << std::endl;

      arg.resize(1);
      arg.push_back(SOURCE_DATA);
      arg.push_back(CHARGE_DATA);
    }

    // Read the data from SOURCE_FILE interpreted as source_types
    std::ifstream source_file(arg[1]);
    source_file >> source;

    // Read the data from CHARGE_FILE interpreted as charge_types
    std::ifstream charge_file(arg[2]);
    charge_file >> charge;

    assert(source.size() == charge.size());
    N = charge.size();
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
  MPI_Bcast(source.data(), sizeof(source_type) * source.size(),
            MPI_CHAR, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(charge.data(), sizeof(charge_type) * charge.size(),
            MPI_CHAR, MASTER, MPI_COMM_WORLD);
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

  if (rank == MASTER) {
    std::cout << "Broadcast - checksum answer is: "
              << std::accumulate(result.begin(), result.end(), result_type())
              << std::endl;
    std::ofstream result_file("data/result.txt");
    result_file << result << std::endl;
  }

  MPI_Finalize();
  return 0;
}
