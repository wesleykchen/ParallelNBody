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

  std::vector<source_type> data;
  std::vector<charge_type> sigma;
  unsigned N;

  if (rank == MASTER) {
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

    // Read the data from PHI_FILE interpreted as source_types
    std::ifstream data_file(arg[1]);
    data_file >> data;

    // Read the data from SIGMA_FILE interpreted as charge_types
    std::ifstream sigma_file(arg[2]);
    sigma_file >> sigma;

    assert(data.size() == sigma.size());
    N = sigma.size();
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
    data  = std::vector<source_type>(N);
    sigma = std::vector<charge_type>(N);
  }

  // Broadcast the data to all processes
  commTimer.start();
  MPI_Bcast(data.data(), sizeof(source_type) * data.size(), MPI_CHAR, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(sigma.data(), sizeof(charge_type) * sigma.size(), MPI_CHAR, MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  // all processors have a chunk to hold their temporary answers
  std::vector<result_type> myphi(idiv_up(N,P));

  // Evaluate computation
  p2p(K,
      data.begin(), data.end(), sigma.begin(),
      data.begin() + calcStart(rank,P,N),
      data.begin() + calcEnd(rank,P,N),
      myphi.begin());

  // Collect results and display
  std::vector<result_type> phi;
  if (rank == MASTER)
    phi = std::vector<result_type>(N);

  commTimer.start();
  MPI_Gather(myphi.data(), sizeof(result_type) * myphi.size(), MPI_CHAR,
             phi.data(), sizeof(result_type) * myphi.size(), MPI_CHAR,
             MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  double time = timer.elapsed();
  printf("[%d] Timer: %e\n", rank, time);
  printf("[%d] CommTimer: %e\n", rank, totalCommTime);

  if (rank == MASTER) {
    result_type check = std::accumulate(phi.begin(), phi.end(), result_type());
    std::cout << "Broadcast - checksum answer is: " << check << std::endl;
    std::ofstream phi_file("data/phi.txt");
    phi_file << phi << std::endl;
  }

  MPI_Finalize();
  return 0;
}
