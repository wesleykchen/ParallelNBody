#include "Util.hpp"

// Scatter version of the n-body algorithm

#include "kernel/Laplace.kern"
#include "meta/kernel_traits.hpp"

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int P;
  MPI_Comm_size(MPI_COMM_WORLD, &P);

  typedef LaplacePotential kernel_type;
  kernel_type K;

  // Define source_type, target_type, charge_type, result_type
  typedef kernel_type::source_type source_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::result_type result_type;

  std::vector<source_type> source;
  std::vector<charge_type> charge;
  unsigned N;

  Clock timer;
  Clock commTimer;
  double totalCommTime = 0;

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

    // Read the data from SOURCE_FILE interpreted as Points
    std::ifstream source_file(arg[1]);
    source_file >> source;

    // Read the data from CHARGE_FILE interpreted as doubles
    std::ifstream charge_file(arg[2]);
    charge_file >> charge;

    assert(source.size() == charge.size());
    N = charge.size();
    std::cout << "N = " << N << std::endl;
    std::cout << "P = " << P << std::endl;

    // Pad with garbage values if needed
    source.resize(P * idiv_up(N,P));
    charge.resize(P * idiv_up(N,P));
  }

  // Broadcast the size of the problem to all processes
  timer.start();
  commTimer.start();
  MPI_Bcast(&N, sizeof(N), MPI_CHAR, MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  // TODO: Generalize by excluding the garbage values
  if (N % P != 0) {
    printf("Quitting. The number of processors must divide the total number of tasks.\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    exit(0);
  }

  // Allocate memory all processes
  std::vector<result_type> phiI(idiv_up(N,P));
  std::vector<source_type> xI(idiv_up(N,P));
  std::vector<source_type> xJ(idiv_up(N,P));
  std::vector<charge_type> sigmaJ(idiv_up(N,P));

  // Scatter the data to all processes
  commTimer.start();
  MPI_Scatter(source.data(), sizeof(source_type) * xI.size(), MPI_CHAR,
              xI.data(), sizeof(source_type) * xI.size(), MPI_CHAR,
              MASTER, MPI_COMM_WORLD);
  MPI_Scatter(charge.data(), sizeof(charge_type) * sigmaJ.size(), MPI_CHAR,
              sigmaJ.data(), sizeof(charge_type) * sigmaJ.size(), MPI_CHAR,
              MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  // Calculate the symmetric first block
  p2p(K,
      xI.begin(), xI.end(),
      sigmaJ.begin(), phiI.begin());

  // Copy xI -> xJ
  std::copy(xI.begin(), xI.end(), xJ.begin());

  MPI_Status status;
  for (int shiftCount = 1; shiftCount < P; ++shiftCount) {
    commTimer.start();
    // Add P to prevent negative numbers
    int prev  = (rank - 1 + P) % P;
    int next = (rank + 1 + P) % P;
    MPI_Sendrecv_replace(xJ.data(), sizeof(source_type) * xJ.size(),
                         MPI_CHAR, next, 0, prev, 0,
                         MPI_COMM_WORLD, &status);
    MPI_Sendrecv_replace(sigmaJ.data(), sizeof(charge_type) * sigmaJ.size(),
                         MPI_CHAR, next, 0, prev, 0,
                         MPI_COMM_WORLD, &status);
    totalCommTime += commTimer.elapsed();

    // Calculate the current block
    p2p(K,
        xJ.begin(), xJ.end(), sigmaJ.begin(),
        xI.begin(), xI.end(), phiI.begin());
  }

  std::vector<result_type> phi;
  if (rank == MASTER)
    phi = std::vector<result_type>(P*idiv_up(N,P));

  // Collect results and display
  commTimer.start();
  MPI_Gather(phiI.data(), sizeof(result_type) * phiI.size(), MPI_CHAR,
             phi.data(), sizeof(result_type) * phiI.size(), MPI_CHAR,
             MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  double time = timer.elapsed();
  printf("[%d] Timer: %e\n", rank, time);
  printf("[%d] CommTimer: %e\n", rank, totalCommTime);

  if (rank == MASTER) {
    result_type check = std::accumulate(phi.begin(), phi.end(), result_type());
    std::cout << "Scatter - checksum answer is: " << check << std::endl;
    std::ofstream phi_file("data/phi.txt");
    phi_file << phi << std::endl;
  }

  MPI_Finalize();
  return 0;
}
