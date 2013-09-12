#include "Util.hpp"

#include "Vec.hpp"

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

  typedef Vec<3,double> Point;
  std::vector<Point> data;
  std::vector<double> sigma;
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

    // Read the data from PHI_FILE interpreted as Points
    std::ifstream data_file(arg[1]);
    data_file >> data;

    // Read the data from SIGMA_FILE interpreted as doubles
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
  MPI_Bcast(&N, 1, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  // Allocate memory on all other processes
  if (rank != MASTER) {
    data = std::vector<Point>(N);
    sigma = std::vector<double>(N);
  }

  // Broadcast the data to all processes
  unsigned count = Point::size() * N;
  commTimer.start();
  MPI_Bcast(data.data(), count, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(sigma.data(), N, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  // all processors have a chunk to hold their temporary answers
  std::vector<double> myphi(idiv_up(N,P));

  // evaluate computation
  block_eval(data.begin(), data.end(), sigma.begin(),
             data.begin() + calcStart(rank,P,N),
             data.begin() + calcEnd(rank,P,N),
             myphi.begin());

  // Collect results and display
  std::vector<double> phi;
  if (rank == MASTER)
    phi = std::vector<double>(N);

  commTimer.start();
  MPI_Gather(myphi.data(), idiv_up(N,P), MPI_DOUBLE,
             phi.data(), idiv_up(N,P), MPI_DOUBLE,
             MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  double time = timer.elapsed();
  printf("[%d] Timer: %e\n", rank, time);
  printf("[%d] CommTimer: %e\n", rank, totalCommTime);

  if (rank == MASTER) {
    double checkSum = std::accumulate(phi.begin(), phi.end(), 0.0);
    std::cout << "Broadcast - checksum answer is: " << checkSum << std::endl;
    std::ofstream phi_file("data/phi.txt");
    phi_file << phi << std::endl;
  }

  MPI_Finalize();
  return 0;
}
