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
  int numtasks;
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  if (argc < 2) {
    if (rank == MASTER) {
      std::cerr << "Usage: " << argv[0] << " PHI_FILE SIGMA_FILE" << std::endl;
      //exit(1);
      // XXX: Remove
      std::cerr << "Using default " << PHIDATA << " " << SIGMADATA << std::endl;
    }
    argc = 3;
    char** new_argv = new char*[3];
    new_argv[0] = argv[0];
    new_argv[1] = (char*)& PHIDATA;
    new_argv[2] = (char*)& SIGMADATA;
    argv = new_argv;
  }

  typedef Vec<3,double> Point;
  std::vector<Point> data;
  std::vector<double> sigma;
  unsigned N;

  Clock timer;
  Clock commTimer;
  double totalCommTime = 0;

  // Read data on MASTER
  if (rank == MASTER) {
    // Read the data from PHI_FILE interpreted as Points
    std::ifstream data_file(argv[1]);
    data_file >> data;

    // Read the data from SIGMA_FILE interpreted as doubles
    std::ifstream sigma_file(argv[2]);
    sigma_file >> sigma;

    assert(data.size() == sigma.size());
    N = sigma.size();
    std::cout << "N = " << N << std::endl;
  }

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
  MPI_Bcast(&data[0], count, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&sigma[0], N, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  // all processors have a chunk to hold their temporary answers
  std::vector<double> myphi(idiv_up(N,numtasks));

  // evaluate computation
  block_eval(data.begin(), data.end(), sigma.begin(),
             data.begin() + calcStart(rank,numtasks,N),
             data.begin() + calcEnd(rank,numtasks,N),
             myphi.begin());

  // Collect results and display
  std::vector<double> phi;
  if (rank == MASTER)
    phi = std::vector<double>(N);

  commTimer.start();
  MPI_Gather(&myphi[0], idiv_up(N,numtasks), MPI_DOUBLE,
             &phi[0], idiv_up(N,numtasks), MPI_DOUBLE,
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
