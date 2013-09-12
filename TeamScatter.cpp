#include "Util.hpp"
#include "Vec.hpp"

// Team Scatter version of the n-body algorithm

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int P;
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  unsigned teamsize;

  typedef Vec<3,double> Point;
  std::vector<Point> data;
  std::vector<double> sigma;
  unsigned N;

  if (rank == MASTER) {
    std::vector<std::string> arg(argv, argv + argc);

    // Parse optional command line args
    teamsize = 1;   // Default value
    for (unsigned i = 1; i < arg.size(); ++i) {
      if (arg[i] == "-c") {
        if (i+1 < arg.size()) {
          teamsize = string_to_<unsigned>(arg[i+1]);
          arg.erase(arg.begin() + i, arg.begin() + i + 2);  // Erase these two
          --i;                                              // Reset index
        } else {
          std::cerr << "-c option requires one argument." << std::endl;
          return 1;
        }
      }
    }

    if (arg.size() < 2) {
      std::cerr << "Usage: " << arg[0] << " PHI_FILE SIGMA_FILE [-c TEAMSIZE]" << std::endl;
      //exit(1);
      // XXX: Remove
      std::cerr << "Using default " << PHIDATA << " " << SIGMADATA << std::endl;

      arg.resize(1);  // keep only the executable name
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
    std::cout << "Teamsize = " << teamsize << std::endl;
  }

  Clock timer;
  Clock commTimer;
  double totalCommTime = 0;

  timer.start();

  // Broadcast the size of the problem to all processes
  commTimer.start();
  MPI_Bcast(&N, 1, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  // Broadcast the teamsize to all processes
  commTimer.start();
  MPI_Bcast(&teamsize, 1, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  unsigned num_teams = P / teamsize;
  // Determine coordinates in processor team grid
  unsigned team  = rank / teamsize;
  unsigned trank = rank % teamsize;

  // TODO: How to generalize?
  if (N % P != 0) {
    printf("Quitting. The number of processors must divide the number of points\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    exit(0);
  }

  if (P % teamsize != 0) {
    printf("Quitting. The teamsize (c) must divide the total number of processors (p).\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    exit(0);
  }

  if (teamsize * teamsize > unsigned(P)) {
    printf("Quitting. The teamsize ^ 2 (c^2) must be less than or equal to the number of  processors (p).\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    exit(0);
  }

  // Split comm into row and column communicators
  MPI_Comm team_comm;
  MPI_Comm_split(MPI_COMM_WORLD, team, rank, &team_comm);
  MPI_Comm row_comm;
  MPI_Comm_split(MPI_COMM_WORLD, trank, rank, &row_comm);

  // Declare data for the block computations
  std::vector<Point> xI(idiv_up(N,num_teams));
  std::vector<Point> xJ(idiv_up(N,num_teams));
  std::vector<double> sigmaJ(idiv_up(N,num_teams));

  unsigned teamDataCount = Point::size() * N / num_teams;

  // Scatter data from master to team leaders
  if (trank == MASTER) {
    commTimer.start();
    MPI_Scatter(&data[0], teamDataCount, MPI_DOUBLE,
                &xI[0], teamDataCount, MPI_DOUBLE,
                MASTER, row_comm);
    MPI_Scatter(&sigma[0], teamDataCount / Point::size(), MPI_DOUBLE,
                &sigmaJ[0], teamDataCount / Point::size(), MPI_DOUBLE,
                MASTER, row_comm);
    totalCommTime += commTimer.elapsed();
  }

  // Team leaders broadcast to team
  commTimer.start();
  MPI_Bcast(&xI[0], teamDataCount, MPI_DOUBLE,
            MASTER, team_comm);
  MPI_Bcast(&sigmaJ[0], teamDataCount / Point::size(), MPI_DOUBLE,
            MASTER, team_comm);
  totalCommTime += commTimer.elapsed();

  // Copy xI -> xJ
  std::copy(xI.begin(), xI.end(), xJ.begin());

  // Perform initial offset by teamrank
  commTimer.start();
  // % is implementation defined, adding num_teams to prevent negative numbers
  int prev = (team - trank + num_teams) % num_teams;
  int next = (team + trank + num_teams) % num_teams;
  MPI_Status status;
  MPI_Sendrecv_replace(&xJ[0], teamDataCount, MPI_DOUBLE,
                       next, 0, prev, 0,
                       row_comm, &status);
  MPI_Sendrecv_replace(&sigmaJ[0], teamDataCount / Point::size(), MPI_DOUBLE,
                       next, 0, prev, 0,
                       row_comm, &status);
  totalCommTime += commTimer.elapsed();

  // Hold accumulated answers
  std::vector<double> phiI(idiv_up(N,num_teams));

  // Calculate the current block
  block_eval(xJ.begin(), xJ.end(), sigmaJ.begin(),
             xI.begin(), xI.end(), phiI.begin());

  // Looping process to shift the data between the teams
  int ceilPC2 = idiv_up(P, teamsize*teamsize);
  for (int shiftCount = 1; shiftCount < ceilPC2; ++shiftCount) {
    commTimer.start();
    // % is implementation defined, adding num_teams to prevent negative numbers
    int prev = (team - teamsize + num_teams) % num_teams;
    int next = (team + teamsize + num_teams) % num_teams;
    MPI_Sendrecv_replace(&xJ[0], teamDataCount, MPI_DOUBLE,
                         next, 0, prev, 0,
                         row_comm, &status);
    MPI_Sendrecv_replace(&sigmaJ[0], teamDataCount / Point::size(), MPI_DOUBLE,
                         next, 0, prev, 0,
                         row_comm, &status);
    totalCommTime += commTimer.elapsed();

    // Compute on the last iteration only if
    // 1) The teamsize divides the number of teams (everyone computes)
    // 2) Your team rank is one of the remainders
    if (shiftCount < ceilPC2-1
        || (num_teams % teamsize == 0 || trank < num_teams % teamsize)) {
      block_eval(xJ.begin(), xJ.end(), sigmaJ.begin(),
                 xI.begin(), xI.end(), phiI.begin());
    }
  }

  // Allocate teamphiI on team leaders
  std::vector<double> teamphiI(1);
  if (trank == MASTER)
    teamphiI = std::vector<double>(idiv_up(N,num_teams));

  // Reduce answers to the team leader
  commTimer.start();
  MPI_Reduce(&phiI[0], &teamphiI[0], teamDataCount / Point::size(), MPI_DOUBLE,
             MPI_SUM, MASTER, team_comm);
  totalCommTime += commTimer.elapsed();

  // Allocate phi on master
  std::vector<double> phi(1);
  if (rank == MASTER)
    phi = std::vector<double>(P*idiv_up(N,P));

  // Gather team leader answers to master
  if (trank == MASTER) {
    commTimer.start();
    MPI_Gather(&teamphiI[0], teamDataCount / Point::size(), MPI_DOUBLE,
               &phi[0], teamDataCount / Point::size(), MPI_DOUBLE,
               MASTER, row_comm);
    totalCommTime += commTimer.elapsed();
  }

  double time = timer.elapsed();
  printf("[%d] Timer: %e\n", rank, time);
  printf("[%d] CommTimer: %e\n", rank, totalCommTime);

  if (rank == MASTER) {
    double checkSum = std::accumulate(phi.begin(), phi.end(), 0.0);
    std::cout << "Scatter - checksum answer is: " << checkSum << std::endl;
    std::ofstream phi_file("data/phi.txt");
    phi_file << phi << std::endl;
  }

  MPI_Finalize();
  return 0;
}
