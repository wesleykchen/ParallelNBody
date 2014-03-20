#include "Util.hpp"

#include "kernel/Laplace.kern"
#include "meta/kernel_traits.hpp"

// Team Scatter version of the n-body algorithm

int main(int argc, char** argv)
{
  bool checkErrors = true;
  unsigned teamsize = 1;

  // Parse optional command line args
  std::vector<std::string> arg(argv, argv + argc);
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
    if (arg[i] == "-nocheck") {
      checkErrors = false;
      arg.erase(arg.begin() + i, arg.begin() + i + 1);  // Erase this arg
      --i;                                              // Reset index
    }
  }

  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int P;
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  // Scratch status for MPI
  MPI_Status status;

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
    if (arg.size() < 3) {
      std::cerr << "Usage: " << arg[0] << " SOURCE_FILE CHARGE_FILE [-c TEAMSIZE]" << std::endl;
      //exit(1);
      // XXX: Remove
      std::cerr << "Using default " << SOURCE_DATA << " " << CHARGE_DATA << std::endl;

      arg.resize(1);  // keep only the executable name
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
    std::cout << "Teamsize = " << teamsize << std::endl;
  }

  Clock timer;
  Clock commTimer;
  double totalCommTime = 0;

  timer.start();

  // Broadcast the size of the problem to all processes
  commTimer.start();
  MPI_Bcast(&N, sizeof(N), MPI_CHAR, MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  // Broadcast the teamsize to all processes
  commTimer.start();
  MPI_Bcast(&teamsize, sizeof(teamsize), MPI_CHAR, MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

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

  unsigned num_teams = P / teamsize;
  // Determine coordinates in processor team grid
  unsigned team  = rank / teamsize;
  unsigned trank = rank % teamsize;

  // Split comm into row and column communicators
  MPI_Comm team_comm;
  MPI_Comm_split(MPI_COMM_WORLD, team, rank, &team_comm);
  MPI_Comm row_comm;
  MPI_Comm_split(MPI_COMM_WORLD, trank, rank, &row_comm);

  // Declare data for the block computations
  std::vector<source_type> xI(idiv_up(N,num_teams));
  std::vector<source_type> xJ(idiv_up(N,num_teams));
  std::vector<charge_type> cJ(idiv_up(N,num_teams));

  // Scatter data from master to team leaders
  if (trank == MASTER) {
    commTimer.start();
    MPI_Scatter(source.data(), sizeof(source_type) * xI.size(), MPI_CHAR,
                xI.data(), sizeof(source_type) * xI.size(), MPI_CHAR,
                MASTER, row_comm);
    MPI_Scatter(charge.data(), sizeof(charge_type) * cJ.size(), MPI_CHAR,
                cJ.data(), sizeof(charge_type) * cJ.size(), MPI_CHAR,
                MASTER, row_comm);
    totalCommTime += commTimer.elapsed();
  }

  // Team leaders broadcast to team
  commTimer.start();
  MPI_Bcast(xI.data(), sizeof(source_type) * xI.size(), MPI_CHAR,
            MASTER, team_comm);
  MPI_Bcast(cJ.data(), sizeof(charge_type) * cJ.size(), MPI_CHAR,
            MASTER, team_comm);
  totalCommTime += commTimer.elapsed();

  // Copy xI -> xJ
  std::copy(xI.begin(), xI.end(), xJ.begin());

  // Perform initial offset by teamrank
  commTimer.start();

  int prev = (team - trank + num_teams) % num_teams;
  int next = (team + trank) % num_teams;
  MPI_Sendrecv_replace(xJ.data(), sizeof(source_type) * xJ.size(),
                       MPI_CHAR, next, 0, prev, 0,
                       row_comm, &status);
  MPI_Sendrecv_replace(cJ.data(), sizeof(charge_type) * cJ.size(),
                       MPI_CHAR, next, 0, prev, 0,
                       row_comm, &status);
  totalCommTime += commTimer.elapsed();

  // Hold accumulated answers
  std::vector<double> rI(idiv_up(N,num_teams));

  if (trank == MASTER) {
    // If this is the team leader, compute the symmetric diagonal block
    p2p(K,
        xJ.begin(), xJ.end(), cJ.begin(), rI.begin());
  } else {
    // Else, compute the off-diagonal block
    p2p(K,
        xJ.begin(), xJ.end(), cJ.begin(),
        xI.begin(), xI.end(), rI.begin());
  }

  // Looping process to shift the data between the teams
  int ceilPC2 = idiv_up(P, teamsize*teamsize);
  for (int shiftCount = 1; shiftCount < ceilPC2; ++shiftCount) {
    commTimer.start();
    // Add num_teams to prevent negative numbers
    int prev = (team - teamsize + num_teams) % num_teams;
    int next = (team + teamsize + num_teams) % num_teams;
    MPI_Sendrecv_replace(xJ.data(), sizeof(source_type) * xJ.size(),
                         MPI_CHAR, next, 0, prev, 0,
                         row_comm, &status);
    MPI_Sendrecv_replace(cJ.data(), sizeof(charge_type) * cJ.size(),
                         MPI_CHAR, next, 0, prev, 0,
                         row_comm, &status);
    totalCommTime += commTimer.elapsed();

    // Compute on the last iteration only if
    // 1) The teamsize divides the number of teams (everyone computes)
    // 2) Your team rank is one of the remainders
    if (shiftCount < ceilPC2-1
        || (num_teams % teamsize == 0 || trank < num_teams % teamsize)) {
      p2p(K,
          xJ.begin(), xJ.end(), cJ.begin(),
          xI.begin(), xI.end(), rI.begin());
    }
  }

  // Allocate teamrI on team leaders
  std::vector<result_type> teamrI;
  if (trank == MASTER)
    teamrI = std::vector<result_type>(idiv_up(N,num_teams));

  // Reduce answers to the team leader
  commTimer.start();
  // TODO: Generalize
  static_assert(std::is_same<result_type, double>::value,
                "Need result_type == double for now");
  MPI_Reduce(rI.data(), teamrI.data(), rI.size(),
             MPI_DOUBLE, MPI_SUM, MASTER, team_comm);
  totalCommTime += commTimer.elapsed();

  // Allocate result on master
  std::vector<result_type> result;
  if (rank == MASTER)
    result = std::vector<result_type>(P*idiv_up(N,P));

  // Gather team leader answers to master
  if (trank == MASTER) {
    commTimer.start();
    MPI_Gather(teamrI.data(), sizeof(result_type) * teamrI.size(), MPI_CHAR,
               result.data(), sizeof(result_type) * teamrI.size(), MPI_CHAR,
               MASTER, row_comm);
    totalCommTime += commTimer.elapsed();
  }

  double time = timer.elapsed();
  printf("[%d] Timer: %e\n", rank, time);
  printf("[%d] CommTimer: %e\n", rank, totalCommTime);

  // Check the result
  if (rank == MASTER && checkErrors) {
    std::cout << "Computing direct matvec..." << std::endl;

    std::vector<result_type> exact(N);

    // Compute the result with a direct matrix-vector multiplication
    p2p(K, source.begin(), source.end(), charge.begin(), exact.begin());

    print_error(exact, result);
  }

  if (rank == MASTER) {
    std::ofstream result_file("data/result.txt");
    result_file << result << std::endl;
  }

  MPI_Finalize();
  return 0;
}
