#include "Util.hpp"

#include "kernel/Laplace.kern"
#include "meta/kernel_traits.hpp"

// Functions to convert between team/member/iteration (T,C,I) and 2D black indexes (X,Y)
void tci2XY (int t, int c, int i, int T, int C, int &X, int &Y) {
    X = t;
    Y = (t + c + i * C) % T;
}

void XY2tci (int X, int Y, int T, int C, int &t, int &c, int &i) {
    t = X;
    int intermediate = (Y-X) % T;
    // fix from zero
    c = (intermediate + C)% C;
    i = intermediate / C;
}

// Symmetric Team Scatter version of the n-body algorithm

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int P;
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  unsigned teamsize;

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
  std::vector<charge_type> sigmaJ(idiv_up(N,num_teams));

  // New for symmetric, must store sigmaI as well
  std::vector<charge_type> sigmaI(idiv_up(N,num_teams));

  // Scatter data from master to team leaders
  if (trank == MASTER) {
    commTimer.start();
    MPI_Scatter(source.data(), sizeof(source_type) * xI.size(), MPI_CHAR,
                xI.data(), sizeof(source_type) * xI.size(), MPI_CHAR,
                MASTER, row_comm);
    MPI_Scatter(charge.data(), sizeof(charge_type) * sigmaJ.size(), MPI_CHAR,
                sigmaJ.data(), sizeof(charge_type) * sigmaJ.size(), MPI_CHAR,
                MASTER, row_comm);
    totalCommTime += commTimer.elapsed();
  }

  // Team leaders broadcast to team
  commTimer.start();
  MPI_Bcast(xI.data(), sizeof(source_type) * xI.size(), MPI_CHAR,
            MASTER, team_comm);
  MPI_Bcast(sigmaJ.data(), sizeof(charge_type) * sigmaJ.size(), MPI_CHAR,
            MASTER, team_comm);
  totalCommTime += commTimer.elapsed();

  // Copy xI -> xJ
  std::copy(xI.begin(), xI.end(), xJ.begin());

  // Copy sigmaJ -> sigmaI - new for symmetric
  std::copy(sigmaJ.begin(), sigmaJ.end(), sigmaI.begin());

  // Perform initial offset by teamrank
  commTimer.start();
  // Add num_teams to prevent negative numbers
  int prev = (team - trank + num_teams) % num_teams;
  int next = (team + trank + num_teams) % num_teams;
  MPI_Status status;
  MPI_Sendrecv_replace(xJ.data(), sizeof(source_type) * xJ.size(),
                       MPI_CHAR, next, 0, prev, 0,
                       row_comm, &status);
  MPI_Sendrecv_replace(sigmaJ.data(), sizeof(charge_type) * sigmaJ.size(),
                       MPI_CHAR, next, 0, prev, 0,
                       row_comm, &status);
  totalCommTime += commTimer.elapsed();

  // TODO Should this be result type?

  // Hold accumulated answers
  std::vector<result_type> phiI(idiv_up(N,num_teams));

  // Answers for symmetric part - new
  std::vector<result_type> phiJ(idiv_up(N,num_teams));

  // Calculate the current block, new for symm using symP2P off diagonal
  p2p(K,
      xJ.begin(), xJ.end(), sigmaJ.begin(), phiI.begin(),
      xI.begin(), xI.end(), sigmaI.begin(), phiJ.begin());

  // Start adding changes to team scatter

  // Looping process to shift the data between the teams
  // Now fewer times then before
  int totalIterations = idiv_up(num_teams + 1, 2 * teamsize);;

  // Declare space for receiving
  std::vector<result_type> receivedPhiI(idiv_up(N,num_teams));

  // Calculate where to send to and if it needs to send
  int myX. myY = 0;
  tci2XY(team, trank, iteration, num_teams, teamsize, myX, myY);

  int mirrorX = myY;
  int mirrorY = myX;
  int mirrorT, mirrorC, mirrorI;
  XY2tci(mirrorX, mirrorY, num_teams, teamsize, mirrorT, mirrorC, mirrorI);

  if ( trank != 0) {

    // Send/receive data
    // phiJ.size() == xJ.size()
    commTimer.start();
    MPI_Sendrecv(phiJ.data(), sizeof(result_type) * phiJ.size(),
                MPI_CHAR, teamsize*mirrorT + mirrorC, 0,
                receivedPhiI.data(), sizeof(result_type) * phiJ.size(),
                MPI_CHAR, MPI_ANY_SOURCE, 0,
                MPI_COMM_WORLD, &status);
    totalCommTime += commTimer.elapsed();

    // Accumulate into phiI
    for (auto iout = phiI.begin(), iin = receivedPhiI.begin(); iout != phiI.end(); ++iout, ++iin)
      *iout += *iin;
  }

  std::cout<< "Processor: "<< rank << "Team: " << team<< "Trank: " << trank<< "Numteams: " << num_teams<< "Teamsize: " << teamsize<< "I: " << 0 << " " << myX<< " " << myY << std::endl;

  // ENTER LOOP
  for (int iteration = 1; iteration < totalIterations; ++iteration) {
    commTimer.start();

    int prev = (team - teamsize) % num_teams;
    int next = (team + teamsize) % num_teams;
    MPI_Sendrecv_replace(xJ.data(), sizeof(source_type) * xJ.size(),
                         MPI_CHAR, next, 0, prev, 0,
                         row_comm, &status);
    MPI_Sendrecv_replace(sigmaJ.data(), sizeof(charge_type) * sigmaJ.size(),
                         MPI_CHAR, next, 0, prev, 0,
                         row_comm, &status);
    totalCommTime += commTimer.elapsed();

    // Allocate and phiJ space and set to zero here
    // TODO

    // Compute block
    p2p(K,
        xJ.begin(), xJ.end(), sigmaJ.begin(), phiI.begin(),
        xI.begin(), xI.end(), sigmaI.begin(), phiJ.begin());

    // Calculate where to send to and if it needs to send
    tci2XY(team, trank, iteration, num_teams, teamsize, myX, myY);

    mirrorX = myY;
    mirrorY = myX;
    XY2tci(mirrorX, mirrorY, num_teams, teamsize, mirrorT, mirrorC, mirrorI);

    std::cout<< "Processor: "<< rank << "Team: " << team<< "Trank: " << trank<< "Numteams: " << num_teams<< "Teamsize: " << teamsize<< "I: " << iteration<< " " << myX<< " " << myY << std::endl;

    if ( mirrorI != totalIterations - 1) {

      // Send/receive data
      // phiJ.size() == xJ.size()
      commTimer.start();
      MPI_Sendrecv(phiJ.data(), sizeof(result_type) * phiJ.size(),
                   MPI_CHAR, teamsize*mirrorT + mirrorC, 0,
                   receivedPhiI.data(), sizeof(result_type) * phiJ.size(),
                   MPI_CHAR, MPI_ANY_SOURCE, 0,
                   MPI_COMM_WORLD, &status);
      totalCommTime += commTimer.elapsed();

      // Accumulate into phiI
      for (auto iout = phiI.begin(), iin = receivedPhiI.begin(); iout != phiI.end(); ++iout, ++iin)
        *iout += *iin;
    }
  }

  // Allocate teamphiI on team leaders
  std::vector<result_type> teamphiI;
  if (trank == MASTER)
    teamphiI = std::vector<result_type>(idiv_up(N,num_teams));

  // Reduce answers to the team leader
  commTimer.start();
  // TODO: Generalize
  static_assert(std::is_same<result_type, double>::value,
                "Need result_type == double for now");
  MPI_Reduce(phiI.data(), teamphiI.data(), phiI.size(),
             MPI_DOUBLE, MPI_SUM, MASTER, team_comm);
  totalCommTime += commTimer.elapsed();

  // Allocate phi on master
  std::vector<result_type> phi;
  if (rank == MASTER)
    phi = std::vector<result_type>(P*idiv_up(N,P));

  // Gather team leader answers to master
  if (trank == MASTER) {
    commTimer.start();
    MPI_Gather(teamphiI.data(), sizeof(result_type) * teamphiI.size(), MPI_CHAR,
               phi.data(), sizeof(result_type) * teamphiI.size(), MPI_CHAR,
               MASTER, row_comm);
    totalCommTime += commTimer.elapsed();
  }

  double time = timer.elapsed();
  printf("[%d] Timer: %e\n", rank, time);
  printf("[%d] CommTimer: %e\n", rank, totalCommTime);

  if (rank == MASTER) {
    result_type check = std::accumulate(phi.begin(), phi.end(), result_type());
    std::cout << "Symmetric - checksum answer is: " << check << std::endl;
    std::ofstream phi_file("data/phi.txt");
    phi_file << phi << std::endl;
  }

  MPI_Finalize();
  return 0;
}
