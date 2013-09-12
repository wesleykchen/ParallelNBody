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
          // Erase these two elements
          arg.erase(arg.begin() + i, arg.begin() + i + 2);
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

  // Broadcast the size of the problem to all processes
  timer.start();
  commTimer.start();
  MPI_Bcast(&N, 1, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  // Broadcast the teamsize to all processes
  timer.start();
  commTimer.start();
  MPI_Bcast(&teamsize, 1, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  // TODO: How to generalize?
  if (N % P != 0) {
    printf("Quitting. The number of processors must divide the total number of tasks.\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    exit(0);
  }

  if (N % teamsize != 0) {
    printf("Quitting. The teamsize (c) must divide the total number of processors (p).\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    exit(0);
  }

  if (teamsize * teamsize > unsigned(P)) {
    printf("Quitting. The teamsize ^ 2 (c^2) must be less than or equal to the number of  processors (p).\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    exit(0);
  }


  int num_teams = P / teamsize;


  // Determine coordinates in processor team grid
  unsigned team  = rank / teamsize;
  unsigned trank = rank % teamsize;


  // Split comm into row and column communicators

  MPI_Comm team_comm;
  MPI_Comm_split(MPI_COMM_WORLD, team, trank, &team_comm);

  MPI_Comm row_comm;
  MPI_Comm_split(MPI_COMM_WORLD, trank, team, &row_comm);

  // declarations: TODO: does everyone needs this or just master
  std::string x, y, z;
  std::ifstream inputFile(PHIDATA);

  // declare data for original as well as chunk

  std::vector<Point> xI(idiv_up(N,num_teams));
  std::vector<Point> xJ(idiv_up(N,num_teams));
  std::vector<double> sigmaJ(idiv_up(N,num_teams));

  /*
    double xI[NUMPOINTS/num_teams][DATADIM];
    double xJ[NUMPOINTS/num_teams][DATADIM];
    double sigmaJ[NUMPOINTS/num_teams];
  */

  // not all processors need this memory declared - switch to only team leaders later
  // hold team gathered phi
  std::vector<double> teamphiI(idiv_up(N,num_teams));
  //double teamphi[NUMPOINTS/num_teams];

  // hold gathered answers
  std::vector<double> phiI(idiv_up(N,P));
  //double phi[NUMPOINTS];

  unsigned teamDataCount = Point::size() * N / num_teams;
  // int teamDataCount = DATADIM * NUMPOINTS / num_teams;

  /* Start master tasks */
  if (rank == MASTER) {

    commTimer.start();
    // Master scatter to team leaders
    MPI_Scatter(&data[0], teamDataCount, MPI_DOUBLE,
                &xI[0], teamDataCount, MPI_DOUBLE,
                MASTER, row_comm);
    MPI_Scatter(&sigma[0], teamDataCount / Point::size(), MPI_DOUBLE,
                &sigmaJ[0], teamDataCount / Point::size(), MPI_DOUBLE,
                MASTER, row_comm);

    // printf("Processor %d of team %d  received first xI of %e\n", rank, team, xI[0][0]);
    totalCommTime += commTimer.elapsed();

  } else if (trank == MASTER){

    commTimer.start();
    // Receive scattered data
    MPI_Scatter(NULL, 0, MPI_DOUBLE,
                &xI[0], teamDataCount, MPI_DOUBLE,
                MASTER, row_comm);

    MPI_Scatter(NULL, 0, MPI_DOUBLE,
                &sigmaJ[0], teamDataCount / Point::size(), MPI_DOUBLE,
                MASTER, row_comm);

    // printf("Processor %d of team %d  received first xI of %e\n", rank, team, xI[0][0]);
    totalCommTime += commTimer.elapsed();
  } else {
    // printf("Processor %d of team %d is happy\n", rank, team);
    // do nothing
  }

  commTimer.start();
  // Team leaders broadcast to team
  MPI_Bcast(&xI[0], teamDataCount, MPI_DOUBLE, MASTER, team_comm);
  MPI_Bcast(&sigmaJ[0], teamDataCount / Point::size(), MPI_DOUBLE, MASTER, team_comm);

  totalCommTime += commTimer.elapsed();
  // printf("At first, team %d processor %d has value of %e\n", team, trank, xI[0][0]);

  // Copy xI -> xJ
  std::copy(xI.begin(), xI.end(), xJ.begin());

  // printf("Processor %d will now work with data starting this value %e\n with a xI value of %e and a sigmaJ of %e\n", rank, xJ[0][0], xI[0][0], sigmaJ[0]);

  MPI_Status status;

  commTimer.start();
  // perform initial offset by teamrank
  // % is implementation defined, adding num_teams to prevent negative numbers
  int prev = (team - trank + num_teams) % num_teams;
  int next = (team + trank + num_teams) % num_teams;
  MPI_Sendrecv_replace(&xJ[0], teamDataCount, MPI_DOUBLE,
                       next, 0, prev, 0,
                       row_comm, &status);
  MPI_Sendrecv_replace(&sigmaJ[0], teamDataCount / Point::size(), MPI_DOUBLE,
                       next, 0, prev, 0,
                       row_comm, &status);

  totalCommTime += commTimer.elapsed();
  // printf("After offset, team %d processor %d starting with xJ of %e\n", team, trank, xJ[0][0]);

  // printf("After shift #0, team %d processor %d xI start: %e\n", team, trank, xJ[0][0]);

  // Calculate the current block
  block_eval(xJ.begin(), xJ.end(), sigmaJ.begin(),
             xI.begin(), xI.end(), phiI.begin());

  // move into the looping process to shift the data between the teams

  int ceilPC2 = (P + teamsize * teamsize - 1) / teamsize / teamsize;

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
    // printf("After shift #%d, team %d processor %d xI start: %e, prev: %d, next: %d\n", shiftCount, team, trank, xJ[0][0], prev, next);

    // test if last iteration to do calculation exceptions
    if ( shiftCount + 1 == ceilPC2) {
      // calculate leftovers to calculate
      if (trank <  num_teams % teamsize) {
        // Calculate the current block
        block_eval(xJ.begin(), xJ.end(), sigmaJ.begin(),
                   xI.begin(), xI.end(), phiI.begin());
      }
      // acount for perfect division
      else if (num_teams % teamsize == 0) {
        // Calculate the current block
        block_eval(xJ.begin(), xJ.end(), sigmaJ.begin(),
                   xI.begin(), xI.end(), phiI.begin());
      }
      else{
        // do nothing
      }
    }
    else {
      // Calculate the current block
      block_eval(xJ.begin(), xJ.end(), sigmaJ.begin(),
                 xI.begin(), xI.end(), phiI.begin());
    }
    //myanswer += calculate(xI, xJ, rank, P);
  }


  std::vector<double> phi;
  if (rank == MASTER)
    phi = std::vector<double>(P*idiv_up(N,P));

  commTimer.start();
  // Sum reduce answers to the team leader
  MPI_Reduce(&phiI[0], &teamphiI[0], teamDataCount / Point::size(), MPI_DOUBLE, MPI_SUM, MASTER, team_comm);

  totalCommTime += commTimer.elapsed();

  // Master gathers team leader results (team Master so to speak)
  if (trank == MASTER) {
    commTimer.start();
    MPI_Gather(&teamphiI[0], teamDataCount / Point::size(), MPI_DOUBLE,
               &phi[0], teamDataCount / Point::size(), MPI_DOUBLE, MASTER, row_comm);
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
