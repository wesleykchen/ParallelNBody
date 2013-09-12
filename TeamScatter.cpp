#include "Util.hpp"
#include "Vec.hpp"

//XXX: Remove
#define TEAMSIZE 4

// Team Scatter version of the n-body algorithm

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int numtasks;
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  if (argc < 3) {
    if (rank == MASTER) {
      std::cerr << "Usage: " << argv[0] << " PHI_FILE SIGMA_FILE TEAMSIZE" << std::endl;
      //exit(1);
      // XXX: Remove
      std::cerr << "Using default " << PHIDATA << " " << SIGMADATA << " " << TEAMSIZE << std::endl;
    }
    // did I edit this block properly? need to allow for TEAMSIZE
    int teamsize = atoi(argv[3]);

    argc = 4;
    char** new_argv = new char*[3];
    new_argv[0] = argv[0];
    new_argv[1] = (char*)& PHIDATA;
    new_argv[2] = (char*)& SIGMADATA;
    //specifically here
    new_argv[3] = atoi(argv[3]);
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
    std::cout << "P = " << numtasks << std::endl;
    std::cout << "Teamsize = " << teamsize << std::endl;
  }

  // Broadcast the size of the problem to all processes
  timer.start();
  commTimer.start();
  MPI_Bcast(&N, 1, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD);
  totalCommTime += commTimer.elapsed();

  // TODO: How to generalize?
  if (N % numtasks != 0) {
    printf("Quitting. The number of processors must divide the total number of tasks.\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    exit(0);
  }

  if (N % teamsize != 0) {
    printf("Quitting. The teamsize (c) must divide the total number of processors (p).\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    exit(0);
  }

  if (teamsize * teamsize > numtasks) {
    printf("Quitting. The teamsize ^ 2 (c^2) must be less than or equal to the number of  processors (p).\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    exit(0);
  }


  int Nteams = numtasks / teamsize;


  // Determine coordinates in processor team grid
  int team = rank / teamsize;
  int trank = rank % teamsize;


  // Split comm into row and column communicators

  MPI_Comm team_comm;
  MPI_Comm_split(MPI_COMM_WORLD, team, trank, &team_comm);

  MPI_Comm row_comm;
  MPI_Comm_split(MPI_COMM_WORLD, trank, team, &row_comm);

  // declarations: TODO: does everyone needs this or just master
  std::string x, y, z;
  std::ifstream inputFile(PHIDATA);

  // declare data for original as well as chunk

  std::vector<Point> xI(N,Nteams);
  std::vector<Point> xJ(N,Nteams));
std::vector<double> sigmaJ(N,Nteams);

/*
  double xI[NUMPOINTS/Nteams][DATADIM];
  double xJ[NUMPOINTS/Nteams][DATADIM];
  double sigmaJ[NUMPOINTS/Nteams];
*/

// not all processors need this memory declared - switch to only team leaders later
// hold team gathered phi
std::vector<double> teamPhiI(N/Nteams));
//double teamphi[NUMPOINTS/Nteams];

// hold gathered answers
std::vector<double> phiI(idiv_up(N,numtasks));
//double phi[NUMPOINTS];

unsigned teamDataCount = Point::size() * N / Nteams;
// int teamDataCount = DATADIM * NUMPOINTS / Nteams;

/* Start master tasks */
if (rank == MASTER) {

  commTimer.start();
  // Master scatter to team leaders
  MPI_Scatter(&data[0], teamDataCount, MPI_DOUBLE,
              &xI[0], teamDataCount, MPI_DOUBLE,
              MASTER, row_comm);
  MPI_Scatter(sigma, teamDataCount / Point::size(), MPI_DOUBLE,
              sigmaJ, teamDataCount / Point::size(), MPI_DOUBLE,
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
// % is implementation defined, adding Nteams to prevent negative numbers
int prev = (team - trank + Nteams) % Nteams;
int next = (team + trank + Nteams) % Nteams;
MPI_Sendrecv_replace(&xJ[0], teamDataCount, MPI_DOUBLE,
                     next, tag1, prev, tag1,
                     row_comm, &status);
MPI_Sendrecv_replace(&sigmaJ[0], teamDataCount / Point::size(), MPI_DOUBLE,
                     next, tag1, prev, tag1,
                     row_comm, &status);

totalCommTime += commTimer.elapsed();
// printf("After offset, team %d processor %d starting with xJ of %e\n", team, trank, xJ[0][0]);

// printf("After shift #0, team %d processor %d xI start: %e\n", team, trank, xJ[0][0]);

// Calculate the current block
block_eval(xJ.begin(), xJ.end(), sigmaJ.begin(),
           xI.begin(), xI.end(), phiI.begin());

// move into the looping process to shift the data between the teams

int ceilPC2 = (numtasks + teamsize * teamsize - 1) / teamsize / teamsize;

for (int shiftCount = 1; shiftCount < ceilPC2; ++shiftCount) {
  commTimer.start();
  // % is implementation defined, adding Nteams to prevent negative numbers
  int prev = (team - TEAMSIZE + Nteams) % Nteams;
  int next = (team + TEAMSIZE + Nteams) % Nteams;
  MPI_Sendrecv_replace(&xJ[0], teamDataCount, MPI_DOUBLE,
                       next, tag1, prev, tag1,
                       row_comm, &status);
  MPI_Sendrecv_replace(&sigmaJ[0], teamDataCount / Point::size(), MPI_DOUBLE,
                       next, tag1, prev, tag1,
                       row_comm, &status);

  totalCommTime += commTimer.elapsed();
  // printf("After shift #%d, team %d processor %d xI start: %e, prev: %d, next: %d\n", shiftCount, team, trank, xJ[0][0], prev, next);

  // test if last iteration to do calculation exceptions
  if ( shiftCount + 1 == ceilPC2) {
    // calculate leftovers to calculate
    if (trank <  Nteams % TEAMSIZE) {
      // Calculate the current block
      block_eval(xJ.begin(), xJ.end(), sigmaJ.begin(),
                 xI.begin(), xI.end(), phiI.begin());
    }
    // acount for perfect division
    else if (Nteams % TEAMSIZE == 0) {
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
  //myanswer += calculate(xI, xJ, rank, numtasks);
}


std::vector<double> phi;
if (rank == MASTER)
  phi = std::vector<double>(numtasks*idiv_up(N,numtasks));

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
