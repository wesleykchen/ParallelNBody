#include "Util.hpp"
#define TEAMSIZE 1

// Scatter version of the n-body algorithm

void calculate(int numtasks, double xI[][DATADIM], double xJ[][DATADIM],
                 double sigmaJ[], double myanswer[]);

int main(int argc, char** argv)
{
  int numtasks, rank, size;
  int rc = -1;
  int tag1 = 1;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);


  int Nteams = numtasks / TEAMSIZE;
  // int c = TEAMSIZE;

  if (NUMPOINTS % numtasks != 0) {
    printf("Quitting. The number of processors must divide the total number of points.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
    exit(0);
  }

  if (numtasks % TEAMSIZE != 0) {
    printf("Quitting. The teamsize (c) must divide the total number of processors (p).\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
    exit(0);
  }

  if (TEAMSIZE * TEAMSIZE > numtasks) {
    printf("Quitting. The teamsize ^ 2 (c^2) must be less than or equal to the number of  processors (p).\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
    exit(0);
  }

  // Determine coordinates in processor team grid
  int team = rank / TEAMSIZE;
  int trank = rank % TEAMSIZE;


  // Split comm into row and column communicators

  MPI_Comm team_comm;
  MPI_Comm_split(MPI_COMM_WORLD, team, trank, &team_comm);

  MPI_Comm row_comm;
  MPI_Comm_split(MPI_COMM_WORLD, trank, team, &row_comm);

  // declarations: TODO: does everyone needs this or just master
  string x, y, z;
  ifstream inputFile(PHIDATA);

  // declare data for original as well as chunk

  double xI[NUMPOINTS/Nteams][DATADIM];
  double xJ[NUMPOINTS/Nteams][DATADIM];
  double sigmaJ[NUMPOINTS/Nteams];

  // not all processors need this memory declared...

  // hold team gathered phi
  double teamphi[NUMPOINTS/Nteams];

  // hold gathered answers
  double phi[NUMPOINTS];

  int teamDataCount = DATADIM * NUMPOINTS / Nteams;

  Clock timer;
  timer.start();

  Clock commTimer;
  double totalCommTime = 0;
  double tempTime;

  /* Start master tasks */
  if (rank == MASTER) {
    double data[NUMPOINTS][DATADIM];
    double sigma[NUMPOINTS];

    if (inputFile.is_open()) {
      int index = 0;
      while ( true ) {

        getline( inputFile, x, '\t');
        if (inputFile.eof()) {
          break;
        }
        getline( inputFile, y, '\t');
        getline( inputFile, z, '\n');

        //cout << x << "." << y << "." << z << endl;

        data[index][0] = atof(x.c_str());
        data[index][1] = atof(y.c_str());
        data[index][2] = atof(z.c_str());

        ++index;
      }
      inputFile.close();
      inputFile.clear();
    }
    else {
      printf("Unable to open file");
    }

    inputFile.open(SIGMADATA);
    if (inputFile.is_open()) {
      int index = 0;
      while ( true ) {
        getline( inputFile, x, '\n');
        if (inputFile.eof()) {
          break;
        }
        sigma[index] = atof(x.c_str());

        ++index;
      }
      inputFile.close();
    } else {
      cout << "Unable to open sigma file";
    }
    commTimer.start();
    // Master scatter to team leaders
    MPI_Scatter(data, teamDataCount, MPI_DOUBLE,
                xI, teamDataCount, MPI_DOUBLE, MASTER, row_comm);
    MPI_Scatter(sigma, teamDataCount / DATADIM, MPI_DOUBLE,
                sigmaJ, teamDataCount / DATADIM, MPI_DOUBLE, MASTER, row_comm);
    // printf("Processor %d of team %d  received first xI of %e\n", rank, team, xI[0][0]);
    tempTime = commTimer.elapsed();
    totalCommTime += tempTime;
  } else if (trank == MASTER){
    commTimer.start();
    // Receive scattered data
    MPI_Scatter(NULL, 0, MPI_DOUBLE,
                xI, teamDataCount, MPI_DOUBLE, MASTER, row_comm);

    MPI_Scatter(NULL, 0, MPI_DOUBLE,
                sigmaJ, teamDataCount / DATADIM, MPI_DOUBLE, MASTER, row_comm);
    // printf("Processor %d of team %d  received first xI of %e\n", rank, team, xI[0][0]);
    tempTime = commTimer.elapsed();
    totalCommTime += tempTime;
  } else {
    // printf("Processor %d of team %d is happy\n", rank, team);
    // do nothing
  }

  commTimer.start();
  // Team leaders broadcast to team
  MPI_Bcast(xI, teamDataCount, MPI_DOUBLE, MASTER, team_comm);
  MPI_Bcast(sigmaJ, teamDataCount / DATADIM, MPI_DOUBLE, MASTER, team_comm);

  tempTime = commTimer.elapsed();
  totalCommTime += tempTime;

  // printf("At first, team %d processor %d has value of %e\n", team, trank, xI[0][0]);

  // initialize xJ with xI

  memcpy(xJ, xI, sizeof(xI));

  // printf("Processor %d will now work with data starting this value %e\n with a xI value of %e and a sigmaJ of %e\n", rank, xJ[0][0], xI[0][0], sigmaJ[0]);

  MPI_Status status;

  commTimer.start();
  // perform initial offset by teamrank
  // % is implementation defined, adding Nteams to prevent negative numbers
  int prev = (team - trank + Nteams) % Nteams;
  int next = (team + trank +Nteams) % Nteams;
  MPI_Sendrecv_replace(xJ, teamDataCount, MPI_DOUBLE,
                         next, tag1, prev, tag1,
                         row_comm, &status);
  MPI_Sendrecv_replace(sigmaJ, teamDataCount / DATADIM, MPI_DOUBLE,
                         next, tag1, prev, tag1,
                         row_comm, &status);

  tempTime = commTimer.elapsed();
  totalCommTime += tempTime;
  // printf("After offset, team %d processor %d starting with xJ of %e\n", team, trank, xJ[0][0]);

  double myphi[NUMPOINTS / Nteams];
  // clear this space
  for (int i = 0; i < NUMPOINTS / Nteams; ++i) {
    myphi[i] = 0;
  }


  // printf("After shift #0, team %d processor %d xI start: %e\n", team, trank, xJ[0][0]);
  // calculate using first chunk of data
  calculate(Nteams, xI, xJ, sigmaJ, myphi);

  // move into the looping process to shift the data between the teams

  int ceilPC2 = (numtasks + TEAMSIZE * TEAMSIZE - 1) / TEAMSIZE / TEAMSIZE;
  //int ceilPC2 = 2;
  for (int shiftCount = 1; shiftCount < ceilPC2; ++shiftCount) {
    commTimer.start();
    // % is implementation defined, adding Nteams to prevent negative numbers
    int prev = (team - TEAMSIZE + Nteams) % Nteams;
    int next = (team + TEAMSIZE + Nteams) % Nteams;
    MPI_Sendrecv_replace(xJ, teamDataCount, MPI_DOUBLE,
                       next, tag1, prev, tag1,
                       row_comm, &status);
    MPI_Sendrecv_replace(sigmaJ, teamDataCount / DATADIM, MPI_DOUBLE,
                       next, tag1, prev, tag1,
                       row_comm, &status);
    tempTime = commTimer.elapsed();
    totalCommTime += tempTime;
    // printf("After shift #%d, team %d processor %d xI start: %e, prev: %d, next: %d\n", shiftCount, team, trank, xJ[0][0], prev, next);

    // test if last iteration to do calculation exceptions
    if ( shiftCount + 1 == ceilPC2) {
      // calculate leftovers to calculate
      if (trank <  Nteams % TEAMSIZE) {
        calculate(Nteams, xI, xJ, sigmaJ, myphi);
      }
      // acount for perfect division
      else if (Nteams % TEAMSIZE == 0) {
        calculate(Nteams, xI, xJ, sigmaJ, myphi);
      }
      else{
        // do nothing
      }
    }
    else {
      calculate(Nteams, xI, xJ, sigmaJ, myphi);
    }
    //myanswer += calculate(xI, xJ, rank, numtasks);
  }

  commTimer.start();
  // Sum reduce answers to the team leader
  MPI_Reduce(myphi, teamphi, teamDataCount / DATADIM, MPI_DOUBLE, MPI_SUM, MASTER, team_comm);

  tempTime = commTimer.elapsed();
  totalCommTime += tempTime;


  // Master gathers team leader results (team Master so to speak)
  if (trank == MASTER) {
    commTimer.start();
    MPI_Gather(teamphi, teamDataCount / DATADIM, MPI_DOUBLE,
             phi, teamDataCount / DATADIM, MPI_DOUBLE, MASTER, row_comm);
    tempTime = commTimer.elapsed();
    totalCommTime += tempTime;
  }

  if (rank == MASTER) {
    double checkSum = 0;
    for ( int i = 0; i < NUMPOINTS; ++i) {
      checkSum += phi[i];
    }
    printf("TeamScatter - the checksum is: %e\n", checkSum);
    printPhi(phi);
  }

  double time = timer.elapsed();
  printf("[%d] Timer: %e\n", rank, time);
  printf("[%d] CommTimer: %e\n", rank, totalCommTime);

  MPI_Finalize();
  return 0;
}

void calculate(int Nteams, double xI[][DATADIM], double xJ[][DATADIM], double sigmaJ[], double myanswer[])
{
  // calculate the n-body portion
  int blockSize = NUMPOINTS / Nteams;
  for (int i = 0; i < blockSize; ++i) {
    for (int j = 0; j < blockSize; ++j) {
      myanswer[i] += evaluate(xI[i][0], xI[i][1], xI[i][2],
                         xJ[j][0], xJ[j][1], xJ[j][2]) * sigmaJ[j];
    }
  }
}
