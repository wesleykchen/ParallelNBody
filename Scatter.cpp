#include "Util.hpp"

// Scatter version of the n-body algorithm

void calculate(int numtasks, double xI[][DATADIM], double xJ[][DATADIM],
               double sigmaJ[], double myanswer[]);


int main(int argc, char** argv)
{
  int numtasks, rank, count;
  int rc = -1;
  int tag1 = 1;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  if (NUMPOINTS % numtasks != 0) {
    printf("Quitting. The number of processors must divide the total number of tasks.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
    exit(0);
  }
  std::string x, y, z;
  std::ifstream inputFile(PHIDATA);

  // declare data for original as well as chunk
  double xI[NUMPOINTS/numtasks][DATADIM];
  double xJ[NUMPOINTS/numtasks][DATADIM];
  double sigmaJ[NUMPOINTS/numtasks];

  // hold gathered answers
  double phi[NUMPOINTS];

  count = DATADIM * NUMPOINTS / numtasks;

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

        //std::cout << x << "." << y << "." << z << std::endl;

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
      std::cout << "Unable to open sigma file";
    }

    commTimer.start();
    // Scatter data
    MPI_Scatter(data, count, MPI_DOUBLE,
                xI, count, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(sigma, count / DATADIM, MPI_DOUBLE,
                sigmaJ, count / DATADIM, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
  } else {
    commTimer.start();
    // Receive scattered data
    MPI_Scatter(NULL, 0, MPI_DOUBLE,
                xI, count, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    MPI_Scatter(NULL, 0, MPI_DOUBLE,
                sigmaJ, count / DATADIM, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
  }

  tempTime = commTimer.elapsed();
  totalCommTime += tempTime;
  // initialize xJ with xI

  // memcpy(xJ, xI, sizeof(xI) * c);

  for ( int i = 0; i < NUMPOINTS / numtasks; ++i) {
    for ( int j = 0; j < DATADIM; ++j) {
      xJ[i][j] = xI[i][j];
    }
  }

  //printf("Processor %d will now work with data starting this value %e with a xI value of %e and a sigmaJ of %e\n", rank, xJ[0][0], xI[0][0], sigmaJ[0]);

  double myphi[NUMPOINTS / numtasks];
  // clear this space
  for (int i = 0; i < NUMPOINTS / numtasks; ++i) {
    myphi[i] = 0;
  }

  // calculate using first chunk of data
  calculate(numtasks, xI, xJ, sigmaJ, myphi);

  MPI_Status status;
  for (int shiftCount = 1; shiftCount < numtasks; ++shiftCount) {

    commTimer.start();
    int prev  = (rank - 1) % numtasks;
    int next = (rank + 1) % numtasks;

    MPI_Sendrecv_replace(xJ, count, MPI_DOUBLE,
                                 next, tag1, prev, tag1,
                                 MPI_COMM_WORLD, &status);

    MPI_Sendrecv_replace(sigmaJ, count / DATADIM, MPI_DOUBLE,
                                 next, tag1, prev, tag1,
                                 MPI_COMM_WORLD, &status);

    tempTime = commTimer.elapsed();
    totalCommTime += tempTime;

    //printf("After swap number %d, processor %d will now work with data starting this value %e with a xI value of %e with a sigmaJ value of %e\n", shiftCount, rank, xJ[0][0], xI[0][0], sigmaJ[0]);
    calculate(numtasks, xI, xJ, sigmaJ, myphi);
    //myanswer += calculate(xI, xJ, rank, numtasks);
  }

  // printf("Processor %d computes its answer of %e\n", rank, myanswer);

  commTimer.start();
  // Collect results and display
  MPI_Gather(myphi, NUMPOINTS / numtasks, MPI_DOUBLE,
                     phi, NUMPOINTS / numtasks, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
  //MPI_Reduce(&myanswer, &answer, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
  tempTime = commTimer.elapsed();
  totalCommTime += tempTime;

  if (rank == MASTER) {
    double checkSum = 0;
    for ( int i = 0; i < NUMPOINTS; ++i) {
      checkSum += phi[i];
    }
    printf("Scatter - checksum is: %e\n", checkSum);
    printPhi(phi);
  }

  double time = timer.elapsed();
  printf("[%d] Timer: %e\n", rank, time);
  printf("[%d] CommTimer: %e\n", rank, totalCommTime);

  MPI_Finalize();
  return 0;
}

void calculate(int numtasks, double xI[][DATADIM], double xJ[][DATADIM], double sigmaJ[], double myanswer[])
{
  // calculate the n-body portion
  int blockSize = NUMPOINTS / numtasks;
  for (int i = 0; i < blockSize; ++i) {
    for (int j = 0; j < blockSize; ++j) {
      myanswer[i] += evaluate(xI[i][0], xI[i][1], xI[i][2],
                              xJ[j][0], xJ[j][1], xJ[j][2]) * sigmaJ[j];
    }
  }
}
