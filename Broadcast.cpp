#include "Util.hpp"

// Broadcast version of n-body algorithm, with 3 points

void calculate(int rank, int numtasks, double data[][DATADIM], double sigma[], double myanswer[]);
int calcStart (int rank, int numtasks);
int calcEnd (int rank, int numtasks);

int main(int argc, char** argv)
{
  int numtasks, rank, count;
  int rc = 0 ;

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
  double data[NUMPOINTS][DATADIM];
  double sigma[NUMPOINTS];

  // hold gathered answers
  double phi[NUMPOINTS];

  Clock timer;
  timer.start();

  Clock commTimer;
  double totalCommTime = 0;
  double tempTime;

  /* Start master tasks */
  if (rank == MASTER) {

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
  }

  // for all processors

  count = DATADIM * NUMPOINTS;

  commTimer.start();
  MPI_Bcast(data, count, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
  //printf("rank= %d  Results: %f %f %f %f\n", rank, data[0][0], data[0][1], data[0][2], data[1][0]);

  MPI_Bcast(sigma, NUMPOINTS, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

  tempTime = commTimer.elapsed();
  totalCommTime += tempTime;

  // all processors have a chunk to hold their temporary answers
  double myphi[NUMPOINTS / numtasks];
  // clear this space
  for (int i = 0; i < NUMPOINTS / numtasks; ++i) {
    myphi[i] = 0;
  }

  // evaluate computation
  calculate(rank, numtasks, data, sigma, myphi);

  //printf("Processor %d computes its answer of %e\n", rank, myphi);

  // Collect results and display

  commTimer.start();
  //MPI_Reduce(&myphi, &phi, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
  MPI_Gather(myphi, NUMPOINTS / numtasks, MPI_DOUBLE, &phi, NUMPOINTS / numtasks, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

  tempTime = commTimer.elapsed();
  totalCommTime += tempTime;

  if (rank == MASTER) {
    double checkSum = 0;
    for ( int i = 0; i < NUMPOINTS; ++i) {
      checkSum += phi[i];
    }
    printf("Broadcast - the checksum is: %e\n", checkSum);
    printPhi(phi);
  }

  double time = timer.elapsed();
  printf("[%d] Timer: %e\n", rank, time);
  printf("[%d] CommTimer: %e\n", rank, totalCommTime);

  MPI_Finalize();
  return 0;
}

void calculate(int rank, int numtasks, double data[][DATADIM], double sigma[NUMPOINTS], double myanswer[])
{
  // calculate start position based on own rank, numranks, NUMPOINTS
  int startIndex = calcStart(rank, numtasks);
  int endIndex = calcEnd(rank, numtasks);

  // printf("Processor %d is starting at %d and going through %d\n", rank, startIndex, endIndex - 1);

  // calculate the n-body portion

  for (int i = startIndex; i < endIndex; ++i) {
    for (int j = 0; j < NUMPOINTS; ++j) {
      myanswer[i - startIndex] += evaluate(data[i][0], data[i][1], data[i][2],
                                           data[j][0], data[j][1], data[j][2]) * sigma[j];
    }
  }
}

int calcStart (int rank, int numtasks) {
  int index;
  if (numtasks - rank >= NUMPOINTS - numtasks * (NUMPOINTS / numtasks)) {
    index = rank * (NUMPOINTS / numtasks);
  }
  else {
    index = rank * (NUMPOINTS / numtasks) - numtasks + rank + NUMPOINTS % numtasks;
  }
  return index;
}

int calcEnd (int rank, int numtasks) {
  return calcStart(rank + 1, numtasks);
}
