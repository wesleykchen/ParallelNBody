
// Symmetric Team Scatter version of the n-body algorithm

#include <tuple>

#include "Util.hpp"

#include "kernel/InvSq.kern"
#include "meta/kernel_traits.hpp"
#include "meta/random.hpp"

typedef std::tuple<int,int>      xy_pair;
typedef std::tuple<int,int,int>  itc_tuple;
typedef std::tuple<int,int>      ir_pair;

struct IndexTransformer {
  IndexTransformer(int num_teams, int team_size)
      : T(num_teams), C(team_size) {
  }

  /** Convert from iteration/team/member (i,t,c) to 2D block index (X,Y) */
  xy_pair itc2xy(int i, int t, int c) const {
    int X = t;
    int Y = (t + c + i * C) % T;
    return xy_pair(X, Y);
  }

  /** Transpose in XY */
  xy_pair xyTranspose(xy_pair xy) const {
    return xy_pair(std::get<1>(xy),std::get<0>(xy));
  }

  /** Convert from 2D block index (X,Y) to iteration/team/member (i,t,c) */
  itc_tuple xy2itc(int X, int Y) const {
    int t = X;
    int intermediate = (Y - X + T) % T;
    // fix from zero
    int c = (intermediate + C) % C;
    int i = intermediate / C;
    return itc_tuple(i,t,c);
  }

  /** Convert a TC pair to a process rank */
  int tc2r(int t, int c) const {
    return t * C + c;
  }

  /** Convert ITC tuple to IR pair */
  ir_pair itc2ir(const itc_tuple& itc) const {
    return ir_pair(std::get<0>(itc), tc2r(std::get<1>(itc), std::get<2>(itc)));
  }

  xy_pair ir2xy(int i, int r) const {
    int t = r / C;
    int c = r % C;
    return itc2xy(i,t,c);
  }

  xy_pair ir2xy(const ir_pair& ir) const {
    return ir2xy(std::get<0>(ir), std::get<1>(ir));
  }

 private:
  int T;   // The number of process teams in the computation
  int C;   // The size of the process teams in the computation
};

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

  if (arg.size() < 2) {
    std::cerr << "Usage: " << arg[0] << " NUMPOINTS [-c TEAMSIZE] [-nocheck]" << std::endl;
    exit(1);
  }

  srand(time(NULL));
  unsigned N = string_to_<int>(arg[1]);

  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int P;
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  // Scratch request for MPI
  MPI_Request request;
  // Scratch status for MPI
  MPI_Status status;

  typedef InvSq kernel_type;
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

  if (rank == MASTER) {
    // generate source data
    for (unsigned i = 0; i < N; ++i)
      source.push_back(meta::random<source_type>::get());

    // generate charge data
    for (unsigned i = 0; i < N; ++i)
      charge.push_back(meta::random<charge_type>::get());

    // display metadata
    std::cout << "N = " << N << std::endl;
    std::cout << "P = " << P << std::endl;
    std::cout << "Teamsize = " << teamsize << std::endl;
  }

  Clock timer;
  Clock compTimer;
  Clock splitTimer;
  Clock reduceTimer;
  Clock shiftTimer;
  Clock sendRecvTimer;

  double totalCompTime = 0;
  double totalSplitTime = 0;
  double totalReduceTime = 0;
  double totalShiftTime = 0;
  double totalSendRecvTime = 0;

  timer.start();

  // Broadcast the size of the problem to all processes
  splitTimer.start();
  MPI_Bcast(&N, sizeof(N), MPI_CHAR, MASTER, MPI_COMM_WORLD);
  totalSplitTime += splitTimer.elapsed();

  // Broadcast the teamsize to all processes
  splitTimer.start();
  MPI_Bcast(&teamsize, sizeof(teamsize), MPI_CHAR, MASTER, MPI_COMM_WORLD);
  totalSplitTime += splitTimer.elapsed();

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

  /***********/
  /** SETUP **/
  /***********/

  unsigned num_teams = P / teamsize;
  // Determine coordinates in processor team grid
  unsigned team  = rank / teamsize;
  unsigned trank = rank % teamsize;

  // Split comm into row and column communicators
  MPI_Comm team_comm;
  MPI_Comm_split(MPI_COMM_WORLD, team, rank, &team_comm);
  MPI_Comm row_comm;
  MPI_Comm_split(MPI_COMM_WORLD, trank, rank, &row_comm);

  // Create transformer to help us convert between block idx and proc index
  IndexTransformer transformer(num_teams, teamsize);

  /********************/
  /** INITIALIZATION **/
  /********************/

  // Declare data for the block computations
  std::vector<source_type> xJ(idiv_up(N,num_teams));
  std::vector<charge_type> cJ(idiv_up(N,num_teams));
  std::vector<result_type> rJ(idiv_up(N,num_teams));

  // Scatter data from master to team leaders
  if (trank == MASTER) {
    splitTimer.start();
    MPI_Scatter(source.data(), sizeof(source_type) * xJ.size(), MPI_CHAR,
                xJ.data(), sizeof(source_type) * xJ.size(), MPI_CHAR,
                MASTER, row_comm);
    MPI_Scatter(charge.data(), sizeof(charge_type) * cJ.size(), MPI_CHAR,
                cJ.data(), sizeof(charge_type) * cJ.size(), MPI_CHAR,
                MASTER, row_comm);
    totalSplitTime += splitTimer.elapsed();
  }

  // Team leaders broadcast to team
  splitTimer.start();
  MPI_Bcast(xJ.data(), sizeof(source_type) * xJ.size(), MPI_CHAR,
            MASTER, team_comm);
  MPI_Bcast(cJ.data(), sizeof(charge_type) * cJ.size(), MPI_CHAR,
            MASTER, team_comm);
  totalSplitTime += splitTimer.elapsed();

  // Copy xJ -> xI
  std::vector<source_type> xI = xJ;
  // Copy cJ -> cI
  std::vector<charge_type> cI = cJ;
  // Initialize block result rI
  std::vector<result_type> rI(idiv_up(N,num_teams));
  // Declare space for receiving
  std::vector<result_type> temp_rI(idiv_up(N,num_teams));

  // Perform initial offset by teamrank
  shiftTimer.start();
  // Add num_teams to prevent negative numbers
  int src = (team + trank + num_teams) % num_teams;
  int dst = (team - trank + num_teams) % num_teams;
  MPI_Sendrecv_replace(xJ.data(), sizeof(source_type) * xJ.size(), MPI_CHAR,
                       dst, 0, src, 0,
                       row_comm, &status);
  MPI_Sendrecv_replace(cJ.data(), sizeof(charge_type) * cJ.size(), MPI_CHAR,
                       dst, 0, src, 0,
                       row_comm, &status);
  totalShiftTime += shiftTimer.elapsed();

  /**********************/
  /** ZEROTH ITERATION **/
  /**********************/

  int curr_iter = 0;

  int last_iter = idiv_up(num_teams + 1, 2*teamsize);

  if (trank == MASTER) {
    // first iteration masters are computing the symmetric diagonal
    compTimer.start();
    p2p(K, xJ.begin(), xJ.end(), cJ.begin(), rI.begin());
    totalCompTime += compTimer.elapsed();

  } else {

    // Get current block's xy-index
    xy_pair xy = transformer.itc2xy(curr_iter, team, trank);

    // Convert the xy to the transpose
    xy_pair xy_trans = transformer.xyTranspose(xy);

    itc_tuple itc_trans = transformer.xy2itc(std::get<0>(xy_trans), std::get<1>(xy_trans));    

    // calculate needs of sending to transpose
    if (std::get<0>(itc_trans) > curr_iter && std::get<0>(itc_trans) != last_iter - 1) {

      // Set rJ to zero
      rJ.assign(rJ.size(), result_type());

      // Compute symmetric off-diagonal
      compTimer.start();
      p2p(K,
          xJ.begin(), xJ.end(), cJ.begin(), rJ.begin(),
          xI.begin(), xI.end(), cI.begin(), rI.begin());
      totalCompTime += compTimer.elapsed();

      int send_dest = transformer.tc2r(std::get<1>(itc_trans),std::get<2>(itc_trans)); 

      //std::cout << "Processor: " << rank << "sending to: " << send_dest << std::endl;
      // Send to proper rank
      sendRecvTimer.start();
        MPI_Isend(rJ.data(), sizeof(result_type) * rJ.size(), MPI_CHAR,
                send_dest, 0,
                MPI_COMM_WORLD, &request);
      totalSendRecvTime += sendRecvTimer.elapsed();
    }
    else {
      // Set rJ to zero
      rJ.assign(rJ.size(), result_type());

      // Compute asymmetric off-diagonal
      compTimer.start();
      p2p(K,
          xJ.begin(), xJ.end(), cJ.begin(),
          xI.begin(), xI.end(), rI.begin());
      totalCompTime += compTimer.elapsed();
    }
  }

  //MPI_Barrier(MPI_COMM_WORLD); 
  if (trank != MASTER) {

    xy_pair recv_dest_xy = transformer.xyTranspose(transformer.ir2xy(xy_pair(num_teams/teamsize - curr_iter - 1, rank)));

    itc_tuple recv_dest_itc = transformer.xy2itc(std::get<0>(recv_dest_xy),std::get<1>(recv_dest_xy));
    int recv_dest = transformer.tc2r(std::get<1>(recv_dest_itc), std::get<2>(recv_dest_itc));

    if (num_teams / teamsize - curr_iter - 1 != last_iter - 1) {

      //std::cout << "Processor: " << rank << "receiving from: " << recv_dest << std::endl;

      // Receive
      sendRecvTimer.start();
      MPI_Recv(temp_rI.data(), sizeof(result_type) * temp_rI.size(), MPI_CHAR,
               recv_dest, 0,
               MPI_COMM_WORLD, &status);
      totalSendRecvTime += sendRecvTimer.elapsed();

      // Accumulate to current answer
      for (auto r = rI.begin(), tr = temp_rI.begin(); r != rI.end(); ++r, ++tr)
        *r += *tr;
    }
  }
  /********************/
  /** ALL ITERATIONS **/
  /********************/

  for (++curr_iter; curr_iter < last_iter; ++curr_iter) {

    // Shift data to the next process to compute the next block
    shiftTimer.start();
    int src = (team + teamsize + num_teams) % num_teams;
    int dst = (team - teamsize + num_teams) % num_teams;
    MPI_Sendrecv_replace(xJ.data(), sizeof(source_type) * xJ.size(), MPI_CHAR,
                         dst, 0, src, 0,
                         row_comm, &status);
    MPI_Sendrecv_replace(cJ.data(), sizeof(charge_type) * cJ.size(), MPI_CHAR,
                         dst, 0, src, 0,
                         row_comm, &status);
    totalShiftTime += shiftTimer.elapsed();

    // Get current block's xy-index
    xy_pair xy = transformer.itc2xy(curr_iter, team, trank);

    // Convert the xy to the transpose
    xy_pair xy_trans = transformer.xyTranspose(xy);

    itc_tuple itc_trans = transformer.xy2itc(std::get<0>(xy_trans), std::get<1>(xy_trans));    

    //std::cout << "On iteration: " << curr_iter << std::endl;

    // calculate needs of sending to transpose and computing
    if (std::get<0>(itc_trans) > curr_iter && std::get<0>(itc_trans) != last_iter - 1) {

      // Set rJ to zero
      rJ.assign(rJ.size(), result_type());

      // Compute symmetric off-diagonal
      compTimer.start();
      p2p(K,
          xJ.begin(), xJ.end(), cJ.begin(), rJ.begin(),
          xI.begin(), xI.end(), cI.begin(), rI.begin());
      totalCompTime += compTimer.elapsed();

      int send_dest = transformer.tc2r(std::get<1>(itc_trans),std::get<2>(itc_trans));

      //std::cout << "Processor: " << rank << "sending to: " << send_dest << std::endl;

      // Send to proper rank
      sendRecvTimer.start();
        MPI_Isend(rJ.data(), sizeof(result_type) * rJ.size(), MPI_CHAR,
                send_dest, 0,
                MPI_COMM_WORLD, &request);
      totalSendRecvTime += sendRecvTimer.elapsed();
    } else {
      // Set rJ to zero
      rJ.assign(rJ.size(), result_type());

      // Compute asymmetric off-diagonal
      compTimer.start();
      p2p(K,
          xJ.begin(), xJ.end(), cJ.begin(),
          xI.begin(), xI.end(), rI.begin());
      totalCompTime += compTimer.elapsed();
     
      if (curr_iter == last_iter -1 ) {
        break;
      }

    }


    int iPrimeOffset = ((trank == 0) ? 0 : 1);
    xy_pair recv_dest_xy = transformer.xyTranspose(transformer.ir2xy(xy_pair(num_teams/teamsize - curr_iter - iPrimeOffset, rank)));

    itc_tuple recv_dest_itc = transformer.xy2itc(std::get<0>(recv_dest_xy),std::get<1>(recv_dest_xy));
    int recv_dest = transformer.tc2r(std::get<1>(recv_dest_itc), std::get<2>(recv_dest_itc));

    if (num_teams/teamsize - curr_iter - iPrimeOffset != last_iter - 1) {

      //std::cout << "Processor: " << rank << "receiving from: " << recv_dest << std::endl;
      // Receive
      sendRecvTimer.start();
      MPI_Recv(temp_rI.data(), sizeof(result_type) * temp_rI.size(), MPI_CHAR,
               recv_dest, 0,
               MPI_COMM_WORLD, &status);
      totalSendRecvTime += sendRecvTimer.elapsed();

      // Accumulate to current answer
      for (auto r = rI.begin(), tr = temp_rI.begin(); r != rI.end(); ++r, ++tr)
        *r += *tr;
    }
      
  }  //  end for iteration

  // Reduce answers to the team leader
  reduceTimer.start();

  // TODO: Generalize
  static_assert(std::is_same<result_type, double>::value,
                "Need result_type == double for now");
  MPI_Reduce(rI.data(), temp_rI.data(), rI.size(), MPI_DOUBLE,
             MPI_SUM, MASTER, team_comm);
  totalReduceTime += reduceTimer.elapsed();

  // Allocate result on master
  std::vector<result_type> result;
  if (rank == MASTER) {
    result = std::vector<result_type>(P*idiv_up(N,P));
  }

  // Gather team leader answers to master
  if (trank == MASTER) {
    reduceTimer.start();
    MPI_Gather(temp_rI.data(), sizeof(result_type) * temp_rI.size(), MPI_CHAR,
               result.data(), sizeof(result_type) * temp_rI.size(), MPI_CHAR,
               MASTER, row_comm);
    totalReduceTime += reduceTimer.elapsed();
  }

  double time = timer.elapsed();
  printf("[%d] Timer: %e\n", rank, time);
  printf("[%d] CommTimer: %e\n", rank, totalSendRecvTime + totalReduceTime
                                       + totalSplitTime + totalShiftTime);
  printf("[%d] CompTimer: %e\n", rank, totalCompTime);

  // Check the result
  if (rank == MASTER && checkErrors) {
    std::cout << "Computing direct matvec..." << std::endl;

    std::vector<result_type> exact(N);

    // Compute the result with a direct matrix-vector multiplication
    compTimer.start();
    p2p(K, source.begin(), source.end(), charge.begin(), exact.begin());
    double directCompTime = compTimer.elapsed();

    print_error(exact, result);
    std::cout << "DirectCompTime: " << directCompTime << std::endl;

    std::ofstream exact_file("data/exact.txt");
    exact_file << exact << std::endl;
  }

  if (rank == MASTER) {
    std::ofstream result_file("data/result.txt");
    result_file << result << std::endl;
  }

  MPI_Finalize();
  return 0;
}
