#include <tuple>
#include <deque>

#include "Util.hpp"

#include "kernel/NonParaBayesian.kern"
#include "meta/kernel_traits.hpp"


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


// Symmetric Team Scatter version of the n-body algorithm

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
  // Scratch request for MPI
  MPI_Request request;
  // Scratch status for MPI
  MPI_Status status;

  typedef NonParaBayesian kernel_type;
  kernel_type K(1,1);

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
      // exit(1);
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



  // TODO: Remove, just for proof of concept and testing

  // Compute the rank and iteration of the our symmetric blocks
  std::vector<ir_pair> iter_rank_trans;

  for (unsigned iter = 0; iter < idiv_up(num_teams, teamsize); ++iter) {
    // Get our future block coordinates
    xy_pair xy = transformer.itc2xy(iter, team, trank);
    // Get the ITC of the process responsible for the symmtric block
    itc_tuple itc = transformer.xy2itc(std::get<1>(xy), std::get<0>(xy));
    // Convert ITC to IR
    iter_rank_trans.push_back(transformer.itc2ir(itc));
  }

  // Copy and sort by iteration
  std::deque<ir_pair> iter_rank_deque(iter_rank_trans.begin(), iter_rank_trans.end());
  // Just for now... make sure this is sorted by iteration
  // TODO: Avoid this sort by computing them in the right order?
  std::sort(iter_rank_deque.begin(), iter_rank_deque.end());


  /********************/
  /** INITIALIZATION **/
  /********************/

  // Declare data for the block computations
  std::vector<source_type> xJ(idiv_up(N,num_teams));
  std::vector<charge_type> cJ(idiv_up(N,num_teams));
  std::vector<result_type> rJ(idiv_up(N,num_teams));

  // Scatter data from master to team leaders
  if (trank == MASTER) {
    commTimer.start();
    MPI_Scatter(source.data(), sizeof(source_type) * xJ.size(), MPI_CHAR,
                xJ.data(), sizeof(source_type) * xJ.size(), MPI_CHAR,
                MASTER, row_comm);
    MPI_Scatter(charge.data(), sizeof(charge_type) * cJ.size(), MPI_CHAR,
                cJ.data(), sizeof(charge_type) * cJ.size(), MPI_CHAR,
                MASTER, row_comm);
    totalCommTime += commTimer.elapsed();
  }

  // Team leaders broadcast to team
  commTimer.start();
  MPI_Bcast(xJ.data(), sizeof(source_type) * xJ.size(), MPI_CHAR,
            MASTER, team_comm);
  MPI_Bcast(cJ.data(), sizeof(charge_type) * cJ.size(), MPI_CHAR,
            MASTER, team_comm);
  totalCommTime += commTimer.elapsed();

  // Copy xJ -> xI
  std::vector<source_type> xI = xJ;
  // Copy cJ -> cI
  std::vector<charge_type> cI = cJ;
  // Initialize block result rI
  std::vector<result_type> rI(idiv_up(N,num_teams));
  // Declare space for receiving
  std::vector<result_type> temp_rI(idiv_up(N,num_teams));

  // Perform initial offset by teamrank
  commTimer.start();
  // Add num_teams to prevent negative numbers
  int src = (team + trank + num_teams) % num_teams;
  int dst = (team - trank + num_teams) % num_teams;
  MPI_Sendrecv_replace(xJ.data(), sizeof(source_type) * xJ.size(), MPI_CHAR,
                       dst, 0, src, 0,
                       row_comm, &status);
  MPI_Sendrecv_replace(cJ.data(), sizeof(charge_type) * cJ.size(), MPI_CHAR,
                       dst, 0, src, 0,
                       row_comm, &status);
  totalCommTime += commTimer.elapsed();

  /**********************/
  /** ZEROTH ITERATION **/
  /**********************/

  int curr_iter = 0;

  if (trank == MASTER) {
    // If this is the team leader, compute the symmetric diagonal block
    p2p(K, xJ.begin(), xJ.end(), cJ.begin(), rI.begin());

    // The front of the list should be this process
    assert(iter_rank_deque.front() == ir_pair(curr_iter, rank));
    iter_rank_deque.pop_front();
  } else {
    // Compute and send if this IR < transpose IR
    if (std::tie(curr_iter, rank) < iter_rank_trans[curr_iter]) {
      // Set rJ to zero
      rJ.assign(rJ.size(), result_type());

      // Compute
      p2p(K,
          xJ.begin(), xJ.end(), cJ.begin(), rJ.begin(),
          xI.begin(), xI.end(), cI.begin(), rI.begin());

      // Send
      commTimer.start();
      MPI_Isend(rJ.data(), sizeof(result_type) * rJ.size(), MPI_CHAR,
                std::get<1>(iter_rank_trans[curr_iter]), 0,
                MPI_COMM_WORLD, &request);
      totalCommTime += commTimer.elapsed();

      // Find and remove from receive list  YUCK HACK
      // Get our block coordinates
      xy_pair xy = transformer.itc2xy(curr_iter, team, trank);
      // Get the ITC of the process responsible for the symmtric block
      itc_tuple itc = transformer.xy2itc(std::get<1>(xy), std::get<0>(xy));
      // Convert ITC to IR
      ir_pair ir = transformer.itc2ir(itc);
      // Remove -- YUCK!
      auto it = std::find(iter_rank_deque.begin(), iter_rank_deque.end(), ir);
      assert(it != iter_rank_deque.end());
      iter_rank_deque.erase(it);
      assert(std::is_sorted(iter_rank_deque.begin(), iter_rank_deque.end()));
    }
  }

  // Receive all symmetric blocks computed at this iteration
  while (!iter_rank_deque.empty()
         && std::get<0>(iter_rank_deque.front()) == curr_iter) {
    // Recv
    MPI_Recv(temp_rI.data(), sizeof(result_type) * temp_rI.size(), MPI_CHAR,
             std::get<1>(iter_rank_deque.front()), 0,
             MPI_COMM_WORLD, &status);
    // Accumulate
    for (auto r = rI.begin(), tr = temp_rI.begin(); r != rI.end(); ++r, ++tr)
      *r += *tr;
    // Pop
    iter_rank_deque.pop_front();
  }

  /********************/
  /** ALL ITERATIONS **/
  /********************/

  int last_iter = idiv_up(num_teams + 1, 2*teamsize);

  for (++curr_iter; curr_iter < last_iter; ++curr_iter) {
    MPI_Barrier(MPI_COMM_WORLD);  // To make sure it's not an rJ race?

    // Shift data to the next process to compute the next block
    commTimer.start();
    int src = (team + teamsize + num_teams) % num_teams;
    int dst = (team - teamsize + num_teams) % num_teams;
    MPI_Sendrecv_replace(xJ.data(), sizeof(source_type) * xJ.size(), MPI_CHAR,
                         dst, 0, src, 0,
                         row_comm, &status);
    MPI_Sendrecv_replace(cJ.data(), sizeof(charge_type) * cJ.size(), MPI_CHAR,
                         dst, 0, src, 0,
                         row_comm, &status);
    totalCommTime += commTimer.elapsed();

    // Compute and send if this IR < transpose IR
    assert(std::tie(curr_iter, rank) != iter_rank_trans[curr_iter]);
    if (std::tie(curr_iter, rank) < iter_rank_trans[curr_iter]) {
      // Set rJ to zero
      rJ.assign(rJ.size(), result_type());

      // Compute
      p2p(K,
          xJ.begin(), xJ.end(), cJ.begin(), rJ.begin(),
          xI.begin(), xI.end(), cI.begin(), rI.begin());

      // Send
      commTimer.start();
      MPI_Isend(rJ.data(), sizeof(result_type) * rJ.size(), MPI_CHAR,
                std::get<1>(iter_rank_trans[curr_iter]), 0,
                MPI_COMM_WORLD, &request);
      totalCommTime += commTimer.elapsed();

      // Find and remove from receive list  YUCK HACK
      // Get our block coordinates
      xy_pair xy = transformer.itc2xy(curr_iter, team, trank);
      // Get the ITC of the process responsible for the symmtric block
      itc_tuple itc = transformer.xy2itc(std::get<1>(xy), std::get<0>(xy));
      // Convert ITC to IR
      ir_pair ir = transformer.itc2ir(itc);
      // Remove -- YUCK!
      auto it = std::find(iter_rank_deque.begin(), iter_rank_deque.end(), ir);
      assert(it != iter_rank_deque.end());
      iter_rank_deque.erase(it);
      assert(std::is_sorted(iter_rank_deque.begin(), iter_rank_deque.end()));
    }

    // Receive all symmetric blocks computed at this iteration
    while (!iter_rank_deque.empty()
           && std::get<0>(iter_rank_deque.front()) == curr_iter) {
      // Recv
      MPI_Recv(temp_rI.data(), sizeof(result_type) * temp_rI.size(), MPI_CHAR,
               std::get<1>(iter_rank_deque.front()), 0,
               MPI_COMM_WORLD, &status);
      // Accumulate
      for (auto r = rI.begin(), tr = temp_rI.begin(); r != rI.end(); ++r, ++tr)
        *r += *tr;
      // Pop
      iter_rank_deque.pop_front();
    }
  }  //  end for iteration

  // Make sure all the blocks are accounted for
  assert(iter_rank_deque.empty());

  // Reduce answers to the team leader
  commTimer.start();
  // TODO: Generalize
  static_assert(std::is_same<result_type, double>::value,
                "Need result_type == double for now");
  MPI_Reduce(rI.data(), temp_rI.data(), rI.size(), MPI_DOUBLE,
             MPI_SUM, MASTER, team_comm);
  totalCommTime += commTimer.elapsed();

  // Allocate result on master
  std::vector<result_type> result;
  if (rank == MASTER)
    result = std::vector<result_type>(P*idiv_up(N,P));

  // Gather team leader answers to master
  if (trank == MASTER) {
    commTimer.start();
    MPI_Gather(temp_rI.data(), sizeof(result_type) * temp_rI.size(), MPI_CHAR,
               result.data(), sizeof(result_type) * temp_rI.size(), MPI_CHAR,
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
