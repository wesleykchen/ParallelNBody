#pragma

#include <mpi.h>

namespace meta {

template <typename T>
struct mpi_type;

template <>
struct mpi_type<double> {
  constexpr static MPI_Datatype value = MPI_DOUBLE;
};

template <>
struct mpi_type<float> {
  constexpr static MPI_Datatype value = MPI_FLOAT;
};

} // end namespace meta
