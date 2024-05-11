#include <random>
#include <tuple>
#include <lapack.hh>
#include <vector>
#include <cmath>
#include "data_types.h"

template <typename T>
inline std::tuple <sau::my_vector<T>, sau::my_matrix<T>> solve_eigensystem (sau::my_matrix<T> A) {
    int n = A.rows();
    // const int n=aseq (A.rows (), A.cols (), "solve_eigensystem (mat)" ); // throw if mat is not symmetric
    sau::my_vector<double> E (n);
    const int64_t info = lapack::syev (lapack::Job::Vec, lapack::Uplo::Upper,
                                       n, A.access_raw_pointer (),
                                       n, E.access_raw_pointer ()
                                       );
    if (info) throw "LAPACK exited with code " + std::to_string (info);
    return {E,A};
}

double lanczos(sau::my_sparse_matrix<double>& H);
std::vector<double> davidson(sau::my_sparse_matrix<double>& H, int k);