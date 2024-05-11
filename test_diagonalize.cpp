#include "eigen_solver.h"

int main(int argc, char* argv[]){
    int n = std::stoi(argv[1]);
    double sparsity = std::stod(argv[2]);
    int num_eig = std::stoi(argv[3]);
    
    sau::my_sparse_matrix<double> mat1(n, n, sparsity);
    mat1.fill_rand_sym();

    std::vector<double> eig;

    // std::cout << "Symmetric Matrix" << std::endl;
    //mat1.print();

    
    sau::my_matrix<double> dense_mat = mat1;

    
    std::tuple <sau::my_vector<double>, sau::my_matrix<double>> eigen_results = solve_eigensystem(dense_mat);

    std::cout << std::endl << "Eigen Value from Lapack" << std::endl;
    for(int i = 0; i < num_eig; i++) std::cout << std::get<0>(eigen_results)[i] << std::endl;
    
    
    // std::get<0>(eigen_results).print();

    // std::cout << std::endl << "Eigen Vectors" << std::endl;
    // std::get<1>(eigen_results).print();

    std::cout << std::endl << "Sending to Davidson" << std::endl;
    eig = davidson(mat1, num_eig);
    std::cout << "Davidson estimate: " << std::endl;
    for(int i = 0; i < num_eig; i++) std::cout << eig[i] << std::endl;

    return 0;
}