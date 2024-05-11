#include "eigen_solver.h"

void orthogonalize_cols(std::vector<sau::my_vector<double>>& previous_b, sau::my_vector<double>& new_b){
    for(int j = 0; j < previous_b.size(); j++){
        new_b -= previous_b[j] * ((previous_b[j] % new_b) / (previous_b[j] % previous_b[j]));
    }
    new_b.normalize();
}

std::vector<double> davidson(sau::my_sparse_matrix<double>& H, int num_sols){
    double tol = 0.000001;
    int k = H.rows();

    sau::my_matrix<double> K(1, 1); 

    std::vector<sau::my_vector<double>> basis;
    std::vector<sau::my_vector<double>> Hb;

    std::vector<std::vector<double>> eig_vals;
    sau::my_matrix<double> eig_vec;

    sau::my_vector<double> b(H.rows());
    b.fill_random();
    b.normalize();

    basis.push_back(b);

    for(int i = 1; i < num_sols; i++){
        b.fill_random();
        orthogonalize_cols(basis, b);
        basis.push_back(b);
    }

    sau::my_vector<double> temp_Hb;
    int old_n = 0, n = 0;

    int current_num_sols = num_sols;
    std::vector<bool> still_testing;
    for(int i = 0; i < num_sols; i++) still_testing.push_back(true);

    int current_solution = 0;

    for(int i = 0; i < k; i++){
        old_n = n;
        n += current_num_sols;
        K.resize(n, n);

        for(int j = old_n; j < n; j++) Hb.push_back(H * basis[j]);

        for(int j = old_n; j < n; j++){
            for(int m = 0; m < j + 1; m++){
                K(m, j) = K(j, m) = basis[m] % Hb[j];
            }
        }

        std::tuple <sau::my_vector<double>, sau::my_matrix<double>> eigen_results = solve_eigensystem(K);

        std::vector<double> new_eig_value_set;
        for(int j = 0; j < num_sols; j++) new_eig_value_set.push_back(std::get<0>(eigen_results)[j]);
    
        eig_vals.push_back(new_eig_value_set);
        
        eig_vec = std::get<1>(eigen_results);

        std::vector<sau::my_vector<double>> eig_vecs;
        for(int j = 0; j < num_sols; j++) eig_vecs.push_back(eig_vec.extract_column(j));
        
        std::vector<sau::my_vector<double>> new_basis;
        for(int j = 0; j < current_num_sols; j++){
            sau::my_vector<double> temp_basis(H.rows());
            temp_basis.zeros();
            new_basis.push_back(temp_basis);
        }

        int basis_count = 0;
        for(int p = 0; p < num_sols; p++){
            if(still_testing[p]){
                for(int j = 0; j < n; j++){
                    for(int m = 0; m < H.rows(); m++){
                        new_basis[basis_count][m] += eig_vecs[p][j] * (Hb[j][m] - eig_vals.back()[p] * basis[j][m]) / (eig_vals.back()[p] - H.diagonal()[m]);
                    }
                }
                basis_count++;
            }
        }
    
        for(int j = 0; j < current_num_sols; j++){
            if(new_basis[j].distance() == 0){
                new_basis[j].fill_random();
            }
        }
        

        for(int j = 0; j < current_num_sols; j++){
            orthogonalize_cols(basis, new_basis[j]);
            basis.push_back(new_basis[j]);
        }
    
        if(i > 2){
            if(still_testing[current_solution] && std::abs(eig_vals[i][current_solution] - eig_vals[i - 1][current_solution]) < tol){
                still_testing[current_solution] = false;
                current_num_sols--;
                current_solution++;
            }
        }
        if(current_solution >= num_sols){
            std::cout << "Iterations: " << i << std::endl;
            return eig_vals.back();
        }
    }
    std::cout << "Iterations: " << k << std::endl;
    return eig_vals.back();
}

double lanczos(sau::my_sparse_matrix<double>& H){
    double tol = 0.00001;
    sau::my_matrix<double> K(1, 1);
    double eig_prev;
    double eig = MAXFLOAT;
    double m = 0;
    
    sau::my_vector<double> V(H.rows());
    V.fill_random();
    V.normalize();

    sau::my_vector<double> V_prev(V);
    sau::my_vector<double> V_next;

    sau::my_vector<double> W = H * V;
    K(0, 0) = V % W;

    V_next = W - V * K(0,0);

    do{
        V_prev = V; V = V_next; m += 1; eig_prev = eig;
        K.resize(m + 1, m + 1);

        K(m, m - 1) = K(m - 1, m) = V.distance();
        V.normalize();

        W = H * V;
        K(m, m) = V % W;
        
        V_next = W - V * K(m, m) - V_prev * K(m, m - 1);

        std::tuple <sau::my_vector<double>, sau::my_matrix<double>> eigen_results = solve_eigensystem(K);
        eig = std::get<0>(eigen_results)[0];
    } while(std::abs(eig - eig_prev) > tol);

    std::cout << "Lanczos iterations: " << m << std::endl;

    return eig;
}

