#include <bitset>
#include <iostream>
#include <vector>
#include <filesystem>
#include <string>
#include <fstream>

#include "eigen_solver.h"
#include "integrals.h++"

namespace fs = std::filesystem;

// swaps the bits at two positions. Used for constructing CISD configurations.
template <size_t N>
bitset<N> swap(bitset<N> original, int a, int b) 
{
    if(original[a] && !original[b]){
        original.reset(a);
        original.set(b);
        return original;
    } 
    if(!original[a] && original[b]){
        original.set(a);
        original.reset(b);
        return original;
    } 
    return original;
}

// produce all CISD configurations for N2
template <size_t N>
std::vector<bitset<N>> CISD(bitset<N> base){
    std::vector<bitset<N>> possible_configs;
    possible_configs.push_back(base);

    bitset<N> single_excitation;

    for(int i = 0; i < 7; i++){
        for(int j = 7; j < 18; j++){
            // single excitation within alpha orbitals
            single_excitation = swap(base, i, j);
            possible_configs.push_back(single_excitation);

            // excitation within beta orbitals
            for(int i2 = 18; i2 < 25; i2++){
                for(int j2 = 25; j2 < 36; j2++){
                    possible_configs.push_back(swap(single_excitation, i2, j2));
                }
            }
            // second excitation within alpha orbitals
            for(int i2 = 0; i2 < i; i2++){
                for(int j2 = 7; j2 < j; j2++){
                    possible_configs.push_back(swap(single_excitation, i2, j2));
                }
            }
        }
    }

    for(int i = 18; i < 25; i++){
        for(int j = 25; j < 36; j++){
            // single excitation within beta orbitals
            single_excitation = swap(base, i, j);
            possible_configs.push_back(single_excitation);

            // second excitation within beta orbitals
            for(int i2 = 18; i2 < i; i2++){
                for(int j2 = 25; j2 < j; j2++){
                    possible_configs.push_back(swap(single_excitation, i2, j2));
                }
            }
        }
        
    }
    return possible_configs;
}

// determines trailing zeros in an integers binary notation
int trailing_zeros(unsigned long long int n){
    if (n == 0) return 64;
    
    int count = 0;
    while(!(n & 1)){
        count++;
        n >>= 1;
    }
    return count;
}

// extracts filenames from FCIdumps folder
std::vector<std::string> get_inputs(){
    fs::path folderPath = "FCIdumps";
    
    std::vector<std::string> files;
    for (const auto& entry : fs::directory_iterator(folderPath)) {
        // Check if the current entry is a regular file
        if (entry.is_regular_file()) {
            // Print the path of the file
            files.push_back(entry.path().string());
        }
    }
    return files;
}

// executes calculations and protocol for problem 2
int main(){
    bitset<64> base_config; // (1) electronic configuration

    // Fill orbitals for base configuration of N2
    // Alpha spin
    for(int i = 0; i < 7; i++){
        base_config.set(i);
    }
    // beta spin
    for(int i = 18; i < 25; i++){
        base_config.set(i);
    }

    std::vector<bitset<64>> possible_configs = CISD(base_config); // (2) get all electronic configurations

    std::cout << "Total configurations: " << possible_configs.size() << std::endl;

    // filestream for final results
    std::ofstream output("data/results3.csv");

    output << "R,state_0,state_1,state_2,state_3,state_4,state_5,HF" << std::endl;

    auto input_files = get_inputs();

    // loop through files
    for(int i = 0; i < input_files.size(); i++){
        std::cout << std::endl << "CALCULATING VALUES FOR: " << input_files[i] << std::endl;
        integrals integral (18,input_files[i]); // initiate integral class with current file

        // (3) form sparse matrix representing second quantized electronic hamiltonian
        sau::my_sparse_matrix<double> H(possible_configs.size(), possible_configs.size());
        
        for(int i = 0; i < possible_configs.size(); i++){
            for(int j = 0; j < possible_configs.size(); j++){
                // (3) (a) identify non zero matrix elements and 
                auto XOR = possible_configs[i] ^ possible_configs[j];
                int pop_count = XOR.count();

                // calculate Hamiltonian value based on the number of bit changes
                if(pop_count == 0) H.insert(integral.calc_zero_move_matrix_element(possible_configs[j]), i, j);

                if(pop_count == 2){

                    int a = trailing_zeros((XOR & possible_configs[i]).to_ullong());
                    int i_value = trailing_zeros((XOR & possible_configs[j]).to_ullong());
                    H.insert(integral.calc_one_move_matrix_element(possible_configs[j], i_value, a), i, j);
                }
            
                if(pop_count == 4){
                    int b = trailing_zeros((XOR & possible_configs[i]).to_ullong());
                    int a = trailing_zeros((XOR & possible_configs[i]).reset(b).to_ullong());

                    int i_value = trailing_zeros((XOR & possible_configs[j]).to_ullong());
                    int j_value = trailing_zeros((XOR & possible_configs[j]).reset(i_value).to_ullong());

                    bitset<64> this_j = possible_configs[j];
                    H.insert(integral.calc_two_move_matrix_element(this_j, j_value, i_value, a, b), i, j);
                }
            }
        }

        // (4) calculate eigen values using Davidson
        std::cout << "Sending to Davidson" << std::endl;
        auto eig = davidson(H, 6);

        // (5) calculate total energy and output results
        output << input_files[i].substr(17, 4) << ",";
        for(int i = 0; i < eig.size(); i++) output << eig[i] + integral.get_core_energy() << ",";
        output << integral.calc_zero_move_matrix_element(base_config) + integral.get_core_energy() << std::endl;
    }
}
