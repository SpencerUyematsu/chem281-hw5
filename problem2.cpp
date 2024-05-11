#include <bitset>
#include <iostream>
#include <vector>
#include <filesystem>
#include <string>
#include <fstream>

#include "eigen_solver.h"
#include "integrals.h++"


namespace fs = std::filesystem;

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

template <size_t N>
std::vector<bitset<N>> CISD(bitset<N> base){
    std::vector<bitset<N>> possible_configs;
    possible_configs.push_back(base);

    bitset<N> single_excitation;
    // Single excitation
    for(int i = 0; i < 7; i++){
        for(int j = 7; j < 18; j++){
            single_excitation = swap(base, i, j);
            possible_configs.push_back(single_excitation);
            for(int i2 = 18; i2 < 25; i2++){
                for(int j2 = 25; j2 < 36; j2++){
                    possible_configs.push_back(swap(single_excitation, i2, j2));
                }
            }
            for(int i2 = 0; i2 < i; i2++){
                for(int j2 = 7; j2 < j; j2++){
                    possible_configs.push_back(swap(single_excitation, i2, j2));
                }
            }
        }
    }

    for(int i = 18; i < 25; i++){
        for(int j = 25; j < 36; j++){
            single_excitation = swap(base, i, j);
            possible_configs.push_back(single_excitation);
            for(int i2 = 18; i2 < i; i2++){
                for(int j2 = 25; j2 < j; j2++){
                    possible_configs.push_back(swap(single_excitation, i2, j2));
                }
            }
        }
        
    }
    
    return possible_configs;
}

int trailing_zeros(unsigned long long int n){
    if (n == 0) return 64;
    
    int count = 0;
    while(!(n & 1)){
        count++;
        n >>= 1;
    }
    return count;
}

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

    // extract 
    std::vector<bitset<64>> possible_configs = CISD(base_config);
    // for(int i = 0; i < possible_configs.size(); i++) std::cout << possible_configs[i] << std::endl;
    std::cout << possible_configs.size() << std::endl;

    std::cout << possible_configs[0] << std::endl;
    std::cout << possible_configs[1] << std::endl;
    auto result = possible_configs[0] ^ possible_configs[1];
    std::cout << trailing_zeros(possible_configs[1].to_ulong()) << std::endl;

    std::ofstream output("data/results3.csv");

    output << "R,state_0,state_1,state_2,state_3,state_4,state_5,HF" << std::endl;

    auto input_files = get_inputs();

    for(int i = 0; i < input_files.size(); i++){
        std::cout << "CALCULATING VALUES FOR: " << input_files[i] << std::endl;
        integrals integral (18,input_files[i]);

        sau::my_sparse_matrix<double> H(possible_configs.size(), possible_configs.size());

        for(int i = 0; i < possible_configs.size(); i++){
            for(int j = 0; j < possible_configs.size(); j++){
                auto XOR = possible_configs[i] ^ possible_configs[j];
                int pop_count = XOR.count();
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
        std::cout << "Sending to Davidson" << std::endl;
        auto eig = davidson(H, 6);

        output << input_files[i].substr(17, 4) << ",";
        for(int i = 0; i < eig.size(); i++) output << eig[i] + integral.get_core_energy() << ",";
        output << integral.calc_zero_move_matrix_element(base_config) + integral.get_core_energy() << std::endl;
    }
}