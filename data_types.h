#include <valarray>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <omp.h>

#include "utilities.h"

namespace sau {
template <typename T>
class my_vector: public std::valarray <T> {
    private:
        T max_random;
        T min_random;
    
    public:
        my_vector(): std::valarray<T>() {set_max_min();}
        my_vector(int size): std::valarray<T>(size) { zeros(); set_max_min();}
        my_vector(std::valarray<T>& other): std::valarray<T>(other) {set_max_min();}

        void resize(size_t new_size){
            my_vector<T> temp(new_size);

            size_t copy_size = std::min(new_size, this->size());

            for(size_t i = 0; i < copy_size; i++){
                temp[i] = (*this)[i];
            }

            this->swap(temp);
        }

        void fill(T value){
            for(int i = 0; i < this->size(); i++){
                (*this)[i] = value;
            }
        }

        void fill_random(){
            for(int i = 0; i < this->size(); i++){
                (*this)[i] = getRandomValue(min_random, max_random);
            }
        }

        double distance(){
            double distance = 0;

            for(int i = 0; i < this->size(); i++){
                distance += (*this)[i] * (*this)[i];
            }

            return std::sqrt(distance);
        }

        void normalize(){
            double this_distance = distance();

            for(int i = 0; i < this->size(); i++){
                (*this)[i] /= this_distance;
            }
        }

        void zeros(){ fill(0); }

        const void print() const{
            for(int i = 0; i < this->size(); i++){
                std::cout << (*this)[i] << "\t";
            }
            std::cout << std::endl;
        }

        T* access_raw_pointer(){
            return &(*this)[0];
        }

        T operator%(const my_vector<T>& other) const {
            T total = 0;
            std::valarray<T> product = (*this) * other;
            for(int i = 0; i < this->size(); i++){
                total += product[i];
            }
            return total;
        }

        my_vector<T> operator*(const double scalar) const {
            my_vector<T> result(*this); 

            for (size_t i = 0; i < this->size(); ++i) {
                result[i] *= scalar;
            }

            return result;
        }

        my_vector<T> operator+(const double scalar) const {
            my_vector<T> result(*this); 

            for (size_t i = 0; i < this->size(); ++i) {
                result[i] += scalar;
            }

            return result;
        }

        my_vector<T> operator-(const my_vector<T> other) const {
            my_vector<T> difference(*this);

            for (size_t i = 0; i < this->size(); ++i) {
                difference[i] -= other[i];
            }

            return difference;
        }  

        my_vector<T> operator+(const my_vector<T> other) const {
            my_vector<T> sum(*this);

            for (size_t i = 0; i < this->size(); ++i) {
                sum[i] += other[i];
            }

            return sum;
        } 

        my_vector<T> operator/(const my_vector<T> other) const {
            my_vector<T> sum(*this);

            for (size_t i = 0; i < this->size(); ++i) {
                if(other[i] != 0){
                    sum[i] /= other[i];
                }
            }

            return sum;
        } 

        my_vector<T> inverse() const {
            my_vector<T> result(*this->size());

            for (size_t i = 0; i < this->size(); ++i) {
                result[i] = 1/(*this)[i];
            }

            return result;
        }     

        void set_max_min(){
            min_random = -100;
            max_random = 100;
        }
};


template <typename T>
class my_matrix{
    private:
        int num_rows;
        int num_cols;

        std::valarray<T> matrix;

        T max_random;
        T min_random;

    public:
        my_matrix() {
            set_max_min();
        }
        my_matrix(int rows, int cols) : num_rows(rows), num_cols(cols){
            matrix = std::valarray<T>(rows*cols);
            zeros();
            set_max_min();
        }

        my_matrix(my_matrix<T>& other): num_rows(other.num_rows), num_cols(other.num_cols){
            matrix = std::valarray<T>(other.num_rows * other.num_cols);

            for(int i = 0; i < other.num_rows * other.num_cols; i++){
                matrix[i] = other.matrix[i];
            }

            num_cols = other.num_cols;
            num_rows = other.num_rows; 
            set_max_min();
        }

        void resize(int new_rows, int new_cols){
            std::valarray<T> temp(new_rows * new_cols);

            size_t copy_cols = std::min(new_cols, num_cols);
            size_t copy_rows = std::min(new_rows, num_rows);

            for(int i = 0; i < copy_cols; i++){
                std::copy(std::begin(matrix) + i * num_rows, std::begin(matrix) + i * num_rows + copy_rows, std::begin(temp) + i * new_rows);
            }

            matrix.swap(temp);

            num_cols = new_cols;
            num_rows = new_rows;
        }

        void fill(T value){
            for(int i = 0; i < num_rows * num_cols; i++){
                matrix[i] = value;
            }
        }

        void fill_rand_sym(){
            for(int row = 0; row < num_rows; row++){
                for(int col = 0; col < num_cols; col++){
                    T random_value = getRandomValue(min_random, max_random);
                    matrix[col * num_rows + row] = random_value;
                    matrix[row * num_rows + col] = random_value;
                }
            }
        }
        my_vector<T> extract_column(int col) const {
            my_vector<T> column(num_rows);

            std::copy(&matrix[col * num_rows], &matrix[(col + 1) * num_rows], column.access_raw_pointer());

            return column;
        }

        void set_column(my_vector<T> column, int col){
            std::copy(column.access_raw_pointer(), column.access_raw_pointer() + column.size(), &matrix[col * num_rows]);
        }

        void zeros(){ fill(0); }

        const int rows() const{ return num_rows; }

        const int cols() const{ return num_cols; }

        T& operator()(int row, int col){ return matrix[col * num_rows + row]; }
        const T& operator()(int row, int col) const { return matrix[col * num_rows + row]; }

        const void print() const{
            for(int row = 0; row < num_rows; row++){
                for(int col = 0; col < num_cols; col++){
                    std::cout << (*this)(row, col) << "\t";
                }
                std::cout << std::endl;
            }
        }

        my_matrix<T> operator-(const T value) const {
            my_matrix<T> result(*this);
            for(int i = 0; i < num_rows * num_cols; i++){
                result -= value;
            }
            return result;
        }

        my_vector<T> operator*(const my_vector<T>& other){
            my_vector<T> result((*this).rows());

            for(int i = 0; i < (*this).rows(); i++){
                for(int k = 0; k < (*this).cols(); k++){
                    result[i] += (*this)(i, k) * other[k];
                }
            }
            return result;
        }

        my_matrix<T> operator*(const my_matrix<T>& other){
            my_matrix<T> result((*this).rows(), other.cols());

            for(int i = 0; i < (*this).rows(); i++){
                for(int j = 0; j < other.cols(); j++){
                    for(int k = 0; k < (*this).cols(); k++){
                        result(i, j) += (*this)(i, k) * other(k, j);
                    }
                }
            }
            return result;
        }

        T* access_raw_pointer(){
            return &matrix[0];
        }

        void set_max_min(){
            min_random = -100;
            max_random = 100;
        }
};

template<typename T>
struct row_element{
    int col;
    T val;

    row_element(T val, int column) : col(column), val(val) {}; 

    bool operator<(const row_element& other) const {
        return col < other.col;
    }
};

template<typename T>
class my_sparse_matrix{
    private:
        int num_rows;
        int num_cols;

        double sparsity;

        T max_random;
        T min_random;

        std::vector<std::vector<row_element<T>>> matrix;
        my_vector<T> diagonol;
    
    public:
        my_sparse_matrix(double sparsity = 0.1) : sparsity(sparsity) { set_max_min(); }
        my_sparse_matrix(int rows, int cols, double sparsity = 0.1) : num_rows(rows), num_cols(cols), sparsity(sparsity){
            matrix.resize(rows);
            diagonol.resize(rows);
            diagonol.zeros();
            set_max_min();
        }

        void resize(int rows, int cols){
            matrix.resize(rows);
            diagonol.resize(rows);
            
            if(cols < num_cols){
            for(int row = 0; row < rows; row++){
                auto end = std::remove_if(matrix[row].begin(), matrix[row].end(), [cols](row_element<T> this_element){
                    return this_element.col >= cols;
                });
                matrix[row].erase(end, matrix[row].end());
            }
            }
            num_rows = rows;
            num_cols = cols;
        }

        void insert(T value, int row, int col){
            bool value_changed = false;
            if(row < num_rows && col < num_cols){
                for(int i = 0; i < matrix[row].size(); i++){
                    if(matrix[row][i].col == col){
                        matrix[row][i].val = value;
                        value_changed = true;
                    }
                }

                if(!value_changed){
                    matrix[row].push_back(row_element(value, col));
                    std::sort(matrix[row].begin(), matrix[row].end());
                }

                if(row == col){
                    diagonol[row] = value;
                }
            }
        }

        void fill_rand_sym(){
            for(int i = 0; i < (sparsity / 2) * num_rows * num_cols; i++){
                int rand_row = getRandomValue(0, num_rows);
                int rand_col = getRandomValue(0, num_cols);

                T random_value = getRandomValue(min_random, max_random);
                insert(random_value, rand_col, rand_row);
                if(rand_col != rand_row){
                    insert(random_value, rand_row, rand_col);
                }
            }
        }

        const int rows() const{ return num_rows; }

        const int columns() const{ return num_cols; }

        T retrieve(int row, int col) const {
            for(auto element : matrix[row]){

                if(element.col == col){
                    return element.val;
                }
            }
            return 0;
        }

        T operator()(int row, int col) const {
            return retrieve(row, col);
        }

        const my_vector<T>& diagonal() const{
            return diagonol;
        }

        operator my_matrix<T>() const {
            my_matrix<T> new_matrix(num_rows, num_cols);

            for(int row = 0; row < num_rows; row++){
                for(auto element : matrix[row]){
                    new_matrix(row, element.col) = element.val;
                }
            }
            return new_matrix;
        }

        my_vector<T> operator*(const my_vector<T>& vector) const{
            my_vector<T> product(num_rows);
            product.zeros();

            for(int row = 0; row < num_rows; row++){
                for(auto value : matrix[row]){
                    product[row] += value.val * vector[value.col];
                }
            }
            return product;
        }

        my_vector<T> omp_mvm(const my_vector<T>& vector) const{
            my_vector<T> product(num_rows);
            product.zeros();

            #pragma omp parallel for
            for(int row = 0; row < num_rows; row++){
                for(auto value : matrix[row]){
                    product[row] += value.val * vector[value.col];
                }
            }
            return product;
        }

        void print(){
            my_matrix<T> this_matrix = (*this);
            this_matrix.print();
        }

        void print_raw(){
            for(int row = 0; row < num_rows; row++){
                std::cout << "Row: " << row << "\t";
                for(auto element : matrix[row]){
                    std::cout << "(val " << element.val << " col " << element.col << ")\t";
                }
                std::cout << std::endl;
            }
        }

        void set_max_min(){
            min_random = -100;
            max_random = 100;
        }
};
}