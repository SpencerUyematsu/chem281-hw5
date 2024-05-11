#include <random>

template<typename T>
T getRandomValue(T minVal, T maxVal) {
    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Define the distribution based on the type T
    if constexpr (std::is_integral_v<T>) {
        std::uniform_int_distribution<T> dist(minVal, maxVal);
        return dist(gen); // Generate and return a random integer
    } else if constexpr (std::is_floating_point_v<T>) {
        std::uniform_real_distribution<T> dist(minVal, maxVal);
        return dist(gen); // Generate and return a random double
    } else {
        throw std::invalid_argument("Unsupported type for getRandomValue");
    }
}
