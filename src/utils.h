#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <string>

#include "structures.h"

template <typename T>
std::string number_to_string(T number)
{
    std::ostringstream ss;
    ss << number;
    return ss.str();
}

// trim from start
inline std::string& ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}
// trim from end
inline std::string& rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}
// trim from both ends
inline std::string& trim(std::string &s) {
    return ltrim(rtrim(s));
}

static char *copy_str(char const *src) {
  const size_t len = strlen(src) + 1;
  char *dest = (char*)malloc(len * sizeof(char));
  return strncpy(dest, src, len);
}

template <typename T>
T* alloc(int N)
{
    return new T[N];
}

template <typename T>
void dealloc(T* array)
{
    delete[] array;
}

inline bool detectLittleEndian() {
    int n = 1;
    // little endian if true
    return (*(char *)&n == 1);
}

inline void clean_level(int nel, float* volumes,
        int* neighbours, float* variables, float* old_variables,
        float* fluxes, float* step_factors, edge_neighbour* edges,
        edge* edge_mx, edge* edge_my, edge* edge_mz, edge* edge_p, edge* edge_pe, float* coords)
{
    dealloc<float>(volumes);
    dealloc<int>(neighbours);

    dealloc<float>(variables);
    dealloc<float>(old_variables);
    dealloc<float>(fluxes);
    dealloc<float>(step_factors);

    dealloc<edge_neighbour>(edges);
    dealloc<edge>(edge_mx);
    dealloc<edge>(edge_my);
    dealloc<edge>(edge_mz);
    dealloc<edge>(edge_p);
    dealloc<edge>(edge_pe);

    dealloc<float>(coords);
}

// template<typename T>
// bool isnan(const T &v) {
//     return v!=v;
// }

// template<typename T>
// bool isinf( const T &value )
// {
//     T max_value = std::numeric_limits<T>::max();
//     T min_value = - max_value;
 
//     return ! ( min_value <= value && value <= max_value );
// }

#endif
