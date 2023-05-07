//
// Created by kaixin on 5/4/23.
//

#ifndef DESCOL_SET_OPERATION_H
#define DESCOL_SET_OPERATION_H

#include <immintrin.h>
#include <cstdint>

int intersect_simd4x(int *set_a, int size_a, int *set_b, int size_b, int *set_c);

#endif //DESCOL_SET_OPERATION_H
