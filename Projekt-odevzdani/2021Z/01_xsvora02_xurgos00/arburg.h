// Implementation from https://c.mql5.com/3/133/Tutorial_on_Burg_smethod_algorithm_recursion.pdf

#ifndef __ARBURG_H__
#define __ARBURG_H__

#include <math.h>
#include <vector>

void BurgAlgorithm( std::vector<double> &coeffs, const std::vector<double> &x );

#endif  // __ARBURG_H__
