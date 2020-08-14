#pragma once
#include <cmath>

#define M_PI1         3.141592653589793238462643383279502884L /* pi */

namespace Parameters {

	//Potential parameters
	static double constexpr omega = 1.0;
	static double constexpr distance = 7.0;
	static int constexpr nmax = 1;
	static double Vb = pow(omega, 2) * pow(distance, 2) / 32.0;

};