#include "helperfunctions.h"



bool arefloatsequal(const float &A, const float &B, const float &maxRelDiff) noexcept
{
	// Calculate the difference.
	const float diff = fabs(A - B);
	// Find the largest of A/B abs. and assign to largest
	const float largest = (fabs(B) > fabs(A)) ? fabs(B) : fabs(A);
	if (diff <= largest * maxRelDiff) {
		return true;
	}
	return false;
}


float clamped_acos(const float & value) noexcept
{
	if (value <= -1.0) {
		return static_cast<float>(PI);
	}
	else if (value >= 1.0) {
		return static_cast<float>(0);
	}
	else {
		return acos(value);
	}
}

std::array<double, 3> crossproduct(const std::array<double, 3>& v1, const std::array<double, 3>& v2) noexcept
{
	std::array<double, 3> v3;
	v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
	v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
	v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
	return v3;
}


void progressbar(const size_t &progresscounter, const size_t &sizeporevec)
{
	oldprocessed = processed;
	processed = std::roundf(static_cast<float>(progresscounter) / static_cast<float>(sizeporevec) * 100); // fertig in prozent (0.xy)
	if (oldprocessed != processed) { // new percent values only
		std::cout << "\rprocessed " << processed << "% of pores!";
	}
}
