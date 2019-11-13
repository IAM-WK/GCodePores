#include "helperfunctions.h"



bool arefloatsequal(const float &A, const float &B, const float &maxRelDiff){
	// Calculate the difference.
	const float diff = fabs(A - B);
	// Find the largest of A/B abs. and assign to largest
	const float largest = (fabs(B) > fabs(A)) ? fabs(B) : fabs(A);
	if (diff <= largest * maxRelDiff) {
		return true;
	}
	return false;
}


float clamped_acos(const float & value)
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
