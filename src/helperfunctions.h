#pragma once
#include <cmath>
#include <array>

#include "PathBase.h"


// this function helps compare two floats considering numerical noise
bool arefloatsequal(const float &A, const float &B, const float &maxRelDiff = FLT_EPSILON);

// calculate acos with clamped values at +-1
float clamped_acos(const float &value);

////*************** generate progress information **************
void progressbar(float &processed, float &oldprocessed, const size_t &progresscounter, const size_t &sizeporevec);

const double PI = 3.141592653589793238463;