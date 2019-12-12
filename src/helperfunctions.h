#pragma once
#include <cmath>
#include <array>

#include "PathBase.h"


// this function helps compare two floats considering numerical noise
bool arefloatsequal(const float &A, const float &B, const float &maxRelDiff = FLT_EPSILON) noexcept;

// calculate acos with clamped values at +-1
float clamped_acos(const float &value) noexcept;

// calculate cross product of two arrays
std::array<double,3> crossproduct(const std::array<double,3> &v1, const std::array<double, 3> &v2) noexcept;

////*************** generate progress information **************
void progressbar(const size_t &progresscounter, const size_t &sizeporevec);

constexpr double PI = 3.141592653589793238463;

// for progressbar
static float processed = -1, oldprocessed = -1; // for progress bar