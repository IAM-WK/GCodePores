#pragma once
#include <vector>
#include <string>
#include <array>

#include "PathBase.h"
#include "helperfunctions.h"

// this class can be used to interpolate a GCode path to generate points at a given distance

class GCodeInterpolation: public PathBase
{

public:

	// this constructor generates a vector "GCodeCoordinateVector" containing the path of the printer out of a GCode file including the direction
	GCodeInterpolation();

	void setDistance(const float &distanceVal);
	void setInput(const PathVector &GCodeVector);
	// this function normalises the vector to contain only coordinates at certain points at a specified distance to the last coordinate
	void Execute();

	// returns CoordinateVector: NormalisedGCode --> this will return empty vec if badly called
	// maybe use enum datatype instead of rep. string
	PathBase::PathVector GetOutput() const;
	

private:

	
	// stores normalised coordinates (either GCodeKOS / ImageKOS)
	PathBase::PathVector NormalisedGCode;
	PathBase::PathVector GCodeVec;

	// vars to store settings
	float distance = -1.f;

	// this function generates values interpolated with const distance in normalised vector
	void interpolatepath(ChunkVector::const_iterator &CCoordIt, const PathVector::const_iterator &CChunkIt);

};
