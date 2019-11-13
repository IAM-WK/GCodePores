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
	float distance;

	float carry_distance; // will store carry from last generated coord to next coord in PrintVector

	// this function generates values interpolated with const distance in normalised vector
	void interpolatepath(ChunkVector::const_iterator &CCoordIt, const PathVector::const_iterator &CChunkIt, const float &searchdistance);
	
	// this function accumulates carry until distance is big enough, if points are too close
	void skipcoords(const ChunkVector::const_iterator &CCoordIt, const PathVector::const_iterator &CChunkIt);

};
