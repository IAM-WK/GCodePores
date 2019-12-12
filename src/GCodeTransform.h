#pragma once
#include <vector>
#include <string>
#include <array>

#include "PathBase.h"
#include "helperfunctions.h"

// this class can be used to do a geometric transformation of a path vector

class GCodeTransform : public PathBase
{

public:

	GCodeTransform(); //CTOR
	
	// transform origin of GCode to origin of printing in CT Data, eg transform vector + rotation
	void setOrigin(const std::vector<float> &ImageOriginVal) noexcept;
	void setAngles(const float &ImageAngleGammaVal, const float &ImageAngleBetaVal, const float &ImageAngleAlphaVal) noexcept;
	void setImagedimensions(const float &y_image_lengthVal, const float &z_image_lengthVal) noexcept;
	void setInput(const PathVector &GCodeVector) noexcept;
	void Execute();
	
	// returns CoordinateVector: ImageCoordinateGCode --> this will return empty vec if badly called
	// maybe use enum datatype instead of rep. string
	PathBase::PathVector GetOutput() const;

private:
	
	std::vector<float> ImageOrigin;
	float ImageAngleGamma = 0.f, ImageAngleBeta = 0.f, ImageAngleAlpha = 0.f, y_image_length = 0.f, z_image_length = 0.f;

	// stores transformed coordinates (imagecoordinatesystem)
	PathBase::PathVector ImageCoordinateGCode;

};
