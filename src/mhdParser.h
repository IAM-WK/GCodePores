#pragma once
#include <fstream>
#include <array>
#include <string>

class mhdParser
{
	// small parser to read dimensions of µCT image from .mhd file
	// itk library can read dimensions easily, but this doesn't introduce dependencies
	// the dimensions are needed for coordinate system transformation


public:
	void setFilename(const std::string& Filename);
	void Execute();
	std::array<float, 3> imagedimensionmm = { -1.f,-1.f,-1.f };
	std::array<float, 3> pixelwidth = { -1.f,-1.f,-1.f};

private:
	std::string mhdfilename;
	std::array<float, 3> imagedimensionvoxels = { -1.f,-1.f,-1.f };
};