#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <array>

#include "PathBase.h"
#include "helperfunctions.h"

// this class can be used to parse a GCode file to a vector and read printing parameters from the file

class GCodeParser: public PathBase
{

public:

	// this constructor generates a parser
	GCodeParser();

	void setFilename(const std::string &Filename);
	void setSkirtoffset(const size_t &skirtoffsetval);
	// setBoundarystrings will only use a nonempty startstring, if you want to use startbyskirt give empty string
	void setBoundarystrings(const std::string &startstringval, const std::string &endstringval);
	void setCompatibility(const bool &freeformcompatval, const bool &slmcompatval);
	void Execute();
	
	// returns CoordinateVector
	PathBase::PathVector GetOutput() const;
	
	void getGCodeparams(float &extr_width, float &lay_height) const;

private:

	// stores coordinates read from GCodeFile; 
	// inner vectors contains contiguous print paths between travels, outer vector as aggregate structure
	PathBase::PathVector PrintVector;

	// vars to store settings
	std::string GCodeFilename, startstring, endstring;
	size_t skirtoffset;
	bool freeformcompat;
	bool slmcompat;
	bool startbystring = false;


	// stores metadate from GCodeFile
	size_t skirtnum = 0;	// number of skirts surrounding the object; used to estimate object start
	size_t numoflines = 0;	// length of GCodefile; used to estimate vector length
	float extrusion_width = 0.0f;
	float layer_height = 0.0f;
	float printpathlength = 0.0f; // length of printed path

	// this reads various information from GCodefile before it is processed
	void readinfos(std::ifstream &GCodeFile, std::string &line_string, const std::string &coordchars);

};
