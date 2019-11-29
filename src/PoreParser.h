#pragma once
#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <string>


class PoreParser
{

public:

	// number of colums in morpholibj output
	// #3 is sphericity
	static const unsigned int dimtuple = 26;
	typedef std::array<std::string, PoreParser::dimtuple> PoreInfoHeaderType;
	typedef std::array<float, PoreParser::dimtuple> PoreTupleType;
	typedef std::vector<PoreTupleType> PoreVecType;

	// this constructor generates a parser
	PoreParser();

	void setFilename(const std::string &Filename);
	void setSpherthres(const float &sphthres);
	void setVolthres(const unsigned int &volthres);
	void setARthres(const float &ARthres);
	void Execute();

	std::vector< std::array< float, 26>> GetOutput() const;
	unsigned int sphericitycol = 0, volumecol = 0, Xcol = 0, Ycol = 0, Zcol = 0, azimuthcol = 0, elevationcol = 0, rad12col = 0, rad13col = 0;
	unsigned int porevolume = 0, filteredporesvolume = 0;

	// histogram of pore orientation (uncorrelated with GCode)

	std::vector <unsigned int> poreorientationhisto_azim_overall;
	std::vector <unsigned int> poreorientationhisto_azim_filtered;

	std::vector <unsigned int> poreorientationhisto_elev_overall;
	std::vector <unsigned int> poreorientationhisto_elev_filtered;

	std::vector <unsigned int> poreazimuth_histo_overall_volume;
	std::vector <unsigned int> poreazimuth_histo_filtered_volume;

	std::vector <unsigned int> poreelevation_histo_overall_volume;
	std::vector <unsigned int> poreelevation_histo_filtered_volume;

	// count of keyhole pores vs gaspores#
	unsigned int gaspores = 0;
	unsigned int khpores = 0;


private:

	// list read from file
	PoreVecType porelist;
	// header from porefile
	PoreInfoHeaderType tupleheader;

	std::string porefilename;
	// sphericity threshold to sort out pores above
	float sphericitythres = -1.f, aspectratiothres = -1.f;
	unsigned int volumethres = 0;

};
