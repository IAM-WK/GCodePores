#pragma once
#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <string>


class PoreParser
{

public:


	// this constructor generates a parser
	PoreParser();

	void setFilename(const std::string &Filename);
	void setSpherthres(const float &sphthres);
	void setVolthres(const unsigned int &volthres);
	void setARthres(const float &ARthres);
	void Execute();

	std::vector< std::array< float, 26>> GetOutput() const;
	unsigned int sphericitycol, volumecol, Xcol, Ycol, Zcol, azimuthcol, elevationcol,rad12col, rad13col;
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
	// number of colums in morpholibj output
	// #3 is sphericity
	static const unsigned int dimtuple = 26;
	// list read from file
	std::vector< std::array<float, PoreParser::dimtuple>> porelist;
	// header from porefile
	std::array<std::string, PoreParser::dimtuple> tupleheader;

	std::string porefilename;
	// sphericity threshold to sort out pores above
	float sphericitythres, aspectratiothres;
	unsigned int volumethres;

};
