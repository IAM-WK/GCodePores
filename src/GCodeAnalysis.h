#pragma once
#include <vector>
#include <cmath>
#include <numeric>
#include "PathBase.h"
#include "helperfunctions.h"

class GCodeAnalysis: public PathBase
{
public:
	 
	typedef float flttype;

	void setvtkfilename(const std::string &vtkfilenamearg);

	// this calculates the angle between the line from the last point to the current point and the line from the current point to the next point
	void calcdeltaangle(const PathBase::PathVector &PathVec, const bool &usedeg, const bool &filterangles, const float &filterdist);

	//void calccurvatures(const PathBase::PathVector &PathVec);
	
	void classifypoints(const PathBase::PathVector &PathVec, const float &tol_loopclosure);
	
	// feedrate to field vector 
	void findfeedrate(const PathBase::PathVector &PathVec);


	// findreversepoint(float range) -> markiert punkte an denen auf 0 gebremst wird innerhalb eines radius (je nach dist untersch punkte in nähe)

	// PathBase::PathVector GetOutput() const;
	void calcpathlength(const PathBase::PathVector &PathVec);

	void setlayerheighttol(const float &layer_height_tol_val);
	void calclayerheights(const PathBase::PathVector& PathVec);


	// calculated vectors - field over vtk lines
	std::vector<flttype> deltaangles;
	std::vector<flttype> deltaangles_filtered;
	std::vector<flttype> curvatures; 
	std::vector<int> pathclassification;
	std::vector<int> pathclassificationwithchunks;
	std::vector<flttype> lengthofchunks; // stores length for every chunk (size() is no of chunks)
	std::vector<flttype> lengthofchunks_points; // stores length of every printpath to each point
	std::vector<flttype> feedrate;

	// chunk i has layer_height[i] : layer height is average of point height -> deal with slightly rotated layers
	std::vector<flttype> layer_height;
	// vec with all unique heights
	std::vector<flttype> heights_avail;


private:
	
	std::string vtkfilename;

	// helper to sort index / coord pairs (used for perim no classification)
	static bool paircomparator(const std::pair<float, size_t> &l, const std::pair<float, size_t> &r);

	// to decide if layer height is unique
	float layer_height_tol;
};