#pragma once
#include <vector>
#include <algorithm>
#include <limits>
#include "PathBase.h"
#include "PoreParser.h"
#include "helperfunctions.h"

class PorePathCorrelate
{

public:

	// this constructor generates a parser
	PorePathCorrelate(PathBase::PathVector* const pathVecPtr, std::vector<int>* const pathclassificationwithchunks, std::vector<float>* const layerheight, std::vector<float>* const heightsavail, std::vector<float>* const lengthofchunks, PoreParser::PoreVecType* const poreVecPtr);
	PorePathCorrelate() = delete;

	void setInputSettings(const float& layerheighttol, const unsigned int& layersconsidered, const unsigned int& nextneighbours, const float& classwidth, const float& pathlengthclasswidth, const double& pixelvolume);
	void setPoreCols(const unsigned int& xcol, const unsigned int& ycol, const unsigned int& zcol, const unsigned int& azimuthcol, const unsigned int& volumecol);

	void Execute();



	// statistic vectors
	// float: deltaangle, int: class, int: volume
	std::vector<std::tuple<float, int, int>> PoretoPathAngles;
	// vec to store how many pores a chunk contains; size is number of chunks, start with 0
	std::vector<size_t> PoreinChunkFrequencies;
	// vec to store how many defect voxels a chunk contains; size is number of chunks, start with 0
	std::vector<size_t> PoreinChunkFrequencies_VW;
	// what is the highest class found?
	size_t MaxClass = 0;
	// vec to store number of pores in each class; eg class 5 -> highest entry = 5
	std::vector<size_t> PoreinClassFrequencies;
	// vec to store number of pores in each class; eg class 5 -> highest entry = 5
	std::vector<size_t> PoreinClassFrequencies_VW;
	// vec to store how many pores a layer contains
	std::vector<size_t> PoreinLayerFrequencies;
	// same, but summed volume instead of number
	std::vector<size_t> PoreinLayerFrequencies_VW;
	// classwidth vectors for deltaangles and pathlengthclassification (for description e.g. 10deg, 20deg...)
	std::vector<std::string> classwidthvecangles;
	std::vector<std::string> classwidthvecpathlength;
	// histogram result vectors
	std::vector <unsigned int> PoretoPathAngle_Histo;
	std::vector <unsigned int> PoretoPathAngle_Histo_Perim;
	std::vector <unsigned int> PoretoPathAngle_Histo_Hatching;
	std::vector <unsigned int> PoretoPathAngle_Histo_VW;
	std::vector <unsigned int> PoretoPathAngle_Histo_Perim_VW;
	std::vector <unsigned int> PoretoPathAngle_Histo_Hatching_VW;
	std::vector<double> PoretoPathAngle_Histo_VW_cmm;
	std::vector<double> PoretoPathAngle_Histo_Perim_VW_cmm;
	std::vector<double> PoretoPathAngle_Histo_Hatching_VW_cmm;
	std::vector<float> PoresperPathlength_Histo;
	std::vector<float> PoresperPathlength_Histo_VW;


private:
	PathBase::PathVector* PathVecPointer = 0;
	PoreParser::PoreVecType* PoreVecPointer = 0;
	
	// vec with classifications numbers for chunks
	std::vector<int>* PathClassificationwithChunks = 0;
	// chunk i has layer_height[i] : layer height is average of point height -> deal with slightly rotated layers
	std::vector<float>* LayerHeight = 0;
	// vec with all unique heights
	std::vector<float>* HeightsAvail = 0;
	// vec with lengths of all chunks
	std::vector<float>* LengthofChunks = 0;

	unsigned int XCol = std::numeric_limits<unsigned int>::max();
	unsigned int YCol = std::numeric_limits<unsigned int>::max();
	unsigned int ZCol = std::numeric_limits<unsigned int>::max();
	unsigned int AzimuthCol = std::numeric_limits<unsigned int>::max();;
	unsigned int VolumeCol = std::numeric_limits<unsigned int>::max();;
	float LayerHeightTol = std::numeric_limits<float>::min();
	unsigned int LayersConsidered = std::numeric_limits<unsigned int>::max();;
	unsigned int NextNeighbours = std::numeric_limits<unsigned int>::max();;
	float ClassWidth = std::numeric_limits<float>::min();
	float PathLengthClassWidth = std::numeric_limits<float>::min();
	double PixelVolume = std::numeric_limits<double>::min();

	// stores length of all chunks in the same class
	std::vector<float> pathlengthsinclass;
};
