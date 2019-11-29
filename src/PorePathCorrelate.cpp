#include "PorePathCorrelate.h"


PorePathCorrelate::PorePathCorrelate(PathBase::PathVector* const pathVecPtr, std::vector<int>* const pathclassificationwithchunks, std::vector<float>* const layerheight, std::vector<float>* const heightsavail, std::vector<float>* const lengthofchunks, PoreParser::PoreVecType* const poreVecPtr)
{
	this->PathVecPointer = pathVecPtr;
	this->PathClassificationwithChunks = pathclassificationwithchunks;
	this->LayerHeight = layerheight;
	this->HeightsAvail = heightsavail;
	this->LengthofChunks = lengthofchunks;
	this->PoreVecPointer = poreVecPtr;
}

void PorePathCorrelate::setInputSettings(const float& layerheighttol, const unsigned int& layersconsidered, const unsigned int& nextneighbours, const float& classwidth, const float& pathlengthclasswidth, const double& pixelvolume)
{
	this->LayerHeightTol = layerheighttol;
	this->LayersConsidered = layersconsidered;
	this->NextNeighbours = nextneighbours;
	this->ClassWidth = classwidth;
	this->PathLengthClassWidth = pathlengthclasswidth;
	this->PixelVolume = pixelvolume;
}

void PorePathCorrelate::setPoreCols(const unsigned int& xcol, const unsigned int& ycol, const unsigned int& zcol, const unsigned int& azimuthcol, const unsigned int& volumecol)
{
	this->XCol = xcol;
	this->YCol = ycol;
	this->ZCol = zcol;
	this->AzimuthCol = azimuthcol;
	this->VolumeCol = volumecol;
}

void PorePathCorrelate::Execute()
{

	//************************************************************************
	// initialize data
	// float: deltaangle, int: class, int: volume
	this->PoretoPathAngles.resize(this->PoreVecPointer->size());
	std::fill(this->PoretoPathAngles.begin(), this->PoretoPathAngles.end(), std::make_tuple(150.f, 0, 0));
	// vec to store how many pores a chunk contains; size is number of chunks, start with 0
	this->PoreinChunkFrequencies.resize(this->PathVecPointer->size());
	std::fill(this->PoreinChunkFrequencies.begin(), this->PoreinChunkFrequencies.end(), 0);
	// vec to store how many defect voxels a chunk contains; size is number of chunks, start with 0
	this->PoreinChunkFrequencies_VW.resize(this->PathVecPointer->size());
	std::fill(this->PoreinChunkFrequencies_VW.begin(), this->PoreinChunkFrequencies_VW.end(), 0);
	// what is the highest class found?
	this->MaxClass = *std::max_element(this->PathClassificationwithChunks->begin(), this->PathClassificationwithChunks->end());
	// vec to store number of pores in each class; eg class 5 -> highest entry = 5
	this->PoreinClassFrequencies.resize(this->MaxClass + 1);
	std::fill(this->PoreinClassFrequencies.begin(), this->PoreinClassFrequencies.end(), 0);
	// vec to store number of pores in each class; eg class 5 -> highest entry = 5
	this->PoreinClassFrequencies_VW.resize(this->MaxClass + 1);
	std::fill(this->PoreinClassFrequencies_VW.begin(), this->PoreinClassFrequencies_VW.end(), 0);
	// vec to store how many pores a layer contains
	this->PoreinLayerFrequencies.resize(this->HeightsAvail->size());
	std::fill(this->PoreinLayerFrequencies.begin(), this->PoreinLayerFrequencies.end(), 0);
	// same, but summed volume instead of number
	this->PoreinLayerFrequencies_VW.resize(this->HeightsAvail->size());
	std::fill(this->PoreinLayerFrequencies_VW.begin(), this->PoreinLayerFrequencies_VW.end(), 0);

	// temporary helper vars
	// height of current pore, minimal difference of height to layer height, azimuth angle of neighbouring path, azimuth angle of current pore, delta of this angles
	float current_pore_height, minheightdiff, GCodeAzimuth, PoreAzimuth, angledelta;
	// index of chunk next to pore
	unsigned int chunknum_next;
	// helper var: delta index of layers around nearest layer, temp index of layers around nearest layer
	int deltaindex, height_index;
	
	// end initialisation
	//************************************************************************
	// for progressbar
	size_t progresscounter = 0;
	const size_t sizeporevec = this->PoreVecPointer->size();

	//************************************************************************
	//********************* loop through pores *******************************
#pragma omp parallel for private(current_pore_height, minheightdiff, GCodeAzimuth, PoreAzimuth, angledelta, chunknum_next, deltaindex, height_index)
	for (size_t porenum = 0; porenum < sizeporevec; ++porenum) {

		// unfortunately this is not so easy in parallel ... :(
#pragma omp atomic
		++progresscounter;

#pragma omp critical (outpinfo)
		{	// show and update progressbar
			progressbar(progresscounter, sizeporevec);
		}
		//********************************************************************************************************************************************
		// reduce points to be considered by sorting out relevant layers
		// first find next layer coord for porecenter.z
		current_pore_height = (*this->PoreVecPointer)[porenum][this->ZCol];
		chunknum_next = std::numeric_limits<unsigned int>::max();
		minheightdiff = std::numeric_limits<float>::max();

		// iterate through height of chunks (layer_height)
		// if height pore - height chunk -> minimal ---> this is the next layer resp. next chunk, but other chunks in same layer are also possible
		for (unsigned int j = 0; j < this->LayerHeight->size(); ++j) {
			if (abs(current_pore_height - this->LayerHeight->at(j)) < minheightdiff) {
				chunknum_next = j;
				minheightdiff = abs(current_pore_height - this->LayerHeight->at(j));
			}
		}
		// this pore is in layer with height of LayerHeight[chunknum_next]
		// search entry in heights_avail with this layer height -> index in heights_avail is index for PoreinLayerFrequencies#
		// increment counter PoreinLayerFrequencies - RESULT
		for (unsigned int layernum = 0; layernum < this->HeightsAvail->size(); ++layernum) {
			if (arefloatsequal(this->LayerHeight->at(chunknum_next), this->HeightsAvail->at(layernum), this->LayerHeightTol)) {
#pragma omp atomic
				PoreinLayerFrequencies[layernum] += 1;
#pragma omp critical (calcvolume) // atomic only supports normal incrementation
				{
					PoreinLayerFrequencies_VW[layernum] += static_cast<size_t>((*this->PoreVecPointer)[porenum][this->VolumeCol]);
				}
			}
		}

		// next step is to reduce heights_avail to only the heights we need
		// formula to alternatingly select height +1/-1/+2/-2 ...:
		// deltaindex = static_cast<int>(round((static_cast<float>(deltaindex) * (-2.f) + 1.f) / 2.f));
		deltaindex = 0, height_index = 0;
		std::vector<float> heights_considered;

		for (unsigned int j = 0; j < this->HeightsAvail->size(); ++j) {
			if (arefloatsequal(this->HeightsAvail->at(j), this->LayerHeight->at(chunknum_next), this->LayerHeightTol)) {
				heights_considered.push_back(this->HeightsAvail->at(j)); // nearest layer itself
				// loop with deltaindex to consider neighbouring layers
				for (unsigned int numdeltalayer = 0; numdeltalayer < this->LayersConsidered; ++numdeltalayer) {
					deltaindex = static_cast<int>(round((static_cast<float>(deltaindex) * (-2.f) + 1.f) / 2.f));
					height_index = j + deltaindex;
					if ((height_index >= 0) && (height_index < this->HeightsAvail->size())) {
						heights_considered.push_back(this->HeightsAvail->at(height_index));
					}
				}
				// if height was found, then there should not be anything else to be found...
				break;
			}
		}
		
		// next we need to translate considered heights to chunknums
		// so we generate a vector with all chunknumbers to be considered
		std::vector<unsigned int> chunks_considered;

		for (unsigned int j = 0; j < this->LayerHeight->size(); ++j) {
			for (unsigned int m = 0; m < heights_considered.size(); ++m) {
				if (arefloatsequal(this->LayerHeight->at(j), heights_considered[m], this->LayerHeightTol)) {
					chunks_considered.push_back(j);
				}
			}
		}
		// end reduce neighbouring points
		// chunks considered contains every chunk to be considered
		//********************************************************************************************************************************************

//********************************************************************************************************************************************

//********************************************************************************************************************************************
// find nearest neighbours of pore
//**************************************************************
// stores X,Y,Z,dX,dY,dZ,distance to pore, pathclass, chunknumber -- of neighbourhood path points 
		std::vector<std::array<float, 9> > neighbourhood;

		// loop through remaining path points, calc distance and store to neighbourhood
		for (unsigned int chunknum_index = 0; chunknum_index < chunks_considered.size(); ++chunknum_index) {
			std::array<float, 9> current_neighbour;
			for (unsigned int point_index = 0; point_index < (*this->PathVecPointer)[chunks_considered[chunknum_index]].size(); ++point_index) {
				// transfer location and direction of point
				for (unsigned int a = 0; a < 6; ++a) {
					current_neighbour[a] = (*this->PathVecPointer)[chunks_considered[chunknum_index]][point_index][a];
				}
				// calculate distance pore to point
				current_neighbour[6] = sqrt((((*this->PoreVecPointer)[porenum][this->XCol] - (*this->PathVecPointer)[chunks_considered[chunknum_index]][point_index][0]) * ((*this->PoreVecPointer)[porenum][this->XCol] - (*this->PathVecPointer)[chunks_considered[chunknum_index]][point_index][0])) + (((*this->PoreVecPointer)[porenum][this->YCol] - (*this->PathVecPointer)[chunks_considered[chunknum_index]][point_index][1]) * ((*this->PoreVecPointer)[porenum][this->YCol] - (*this->PathVecPointer)[chunks_considered[chunknum_index]][point_index][1])) + (((*this->PoreVecPointer)[porenum][this->ZCol] - (*this->PathVecPointer)[chunks_considered[chunknum_index]][point_index][2]) * ((*this->PoreVecPointer)[porenum][this->ZCol] - (*this->PathVecPointer)[chunks_considered[chunknum_index]][point_index][2])));
				// find out pathclassification for this neighbouring point
				current_neighbour[7] = static_cast<float>((*this->PathClassificationwithChunks)[chunks_considered[chunknum_index]]);
				// store chunknumber of this very point --> needed for porespermm
				current_neighbour[8] = static_cast<float>(chunks_considered[chunknum_index]);
				neighbourhood.push_back(current_neighbour);

			}
		}
		// sort neighbourhood to pore distance (lambda function)
		std::sort(neighbourhood.begin(), neighbourhood.end(), [](const std::array<float, 9> & a, const std::array<float, 9> & b) {
			return (a[6] < b[6]);
			});

		// nearest neighbours chunknumber is needed to find out pathclassification
		// this classification is stored to freqinclasses -- RESULT
		// neighbourhood[0][7] is nearest neighbours' pathclassification
#pragma omp atomic
		PoreinClassFrequencies[static_cast<size_t>(neighbourhood[0][7])] += 1;
#pragma omp critical (pathclassification_volweight) // atomic only supports normal incrementation
		{
			PoreinClassFrequencies_VW[static_cast<int>(neighbourhood[0][7])] += static_cast<size_t>((*this->PoreVecPointer)[porenum][this->VolumeCol]);
		}
		// we know the nearest chunk to the pore
		// increment number of pores the nearest chunk contains -- RESULT
		// PoreinChunkFrequencies vector contains the number of pores for each chunk
#pragma omp atomic
		//poreinchunkfrequencies[chunknum_next] += 1;
		PoreinChunkFrequencies[static_cast<size_t>(neighbourhood[0][8])] += 1;
#pragma omp critical (poreinchunkfreq_volweighted) // atomic only supports normal incrementation
		{
			PoreinChunkFrequencies_VW[static_cast<size_t>(neighbourhood[0][8])] += static_cast<size_t>((*this->PoreVecPointer)[porenum][this->VolumeCol]);
		}

		//************************************ Orientation Analysis ************************************
		// calculate average direction of GCode in pore neighbourhood
		std::array<float, 3 > neighbour_direction; neighbour_direction.fill(0);
		for (unsigned int neighbour_index = 0; neighbour_index < this->NextNeighbours; ++neighbour_index) {
			for (size_t h = 0; h < 3; ++h) {
				neighbour_direction[h] = neighbour_direction[h] + neighbourhood[neighbour_index][h + 3];
			}
		}
		for (int h = 0; h < 3; ++h) {
			neighbour_direction[h] = neighbour_direction[h] / static_cast<float>(this->NextNeighbours);
		}

		// convert GCode direction to azimuth/elevation system
		// GCode will usually have elevation = 0°
		// so elevation of pores cannot be correlated with GCode
		// azimuth = 180 / pi * atan2(dy, dx);
		GCodeAzimuth = (180.f / static_cast<float>(PI)) * atan2(neighbour_direction[1], neighbour_direction[0]);

		// azimuth should be between 0 - 180, only one sense of direction
		if (GCodeAzimuth < 0.f) {
			GCodeAzimuth = GCodeAzimuth + 180.f;
		}

		// PoreAzimuth 0 - 180
		PoreAzimuth = (*this->PoreVecPointer)[porenum][this->AzimuthCol];
		if (PoreAzimuth < 0.f) {
			PoreAzimuth = PoreAzimuth + 180.f;
		}



		// calculate delta angle between pore and gcode
		angledelta = abs(GCodeAzimuth - PoreAzimuth);
		// delta should be between 0 - 90 : angles between 90 - 180 should be converted to smaller angle
		if (angledelta > 90.f) {
			angledelta = 180.f - angledelta;
		}

		// store delta phi (classified or do that later on?)
		PoretoPathAngles[porenum] = std::make_tuple(angledelta, static_cast<int>(neighbourhood[0][7]), static_cast<int>((*this->PoreVecPointer)[porenum][this->VolumeCol]));
		//std::cout << "calced delta " << angledelta << std::endl;
		//********************************* End Orientation Analysis ***********************************
	}

	// break line after progresscounter is done
	std::cout << "\n";
	// now do a classification on deltaangles
	// 1. sort vec
	// 2. specify class width
	// loop
	//	3. count objects below i*class width
	//	4. if object exceeds i*class width-> ++i (next class)
	//	5. store in a array with 180/class width entries the number of objects for each i
	// 6. write this to a output file as: i*class width/2 || num angles ----> histogram
	// sort deltaangles by angle value;

	std::sort(this->PoretoPathAngles.begin(), this->PoretoPathAngles.end(), [](const std::tuple<float, int, int>& left, const std::tuple<float, int, int>& right) {
		return std::get<0>(left) < std::get<0>(right);
		});

	// helper variables
	unsigned int classnumber = 1, currentclassfrequency = 0, currentclassvolume = 0, perimcurrentclassvolume = 0, hatchingcurrentclassvolume = 0, perimclassnumber = 1, perimcurrentclassfrequency = 0, hatchclassnumber = 1, hatchcurrentclassfrequency = 0;

	for (unsigned int anglenum = 0; anglenum < this->PoretoPathAngles.size(); ++anglenum) {
		// overall 
		if (std::get<0>(this->PoretoPathAngles[anglenum]) < (static_cast<float>(classnumber) * this->ClassWidth)) {
			++currentclassfrequency;
			currentclassvolume += std::get<2>(this->PoretoPathAngles[anglenum]);
		}
		else if (std::get<0>(this->PoretoPathAngles[anglenum]) >= (static_cast<float>(classnumber) * this->ClassWidth)) {
			PoretoPathAngle_Histo.push_back(currentclassfrequency);
			PoretoPathAngle_Histo_VW.push_back(currentclassvolume);
			++classnumber;
			currentclassfrequency = 1; // 1 was found just now
			currentclassvolume = std::get<2>(this->PoretoPathAngles[anglenum]);
		}
		// perim
		if ((std::get<1>(this->PoretoPathAngles[anglenum]) != 0) && (std::get<0>(this->PoretoPathAngles[anglenum]) < (static_cast<float>(perimclassnumber) * this->ClassWidth))) {
			++perimcurrentclassfrequency;
			perimcurrentclassvolume += std::get<2>(this->PoretoPathAngles[anglenum]);
		}
		else if ((std::get<1>(this->PoretoPathAngles[anglenum]) != 0) && (std::get<0>(this->PoretoPathAngles[anglenum]) >= (static_cast<float>(perimclassnumber) * this->ClassWidth))) {
			PoretoPathAngle_Histo_Perim.push_back(perimcurrentclassfrequency);
			PoretoPathAngle_Histo_Perim_VW.push_back(perimcurrentclassvolume);
			++perimclassnumber;
			perimcurrentclassfrequency = 1; // 1 was found just now
			perimcurrentclassvolume = std::get<2>(this->PoretoPathAngles[anglenum]);
		}
		// hatching
		if ((std::get<1>(this->PoretoPathAngles[anglenum]) == 0) && (std::get<0>(this->PoretoPathAngles[anglenum]) < (static_cast<float>(hatchclassnumber) * this->ClassWidth))) {
			++hatchcurrentclassfrequency;
			hatchingcurrentclassvolume += std::get<2>(this->PoretoPathAngles[anglenum]);
		}
		else if ((std::get<1>(this->PoretoPathAngles[anglenum]) == 0) && (std::get<0>(this->PoretoPathAngles[anglenum]) >= (static_cast<float>(hatchclassnumber) * this->ClassWidth))) {
			PoretoPathAngle_Histo_Hatching.push_back(hatchcurrentclassfrequency);
			PoretoPathAngle_Histo_Hatching_VW.push_back(hatchingcurrentclassvolume);
			++hatchclassnumber;
			hatchcurrentclassfrequency = 1; // 1 was found just now
			hatchingcurrentclassvolume = std::get<2>(this->PoretoPathAngles[anglenum]);
		}
	}
	PoretoPathAngle_Histo.push_back(currentclassfrequency);
	PoretoPathAngle_Histo_Perim.push_back(perimcurrentclassfrequency);
	PoretoPathAngle_Histo_Hatching.push_back(hatchcurrentclassfrequency);
	PoretoPathAngle_Histo_VW.push_back(currentclassvolume);
	PoretoPathAngle_Histo_Hatching_VW.push_back(hatchingcurrentclassvolume);
	PoretoPathAngle_Histo_Perim_VW.push_back(perimcurrentclassvolume);
	
	// convert values to mm^3
	PoretoPathAngle_Histo_VW_cmm.resize(PoretoPathAngle_Histo_VW.size());
	PoretoPathAngle_Histo_Perim_VW_cmm.resize(PoretoPathAngle_Histo_Hatching_VW.size());
	PoretoPathAngle_Histo_Hatching_VW_cmm.resize(PoretoPathAngle_Histo_Perim_VW.size());
	// normal for loop...
	for (unsigned int index = 0; index < PoretoPathAngle_Histo_VW_cmm.size(); ++index) {
		PoretoPathAngle_Histo_VW_cmm[index] = PoretoPathAngle_Histo_VW[index] * this->PixelVolume;
		std::cout << "\n" << PoretoPathAngle_Histo_VW_cmm[index];
	}
	for (unsigned int index = 0; index < PoretoPathAngle_Histo_Perim_VW_cmm.size(); ++index) {
		PoretoPathAngle_Histo_Perim_VW_cmm[index] = PoretoPathAngle_Histo_Perim_VW[index] * this->PixelVolume;
	}
	for (unsigned int index = 0; index < PoretoPathAngle_Histo_Hatching_VW_cmm.size(); ++index) {
		PoretoPathAngle_Histo_Hatching_VW_cmm[index] = PoretoPathAngle_Histo_Hatching_VW[index] * this->PixelVolume;
	}
	
	std::cout << "class width = " << this->ClassWidth << "\n";
	std::cout << "class number = " << classnumber << "\n";
	// output to console for control....
	for (int i = 0; i < PoretoPathAngle_Histo.size(); ++i) {
		std::cout << this->ClassWidth * i << " - " << this->ClassWidth * (i + 1) << " : " << PoretoPathAngle_Histo[i] << "\n";
		this->classwidthvecangles.push_back(std::to_string(this->ClassWidth * i) + " - " + std::to_string(this->ClassWidth * (i + 1)));
	}


	//***************************************************************************************************
	// we have a vector with the frequencies of pores for each chunk
	// classify lengths of chunks, so that similar chunklengths go in the same class
	// chunk i has lengthofchunks[i]; poreinchunkfrequencies gives num of pores in them
	// first: which chunks go in which class?
	// number of classes = max lengthofchunks / classwidth +1
	const float pathlengthnumofclasses = (*(std::max_element(this->LengthofChunks->begin(), this->LengthofChunks->end())) / this->PathLengthClassWidth) + 1;
	// stores pores per mm for each pathlengthclass
	this->PoresperPathlength_Histo.resize(static_cast<size_t>(pathlengthnumofclasses));
	std::fill(this->PoresperPathlength_Histo.begin(), this->PoresperPathlength_Histo.end(), 0.f);
	// stores pores per mm for each pathlengthclass - voxel weighted
	this->PoresperPathlength_Histo_VW.resize(static_cast<size_t>(pathlengthnumofclasses));
	std::fill(this->PoresperPathlength_Histo_VW.begin(), this->PoresperPathlength_Histo_VW.end(), 0.f);
	// stores length of all chunks in the same class
	this->pathlengthsinclass.resize(static_cast<size_t>(pathlengthnumofclasses));
	std::fill(this->pathlengthsinclass.begin(), this->pathlengthsinclass.end(), 0.f);


	// classify chunk in a pathlength; iterate through classes and take every chunk that fits in class
	for (unsigned int pathlengthclass = 0; pathlengthclass < this->PoresperPathlength_Histo.size(); ++pathlengthclass) {
		for (unsigned int chunknumber = 0; chunknumber < this->LengthofChunks->size(); ++chunknumber) {
			// check if this chunk fits in current class
			if ((this->LengthofChunks->at(chunknumber) >= static_cast<float>(pathlengthclass) * this->PathLengthClassWidth) && (this->LengthofChunks->at(chunknumber) < static_cast<float>(pathlengthclass + 1) * this->PathLengthClassWidth)) {
				// add pores to pathlengthshisto[pathlenghtsclass] and add length to pathlengthsinclass for normalising at the end
				this->PoresperPathlength_Histo[pathlengthclass] += static_cast<float>(PoreinChunkFrequencies[chunknumber]);
				pathlengthsinclass[pathlengthclass] += this->LengthofChunks->at(chunknumber);
				this->PoresperPathlength_Histo_VW[pathlengthclass] += PoreinChunkFrequencies_VW[chunknumber];
			}
		}
		// finished iterating over all chunks, now normalise pores in chunk to pores per mm
		// make sure there was at least one pore in chunk class
		if (!arefloatsequal(this->PoresperPathlength_Histo[pathlengthclass], 0.f)) {
			this->PoresperPathlength_Histo[pathlengthclass] = this->PoresperPathlength_Histo[pathlengthclass] / pathlengthsinclass[pathlengthclass];
		}
		// make sure there was at least one pore in chunk class
		if (this->PoresperPathlength_Histo_VW[pathlengthclass] != 0) {
			this->PoresperPathlength_Histo_VW[pathlengthclass] = static_cast<float>(static_cast<float>(this->PoresperPathlength_Histo_VW[pathlengthclass]) / pathlengthsinclass[pathlengthclass]);
		}
	}

	// output the pores per mm as a histogram (sorted)
	std::cout << "Histogram: pores per mm over length of the paths" << std::endl;
	for (int i = 0; i < this->PoresperPathlength_Histo.size(); ++i) {
		std::cout << this->PathLengthClassWidth * i << " - " << this->PathLengthClassWidth* (i + 1) << " : " << this->PoresperPathlength_Histo[i] << std::endl;
		classwidthvecpathlength.push_back(std::to_string(this->PathLengthClassWidth * i) + " - " + std::to_string(this->PathLengthClassWidth * (i + 1)));
	}	
}

