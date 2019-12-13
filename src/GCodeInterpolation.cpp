#include "GCodeInterpolation.h"


GCodeInterpolation::GCodeInterpolation()
{
}

void GCodeInterpolation::setDistance(const float & distanceVal)
{
	this->distance = distanceVal;
}

void GCodeInterpolation::setInput(const PathVector &GCodeVector)
{
	// estimate vector length to reserve space ---
	// copy values to new vec
	this->GCodeVec = GCodeVector;
}



void GCodeInterpolation::Execute()
{
	std::vector< std::vector< std::array<float, dim_coords> > >::const_iterator CChunkIt; // to iterate through coordinate chunks between travels
	std::vector< std::array<float, dim_coords> >::const_iterator CCoordIt;

	//float searchdistance = 0; //distance to next generated point from next point in GCodeCoordinateVector
	for (CChunkIt = this->GCodeVec.begin(); CChunkIt != this->GCodeVec.end(); ++CChunkIt) {
		if (CChunkIt->size() == 1) { // means this is only an intermediate travel move -- non relevant for printing path
			continue;
		}
		for (CCoordIt = CChunkIt->begin(); CCoordIt != CChunkIt->end(); ++CCoordIt) {
			//calculate new desired distance
			//searchdistance = this->distance - this->carry_distance;
			if (CCoordIt->at(6) >= this->distance) {
				//do interpolation
				interpolatepath(CCoordIt, CChunkIt);
				//remember carry distance for next point
			}
			else {
				if (CCoordIt == --CChunkIt->end()) { // last move of a printing chunk
					this->NormalisedGCode.back().push_back({ *CCoordIt }); //store point in a new chunk 
				}
				else if (CCoordIt == CChunkIt->begin()) { // last move was a travel
					this->NormalisedGCode.push_back({ *CCoordIt }); //store point in a new chunk 
				}
			}
		}
	}
		return;
} 


void GCodeInterpolation::interpolatepath(ChunkVector::const_iterator &CCoordIt, const PathVector::const_iterator &CChunkIt)
{
	// CCoordIt is mutable to peek at next entry for feedrate lookup
	std::array<float, dim_coords> newpathpoint; //for new values
	
	if (CCoordIt == --CChunkIt->end()) { // last move of a printing chunk
		this->NormalisedGCode.back().push_back({ *CCoordIt }); //store point in a new chunk 
		return;
	}

	if (CCoordIt == CChunkIt->begin()) { // printer has just travelled
		for (size_t j = 0; j < 3; ++j) { // copy entry
			newpathpoint[j] = CCoordIt->at(j);
			newpathpoint[j + 3] = CCoordIt->at(j + 3);
		}
		newpathpoint[6] = this->distance;
		newpathpoint[7] = (++CCoordIt)->at(7); // increment allowed because checked above if we are at the end
		newpathpoint[8] = CCoordIt->at(8); // update BED and power values from coordinate at next point
		newpathpoint[9] = CCoordIt->at(9);
		(--CCoordIt); // decrement to be correct again, use next feedrate, because feedrate is from travel move after travel otherwise
		this->NormalisedGCode.push_back({ newpathpoint }); //store point in a new chunk
	}
	else {
		for (size_t j = 0; j < 3; ++j) { // copy entry
			newpathpoint[j] = CCoordIt->at(j) + this->distance * CCoordIt->at(j + 3);
			newpathpoint[j + 3] = CCoordIt->at(j + 3);
		}
		this->NormalisedGCode.back().push_back(newpathpoint); //store point in existing chunk
	}

	if ((CCoordIt->at(6) - this->distance) < this->distance) { //distance limit reached already
		return;
	}
	for (size_t k = 1; k <= ((CCoordIt->at(6) - this->distance) / this->distance); ++k) { //interpolates between coords
		for (size_t j = 0; j < 3; ++j) { //entry of coordinates 0->X, 1->Y, 2->Z, 3->dX...
			newpathpoint[j] = this->NormalisedGCode.back().back()[j] + this->distance * CCoordIt->at(j + 3);
			newpathpoint[j + 3] = CCoordIt->at(j + 3);
		}
		newpathpoint[6] = this->distance;
		newpathpoint[7] = CCoordIt->at(7);
		newpathpoint[8] = CCoordIt->at(8);
		newpathpoint[9] = CCoordIt->at(9);
		if (CCoordIt == CChunkIt->begin()) { //printer has just travelled
			newpathpoint[7] = (++CCoordIt)->at(7); // increment allowed because checked above if we are at the end
			newpathpoint[8] = CCoordIt->at(8); // update BED and power values from coordinate at next point
			newpathpoint[9] = CCoordIt->at(9);
			(--CCoordIt); // decrement to be correct again, use next feedrate, because feedrate is from trave move after travel otherwise
		}
		this->NormalisedGCode.back().push_back(newpathpoint); //store points
	}
	return;
}

PathBase::PathVector GCodeInterpolation::GetOutput() const
{
	if (this->NormalisedGCode.size() == 0) {
		throw std::runtime_error("Error after path interpolation: Interpolated coordinate vector is empty!");
	}
	return this->NormalisedGCode;
}


