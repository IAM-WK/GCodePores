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

	float searchdistance = 0; //distance to next generated point from next point in GCodeCoordinateVector
	this->carry_distance = 0;
	for (CChunkIt = this->GCodeVec.begin(); CChunkIt != this->GCodeVec.end(); ++CChunkIt) {
		if (CChunkIt->size() == 1) { // means this is only an intermediate travel move -- non relevant for printing path
			continue;
		}
		for (CCoordIt = CChunkIt->begin(); CCoordIt != CChunkIt->end(); ++CCoordIt) {
			//calculate new desired distance
			searchdistance = this->distance - this->carry_distance;
			if (CCoordIt->at(6) >= searchdistance) {
				//do interpolation
				interpolatepath(CCoordIt, CChunkIt, searchdistance);
				//remember carry distance for next point
			}
			else {
				skipcoords(CCoordIt, CChunkIt);
				//remember carry distance
			}
		}
	}
		return;
} 


void GCodeInterpolation::interpolatepath(ChunkVector::const_iterator &CCoordIt, const PathVector::const_iterator &CChunkIt, const float &searchdistance)
{
	// CCoordIt is mutable to peek at next entry for feedrate lookup
	std::array<float, dim_coords> newpathpoint; //for new values
	if (CCoordIt == --CChunkIt->end()) { // printer is about to travel !!!!! ->end() is past the end so decrement!
										 // dont update carry - no printing in this 
		return; // do nothing since there is nothing printed within this distance
	}

	for (size_t j = 0; j < 3; ++j) { // first: generate entry to searchdistance
		newpathpoint[j] = CCoordIt->at(j) + searchdistance * CCoordIt->at(j + 3);
		newpathpoint[j + 3] = CCoordIt->at(j + 3);
	}
	//std::cout << "CCoordIt->at(8);" << CCoordIt->at(8) << std::endl;
	newpathpoint[6] = this->distance;
	newpathpoint[7] = CCoordIt->at(7);
	newpathpoint[8] = CCoordIt->at(8);
	newpathpoint[9] = CCoordIt->at(9);
	if (CCoordIt == CChunkIt->begin()) { //printer has just travelled
		newpathpoint[7] = (++CCoordIt)->at(7); // increment allowed because checked above if we are at the end
		newpathpoint[8] = CCoordIt->at(8); // update BED and power values from coordinate at next point
		newpathpoint[9] = CCoordIt->at(9);
		(--CCoordIt); // decrement to be correct again, use next feedrate, because feedrate is from travel move after travel otherwise
		this->NormalisedGCode.push_back({ newpathpoint }); //store point in a new chunk
	}
	else {
		this->NormalisedGCode.back().push_back(newpathpoint); //store point in existing chunk
	}

	if ((CCoordIt->at(6) - searchdistance) < this->distance) { //distance limit reached already
		this->carry_distance = (CCoordIt->at(6) - searchdistance);
		return;
	}
	for (size_t k = 1; k <= ((CCoordIt->at(6) - searchdistance) / this->distance); ++k) { //interpolates between coords
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
		this->carry_distance = (CCoordIt->at(6) - searchdistance) - k * this->distance; //update carry
	}
	return;
}

void GCodeInterpolation::skipcoords(const ChunkVector::const_iterator &CCoordIt, const PathVector::const_iterator &CChunkIt)
{
	if (CCoordIt == CChunkIt->begin()) { // last move was a travel
		this->NormalisedGCode.push_back({ *CCoordIt }); //store point in a new chunk 
										   // this is needed because interpolation does not recognise move after travel
		this->carry_distance = 0; //reset carry as new point is stored

	}
	else if (CCoordIt != --CChunkIt->end()) { // decrement because end() is already "past the end"
											  //means next move is not a travel move
		this->carry_distance += CCoordIt->at(6);
	}
	// if only a travel move --> just skip it
	return;
}


PathBase::PathVector GCodeInterpolation::GetOutput() const
{
	if (this->NormalisedGCode.size() == 0) {
		throw std::runtime_error("Error after path interpolation: Interpolated coordinate vector is empty!");
	}
	return this->NormalisedGCode;
}


