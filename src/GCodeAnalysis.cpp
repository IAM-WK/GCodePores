#include "GCodeAnalysis.h"

void GCodeAnalysis::setvtkfilename(const std::string & vtkfilenamearg)
{
	this->vtkfilename = vtkfilenamearg;
}



void GCodeAnalysis::calcdeltaangle(const PathBase::PathVector & PathVec, const bool &usedeg, const bool &filterangles, const float & filterdist)
{
	// save disttonext in vec
	// calc filteredangles
	// anglei = anglei-1+anglei+anglei+1, je nach pathlength
	float deltaangle;
	std::vector<float> pathlengths;
	this->deltaangles.push_back(0.0); // this is a dummy value... first value is somehow arbitrary
	pathlengths.push_back(PathVec[0][0][6]);
	for (std::size_t m = 0; m < PathVec.size(); ++m) {
		for (std::size_t i = 0; i < PathVec[m].size(); ++i) { 
			if (m == 0 && i == 0) {
				continue;
			}
			// check for begin of pathitem
			if (i == 0) { // -> begin of pathitem
				deltaangle = clamped_acos((PathVec[m][0][3] * PathVec[m-1].back()[3]) + (PathVec[m][0][4] * PathVec[m-1].back()[4]) + (PathVec[m][0][5] * PathVec[m-1].back()[5]) );
				this->deltaangles.push_back(deltaangle);
				pathlengths.push_back(PathVec[m][0][6]);
			}
			else { // normal calculation: possible to access i-1
				// calc dx2*dx1+dy2*dy1+dz2*dz1 / norm1*norm2 // norm already 1
				deltaangle = clamped_acos((PathVec[m][i][3] * PathVec[m][i-1][3]) + (PathVec[m][i][4] * PathVec[m][i-1][4]) + (PathVec[m][i][5] * PathVec[m][i-1][5]));
				if (!usedeg) { // save values as radians
					this->deltaangles.push_back(deltaangle);
				}
				else if (usedeg) { // save as deg
					deltaangle = (deltaangle / (2.f * static_cast<float>(PI)) ) * 360.0f; // convert to deg
					this->deltaangles.push_back(deltaangle);
				}
				pathlengths.push_back(PathVec[m][i][6]);
			}

		}
	}
	if (filterangles) {
		this->deltaangles_filtered = this->deltaangles; // copy
		float carry_distance = 0, carry_angle = 0, start_distance = 0;
		size_t temp_i = 0;
		int deltaindex = 0;
		// filter angles for a value over a minimum distance
		for (size_t i = 0; i < this->deltaangles.size(); ++i) {
			// skip points until there are enough points to filter
			if (start_distance < filterdist) {
				start_distance += pathlengths[i];
				continue; // skip to next point
			}
			deltaindex = 0;
			while (carry_distance < filterdist) {
				temp_i = i + deltaindex; // take alternatingly one index ahead and after to the angle
				carry_distance += pathlengths[temp_i];
				carry_angle += this->deltaangles[temp_i];
				deltaindex = static_cast<int>(round((static_cast<float>(deltaindex) * (-2.f) + 1.f) / 2.f));
			}
			carry_distance = 0; // reset carry
			this->deltaangles_filtered[i] = carry_angle;
			carry_angle = 0; // reset angle
		}
		// write to vtkfile
		std::cout << "writing deltaangles_filtered to vtk file!\n";
		PathBase::outputcontourvtkfile(this->deltaangles_filtered, this->vtkfilename, "deltaangles_filtered");
	}
	// write to vtkfile
	std::cout << "writing deltaangles to vtk file!\n";
	PathBase::outputcontourvtkfile(this->deltaangles, this->vtkfilename, "deltaangles");
}

//
//void GCodeAnalysis::calccurvatures(const PathBase::PathVector & PathVec)
//{
//	
//}

void GCodeAnalysis::classifypoints(const PathBase::PathVector & PathVec, const float &tol_loopclosure)
{
	// analyse chunks in PathVec
	// WORKS ONLY FOR CONNECTED OBJECTS!!!!!
	// ALTERNATIVE ALGO:
	// 1 FIND CONVEX HULL -> Perim i
	// 2 REMOVE FOUND POINTS
	// 3 GOTO 1 UNTIL ONLY HATCHING LEFT
	// ALTERNATIVE
	// 1 Startpunkt suchen???
	// 2 Punkt suchen?
	// 3. Punkt: rotiere Linie von Interpolationsabstand in Richtung gegen Uhrzeigersinn (je nachdem wierum man um das Objekt läuft) bis sie auf nächsten Punkt trifft
	// 4. von dem Punkt aus: rotiere L um P3 gg Uhrzeigersinn bis man auf P5 trifft 
	// usw. (muss noch umformuliert werden)
	// if abs(startpoint-endpoint) < tol_loopclosure -> perimeter
	// maybe addionally check if direction changes significantly over chunklength 
	// freeformer infill seems to be straight lines only

	// define helper variables
	float loopclosure_x = 0, loopclosure_y = 0, loopclosure_z = 0, loopclosure_ges;
	float current_layer_height = 0.f; bool firstrun = true; // bool to start new layer at first loop
	std::vector< std::vector< std::pair<float,size_t> > > layered_chunks;// contains chunknumber at which a new layer starts
	for (std::size_t i = 0; i < PathVec.size(); ++i) {
		// start new layer for this chunk if height has changed
		if ((!arefloatsequal(current_layer_height, PathVec[i][0][2])) || firstrun) {
			layered_chunks.push_back({ std::make_pair(0.f,i) });
			firstrun = false;
		}
		else {
			layered_chunks.back().push_back(std::make_pair(0.f, i));
		}
		current_layer_height = PathVec[i][0][2];
		loopclosure_x = PathVec[i][(PathVec[i].size() - 1)][0] - PathVec[i][0][0];
		loopclosure_y = PathVec[i][(PathVec[i].size() - 1)][1] - PathVec[i][0][1];
		loopclosure_z = PathVec[i][(PathVec[i].size() - 1)][2] - PathVec[i][0][2];
		loopclosure_ges = sqrt(loopclosure_x * loopclosure_x + loopclosure_y * loopclosure_y + loopclosure_z * loopclosure_z);
		// loopclosure is first metric
		// sometimes there are hatchings with length < tol_loopclosure
		// therefore there is a second metric to sort them out:
		// does the direction ever change? if yes -> is really a perim!
		std::array<float, 3> pathdirection_compare = {PathVec[i][0][3], PathVec[i][0][4], PathVec[i][0][5]}; // initial value of direction for this pathchunk
		bool directionchanges = false;
		for (std::size_t j = 1; j < PathVec[i].size(); ++j) {
			// is x/y-direction the same as in initial? z-direction is omitted since this value shouldn't change
			if (!arefloatsequal(pathdirection_compare[0], PathVec[i][j][3]) || !arefloatsequal(pathdirection_compare[1], PathVec[i][j][4])) {
				// if this ever happens during analysing the chunk -> is really a perim
				// if never -> is a short hatching
				directionchanges = true;
			}
		}

		if (loopclosure_ges < tol_loopclosure && directionchanges) {
			// startpoint and endpoint are closer than tolerance
			// --> detected perim --> set to perimeter class
			// **** Todo: check if direction changes in this chunk. if not --> hatching
			for (std::size_t j = 0; j < PathVec[i].size(); ++j) {
				this->pathclassification.push_back(1);
			}
		}
		else {
			// ends are not at same location
			// detected hatching/infill --> set to hatching class
			for (std::size_t j = 0; j < PathVec[i].size(); ++j) {
				this->pathclassification.push_back(0);
			}
		}
	}

	// class 0 = hatching, 1 = perim
	// wanted: 1 = perim outer (1), 2 = perim 2 ...
	// classify: outer perim has highest abs coords; for every layer
	std::vector<float> maxcoords;
	std::size_t pointnum = 0;
	for (std::size_t k = 0; k < PathVec.size(); ++k) {
		// k iterates chunks
		float maxcoord = 0;
		for (std::size_t n = 0; n < PathVec[k].size(); ++n) {
			/*if (abs(PathVec[k][n][0]) > maxcoord) {
				maxcoord = abs(PathVec[k][n][0]);
			}
			if (abs(PathVec[k][n][1]) > maxcoord) {
				maxcoord = abs(PathVec[k][n][1]);
			}*/
			/*if (sqrt(abs(PathVec[k][n][0])*abs(PathVec[k][n][0])+ abs(PathVec[k][n][1])*abs(PathVec[k][n][1])) > maxcoord) {
				maxcoord = sqrt(abs(PathVec[k][n][0])*abs(PathVec[k][n][0]) + abs(PathVec[k][n][1])*abs(PathVec[k][n][1]));
			}*/
			// This methods have their problems with more than one structure ....
			// TODO: implement something like giftwrapping 
			if (abs(PathVec[k][n][0]) > maxcoord) {
				maxcoord = abs(PathVec[k][n][0]);
			}
			++pointnum;
		}
		if (this->pathclassification[pointnum-(PathVec[k].size()-1)] != 0) { // only if perim : check first point in chunk
			maxcoords.push_back(maxcoord); // biggest coord in this chunk
		}
		else {
			maxcoords.push_back(0); // zero if hatching
		}
	}
	// now we got the biggest coords in every perim
	// save biggest coords to pairs in layered_chunks
	size_t m = 0; // current index in maxcoords
	for (std::size_t k = 0; k < layered_chunks.size(); ++k) {
		for (std::size_t l = 0; l < layered_chunks[k].size(); ++l) {
			layered_chunks[k][l].first = maxcoords[m++];
		}		
	}

	// now we need to sort every inner vector in layered_chunks after maxcoords
	// index order = perim order, sparing zeros in maxcoord(=hatching)
	for (std::size_t layernum = 0; layernum < layered_chunks.size(); ++layernum) {
		// sort chunks in layer k based on maxcoord, std::sort will sort .first then .second
		std::sort(layered_chunks[layernum].begin(), layered_chunks[layernum].end(), GCodeAnalysis::paircomparator);
		// write indizes sorted to pathclassification
		int perim_number = 1; // current perim no
		for (std::size_t chunknum = 0; chunknum < layered_chunks[layernum].size(); ++chunknum) {
 			// highest coord = perim 1 ...
			// convert chunk index to start point index in field
			size_t point_index = 0;
			for (size_t chunkindex = 0; chunkindex < layered_chunks[layernum][chunknum].second; ++chunkindex) {
				point_index += PathVec[chunkindex].size();
			}
			for (size_t point_in_chunk = 0; point_in_chunk < PathVec[layered_chunks[layernum][chunknum].second].size(); ++point_in_chunk) {
				if (this->pathclassification[point_index] != 0) { // operate only in perims
					this->pathclassification[point_index] = perim_number;
					// den chunk muss man jetzt wieder in punkte rueckuebersetzen
					// und dann auch wirklichen allen punkten die perim number zuweisen
				}
				++point_index;
			}
			this->pathclassificationwithchunks.push_back(this->pathclassification[point_index-1]);
			++perim_number;
		}
	}
	
	//// write to vtkfile
	//std::cout << "writing pathclassifications to vtk file!\n";
	//PathBase::outputcontourvtkfile(this->pathclassification, this->vtkfilename, "pathclassification");
}

void GCodeAnalysis::findfeedrate(const PathBase::PathVector& PathVec)
{
	for (unsigned int i = 0; i < PathVec.size(); ++i) {
		for (unsigned int j = 0; j < PathVec[i].size(); ++j) {
			feedrate.push_back(PathVec[i][j][7]);
		}
	}
}

void GCodeAnalysis::calcpathlength(const PathBase::PathVector & PathVec)
{
	

	// calc length of every chunk
	for (unsigned int i = 0; i < PathVec.size(); ++i) {
		float currentlengthofpath = 0;
		for (unsigned int j = 0; j < PathVec[i].size(); ++j) {
			currentlengthofpath += PathVec[i][j][6];
		}
		lengthofchunks.push_back(currentlengthofpath);
	}

	// store length of chunk to every point in PathVec
	for (unsigned int i = 0; i < PathVec.size(); ++i) {
		for (unsigned int j = 0; j < PathVec[i].size(); ++j) {
			lengthofchunks_points.push_back(lengthofchunks[i]);
		}
	}

	// write to vtkfile
	//std::cout << "writing pathlengths to vtk file!\n";
	//PathBase::outputcontourvtkfile(this->chunklengths, this->vtkfilename, "pathlengths");
}

void GCodeAnalysis::setlayerheighttol(const float &layer_height_tol_val)
{
	this->layer_height_tol = layer_height_tol_val;
}

void GCodeAnalysis::calclayerheights(const PathBase::PathVector& PathVec)
{
	float avg_lay_height;
	for (auto& pathchunk : PathVec) {
		avg_lay_height = 0;
		for (std::size_t i = 0; i < pathchunk.size(); ++i) {
			avg_lay_height += pathchunk[i][2];
		}
		avg_lay_height = avg_lay_height / pathchunk.size();
		this->layer_height.push_back(avg_lay_height);
	}
	
	// next we generate a vec with all unique heights
	this->heights_avail.push_back(this->layer_height[0]); // first entry is unique per definition
	for (unsigned int j = 1; j < this->layer_height.size(); ++j) {
		// this requires the chunks to be in a sorted manner!
		// tol is 1/1000 of layer dim
		if (!arefloatsequal(this->layer_height[j], this->heights_avail.back(), this->layer_height_tol)) {
			this->heights_avail.push_back(this->layer_height[j]);
		}
	}
}


bool GCodeAnalysis::paircomparator(const std::pair<float, size_t> &l, const std::pair<float, size_t> &r)
{
	return (l.first > r.first);
}
