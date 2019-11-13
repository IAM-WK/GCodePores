#include "PoreParser.h"

#include <iostream>

PoreParser::PoreParser()
{
	// set some defaults?
}

void PoreParser::setFilename(const std::string & Filename)
{
	porefilename = Filename;
}

void PoreParser::setSpherthres(const float & sphthres)
{
	sphericitythres = sphthres;
}

void PoreParser::setVolthres(const unsigned int & volthres)
{
	volumethres = volthres;
}

void PoreParser::setARthres(const float & ARthres)
{
	aspectratiothres = ARthres;
}

void PoreParser::Execute()
{
	// read in MorphoLibJ Output and store data
	// a selection is needed to clean out meaningless results
	//		eg: check if volume > volthreshold; R1->R3 > radthres; otherwise do not store this row
	std::ifstream PoreFile;
	std::string line_string;
	PoreFile.open(this->porefilename);
	
	const std::string columndelim = ",";
	// get fileheader
	getline(PoreFile, line_string);
	std::size_t delimpos = 0, endpos = std::string::npos;
	unsigned int columnum = 0;
	do {
		endpos = line_string.substr(delimpos + 1, line_string.length()).find(columndelim); // length header
		if (endpos == std::string::npos) { endpos = line_string.length(); }
		tupleheader[columnum] = line_string.substr(delimpos + 1, endpos);
		++columnum;
		delimpos = line_string.find(columndelim, delimpos + 1); // next start delim
	} while (delimpos != std::string::npos && (endpos != line_string.length() || endpos == std::string::npos ));
	// end header reading
	// todo: column numbers are not hardcoded, but determined by header
	// find col numbers of sphericity, volume, x,y,z, azimuth
	for (unsigned int colnum = 0; colnum < dimtuple; ++colnum) {
		if (tupleheader[colnum] == "Sphericity") {
			sphericitycol = colnum;
		}
		if (tupleheader[colnum] == "Volume") {
			volumecol = colnum;
		}
		if (tupleheader[colnum] == "Elli.Center.X") {
			Xcol = colnum;
		}
		if (tupleheader[colnum] == "Elli.Center.Y") {
			Ycol = colnum;
		}
		if (tupleheader[colnum] == "Elli.Center.Z") {
			Zcol = colnum;
		}
		if (tupleheader[colnum] == "Elli.Azim") {
			azimuthcol = colnum;
		}
		if (tupleheader[colnum] == "Elli.Elev") {
			elevationcol = colnum;
		}
		if (tupleheader[colnum] == "Elli.R1/R2") {
			rad12col = colnum;
		}
		if (tupleheader[colnum] == "Elli.R1/R3") {
			rad13col = colnum;
		}

	}
	//std::cout << "sphcol = " << sphericitycol << " vol " << volumecol << "y col" << Ycol << "azimuthcol " << azimuthcol << std::endl;
	float cur_poreazim; //azimut
	// tuple< azimut, volumeinvoxels>
	std::vector<std::tuple<float, int>> poreazimuth_overall;
	std::vector<std::tuple<float, int>> poreazimuth_filtered;

	float cur_poreelev; //elevation
	std::vector<std::tuple<float, int>> poreelevation_overall;
	std::vector<std::tuple<float, int>> poreelevation_filtered;


	std::array< float, PoreParser::dimtuple> linedata;
	while (getline(PoreFile, line_string))
	{
		// parse line
		columnum = 0; // zero column num
		do {
			endpos = line_string.substr(delimpos + 1, line_string.length()).find(columndelim); // length header
			if (endpos == std::string::npos) { endpos = line_string.length(); }
			linedata[columnum] = stof(line_string.substr(delimpos + 1, endpos));
			++columnum;
			delimpos = line_string.find(columndelim, delimpos + 1); // next start delim
		} while (delimpos != std::string::npos && (endpos != line_string.length() || endpos == std::string::npos));


		// std::cout << linedata[3] << std::endl;
		// add to overall pore volume
		porevolume += static_cast<unsigned int>(linedata[volumecol]);
		// add pore azimuth to histo_overall
		// azimuth should be between 0 - 180, only one sense of direction
		cur_poreazim = linedata[azimuthcol];
		if (cur_poreazim < 0.f) {
			cur_poreazim = cur_poreazim + 180.f;
		}
		poreazimuth_overall.push_back(std::make_tuple(cur_poreazim, static_cast<int>(linedata[volumecol])));
		// same for elevation
		cur_poreelev = linedata[elevationcol];
		if (cur_poreelev < 0.f) {
			cur_poreelev = cur_poreelev + 180.f;
		}
		poreelevation_overall.push_back(std::make_tuple(cur_poreelev, static_cast<int>(linedata[volumecol])));


		// only respect pores with sphericity < threshold and volume >= threshold
		//std::cout << " this data is in "<< tupleheader[3] << "value "<< linedata[3] << std::endl;
		if ( (linedata[sphericitycol] <= sphericitythres) && (linedata[volumecol] >= volumethres) && (linedata[rad12col] >= aspectratiothres || linedata[rad13col] >= aspectratiothres) ) {
			// std::cout << "accepted pore" << linedata[3] << std::endl;
			porelist.push_back(linedata);
			filteredporesvolume += static_cast<unsigned int>(linedata[volumecol]);
			poreazimuth_filtered.push_back(std::make_tuple(cur_poreazim, static_cast<int>(linedata[volumecol])));
			poreelevation_filtered.push_back(std::make_tuple(cur_poreelev, static_cast<int>(linedata[volumecol])));
		}

		// classification: gas pore or kh pore?
		if ((linedata[sphericitycol] <= sphericitythres) && (linedata[rad12col] >= aspectratiothres || linedata[rad13col] >= aspectratiothres)) {
			++khpores;
		}
		else {
			++gaspores;
		}

	}
	// finished gathering data
	// fill pore orientation histograms
	// sort pore lists by angle value
	/*std::sort(poreazimuth_overall.begin(), poreazimuth_overall.end(), [](const float &left, const float &right) {return left < right;});
	std::sort(poreazimuth_filtered.begin(), poreazimuth_filtered.end(), [](const float &left, const float &right) {return left < right; });
	std::sort(poreelevation_overall.begin(), poreelevation_overall.end(), [](const float &left, const float &right) {return left < right; });
	std::sort(poreelevation_filtered.begin(), poreelevation_filtered.end(), [](const float &left, const float &right) {return left < right; });
	*/
	// sort pore lists by angle value
	std::sort(poreazimuth_overall.begin(), poreazimuth_overall.end(), [](const std::tuple<float, int> &left, const std::tuple<float, int> &right) {
		return std::get<0>(left) < std::get<0>(right);
	});
	std::sort(poreazimuth_filtered.begin(), poreazimuth_filtered.end(), [](const std::tuple<float, int> &left, const std::tuple<float, int> &right) {
		return std::get<0>(left) < std::get<0>(right);
	});
	std::sort(poreelevation_overall.begin(), poreelevation_overall.end(), [](const std::tuple<float, int> &left, const std::tuple<float, int> &right) {
		return std::get<0>(left) < std::get<0>(right);
	});
	std::sort(poreelevation_filtered.begin(), poreelevation_filtered.end(), [](const std::tuple<float, int> &left, const std::tuple<float, int> &right) {
		return std::get<0>(left) < std::get<0>(right);
	});





	unsigned int classnumber = 1, currentclassfrequency = 0, classwidth = 10, currentclassvolume = 0;
	for (unsigned int anglenum = 0; anglenum < poreazimuth_overall.size(); ++anglenum) {
		// overall azim
		if (std::get<0>(poreazimuth_overall[anglenum]) < (static_cast<float>(classnumber*classwidth)) ) {
			++currentclassfrequency; 
			currentclassvolume += std::get<1>(poreazimuth_overall[anglenum]);
		}
		else if (std::get<0>(poreazimuth_overall[anglenum]) >= (static_cast<float>(classnumber)*classwidth)) {
			poreorientationhisto_azim_overall.push_back(currentclassfrequency);
			++classnumber;
			currentclassfrequency = 1; // 1 was found just now
			poreazimuth_histo_overall_volume.push_back(currentclassvolume);
			currentclassvolume = std::get<1>(poreazimuth_overall[anglenum]);
		}
	}
	poreorientationhisto_azim_overall.push_back(currentclassfrequency);
	poreazimuth_histo_overall_volume.push_back(currentclassvolume);
	classnumber = 1, currentclassfrequency = 0, currentclassvolume = 0;
	for (unsigned int anglenum = 0; anglenum < poreazimuth_filtered.size(); ++anglenum) {
		// fitered azim 
		if (std::get<0>(poreazimuth_filtered[anglenum]) < (static_cast<float>(classnumber*classwidth))) {
			currentclassvolume += std::get<1>(poreazimuth_filtered[anglenum]);
			++currentclassfrequency;
		}
		else if (std::get<0>(poreazimuth_filtered[anglenum]) >= (static_cast<float>(classnumber)*classwidth)) {
			poreorientationhisto_azim_filtered.push_back(currentclassfrequency);
			++classnumber;
			currentclassfrequency = 1; // 1 was found just now
			poreazimuth_histo_filtered_volume.push_back(currentclassvolume);
			currentclassvolume = std::get<1>(poreazimuth_filtered[anglenum]);
		}
	}
	poreazimuth_histo_filtered_volume.push_back(currentclassvolume);
	poreorientationhisto_azim_filtered.push_back(currentclassfrequency);
	//+++//
	classnumber = 1, currentclassfrequency = 0, currentclassvolume = 0;
	for (unsigned int anglenum = 0; anglenum < poreelevation_overall.size(); ++anglenum) {
		// overall elev
		if (std::get<0>(poreelevation_overall[anglenum]) < (static_cast<float>(classnumber*classwidth))) {
			++currentclassfrequency;
			currentclassvolume += std::get<1>(poreelevation_overall[anglenum]);
		}
		else if (std::get<0>(poreelevation_overall[anglenum]) >= (static_cast<float>(classnumber)*classwidth)) {
			poreorientationhisto_elev_overall.push_back(currentclassfrequency);
			++classnumber;
			currentclassfrequency = 1; // 1 was found just now
			poreelevation_histo_overall_volume.push_back(currentclassvolume);
			currentclassvolume = std::get<1>(poreelevation_overall[anglenum]);
		}
	}
	poreorientationhisto_elev_overall.push_back(currentclassfrequency);
	poreelevation_histo_overall_volume.push_back(currentclassvolume);
	classnumber = 1, currentclassfrequency = 0, currentclassvolume = 0;
	for (unsigned int anglenum = 0; anglenum < poreelevation_filtered.size(); ++anglenum) {
		// overall elev
		if (std::get<0>(poreelevation_filtered[anglenum]) < (static_cast<float>(classnumber*classwidth))) {
			++currentclassfrequency;
			currentclassvolume += std::get<1>(poreelevation_filtered[anglenum]);
		}
		else if (std::get<0>(poreelevation_filtered[anglenum]) >= (static_cast<float>(classnumber)*classwidth)) {
			poreorientationhisto_elev_filtered.push_back(currentclassfrequency);
			++classnumber;
			currentclassfrequency = 1; // 1 was found just now
			poreelevation_histo_filtered_volume.push_back(currentclassvolume);
			currentclassvolume = std::get<1>(poreelevation_filtered[anglenum]);

		}
	}
	poreorientationhisto_elev_filtered.push_back(currentclassfrequency);
	poreelevation_histo_filtered_volume.push_back(currentclassvolume);
	//return;
}

std::vector<std::array<float, 26>> PoreParser::GetOutput() const
{
	if (this->porelist.size() == 0) {
		throw std::runtime_error("Error after parsing: No path was found!");
	}
	return this->porelist;
}

