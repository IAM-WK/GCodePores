#include "GCodeParser.h"


GCodeParser::GCodeParser() noexcept
{
	// set some defaults?
}

void GCodeParser::setFilename(const std::string & Filename) noexcept
{
	this->GCodeFilename = Filename;
}

void GCodeParser::setSkirtoffset(const size_t & skirtoffsetval) noexcept
{
	this->skirtoffset = skirtoffsetval;
}

void GCodeParser::setBoundarystrings(const std::string & startstringval, const std::string & endstringval) noexcept
{
	this->startstring = startstringval;
	this->endstring = endstringval;
	if (!this->startstring.empty()) {
		this->startbystring = true;
	}
}

void GCodeParser::setCompatibility(const bool & freeformcompatval, const bool &slmcompatval) noexcept
{
	this->freeformcompat = freeformcompatval;
	this->slmcompat = slmcompatval;
}

void GCodeParser::readinfos(std::ifstream &GCodeFile, std::string &line_string, const std::string &coordchars)
{
	// count linenumber
	
	while (getline(GCodeFile, line_string))
	{
		if ( line_string.substr(0, 8) == "; skirts" || line_string.substr(0, 17) == ";   skirtOutlines") { // implement this variable where to search in line
			// could also search for various information here until everything needed is found or end is reached
			// then readinfo = false --> proceed to level 2 scan
			size_t skirtpos = line_string.find("= "); //position of skirt numbers in line_string
			if (skirtpos == std::string::npos) {
				skirtpos = line_string.find(",");
			}
			size_t endskirtpos = line_string.substr(skirtpos + 2, line_string.length()).find_first_not_of(coordchars);
			this->skirtnum = stoi(line_string.substr(skirtpos + 1, endskirtpos));
		}
		if (line_string.substr(0, 17) == "; extrusion_width" || line_string.substr(0, 17) == ";   extruderWidth") { // implement this variable where to search in line
			// determine extrusion width -> horizontal diameter = ex_w * pixel_width
			size_t exwpos = line_string.find("= ");
			if (exwpos == std::string::npos) {
				exwpos = line_string.find(",");
			}
			size_t endexwpos = line_string.substr(exwpos + 2, line_string.length()).find_first_not_of(coordchars);
			this->extrusion_width = stof(line_string.substr(exwpos + 1, endexwpos));
		}
		if (line_string.substr(0, 14) == "; layer_height" || line_string.substr(0, 15) == ";   layerHeight") { // implement this variable where to search in line
			// determine layer height -> vertical diameter = lay_h * pixel_height
			size_t layheightpos = line_string.find("= ");
			if (layheightpos == std::string::npos) {
				layheightpos = line_string.find(",");
			}
			size_t endlayheightpos = line_string.substr(layheightpos + 2, line_string.length()).find_first_not_of(coordchars);
			this->layer_height = stof(line_string.substr(layheightpos + 1, endlayheightpos));
		}
		++this->numoflines; // count number of lines
	}
	GCodeFile.clear(); // clear EOF flag
	GCodeFile.seekg(0, std::ios::beg); // reset to start
	return;
}



void GCodeParser::Execute(){
	
	//** define helper strings to parse GCode file. Definition of strings should depend on GCode flavour
	std::string comment, extrude, drive, drivetrav, feedratemask, BEDmask, powermask;
	if (!(this->freeformcompat) && !(this->slmcompat)) { // search patterns for slic3r, simplify3d....
		comment = ";", extrude = "E", drive = "G1", drivetrav = "G0", feedratemask = "F";
	}
	else if (!(this->slmcompat)) { // use search patterns for freeformer slicer
		comment = "//", extrude = "M", drive = "G01", drivetrav = "G00", feedratemask = "F";
	}
	else { // use search patterns for ORLAS slicer
		comment = "<", extrude = "M45", drive = "G01", drivetrav = "G0", feedratemask = "F", BEDmask = "G607", powermask = "G600";
		// new: M45 is persistent until M46 is reached x
		// G01 is denoted after Nxxxxx
		// startstring <Code> endstring </Code> (+tabs)
		// how to remove support structures?
	}

	std::ifstream GCodeFile;
	GCodeFile.exceptions(std::ifstream::badbit); // failbit seems to make conflict with reading last line...
	try {
		GCodeFile.open(this->GCodeFilename);
		// workaround for exception
		if (!GCodeFile.is_open()) {
			std::cerr << "Exception opening/reading GCode file!" << std::endl;
		}
	}
	catch (const std::ifstream::failure& e) {
		std::cerr << "Exception opening/reading GCode file: " << e.what() << std::endl;
		throw; // to throw it up to higher level -> should catch this in main and return EXIT_Failure;
		//throw std::runtime_error("problem reading GCode file!"); //
	}

	std::string line_string;
	const std::string coordchars = "0123456789.+-"; //allowed characters in a coordinate specification // +- for robustness

	float xcoord, ycoord, zcoord, dx, dy, dz, disttonext, feedrate = 0, beamexpandervalue = 0, laserpowervalue = 0;
	bool firstlines = true; // true while skipping lines before object
	bool isatravel = true; // printing head does not extrude during this movement
	bool startednewchunk = false; // movement before this one was a travel move
	size_t travelcount = 0, xpos = std::string::npos, ypos = std::string::npos, zpos = std::string::npos, extrudepos = std::string::npos, commentpos = std::string::npos, feedratepos = std::string::npos, beamexpanderpos = std::string::npos, powerpos = std::string::npos, slmnpos = std::string::npos, endslmnpos = std::string::npos;
	std::array<float, dim_coords> origin; origin.fill(0);

	//************** read infos from GCodeFile ahead ****************
	// ignore setup and skirt --> count travels; object starts n(skirts)+1 travels after begin 
	// n skirts is stated at the bottom of the file
	readinfos(GCodeFile, line_string, coordchars);
	this->PrintVector.reserve(this->numoflines); // overestimated guess of size to avoid reallocating
	//***************************************************************
	
	while (getline(GCodeFile, line_string))
	{
		//*************** analyse line ***************************************
		// tries to find coordinate data
		commentpos = line_string.find(comment);
		extrudepos = line_string.find(extrude); //position of extrusion coordinate in the line_string, needed to avoid jumps
		xpos = line_string.find("X"); //position of X coordinate in the line_string
		ypos = line_string.find("Y");
		zpos = line_string.find("Z");
		feedratepos = line_string.find(feedratemask);
		if (this->slmcompat) {
			beamexpanderpos = line_string.find(BEDmask);
			powerpos = line_string.find(powermask);
		}
		//********************************************************************

		//*************** check if line is only a travelmove *****************
		// if found any x/y/z and no e --> travel    (+ xyz have to occur before comment sign)
		if ((line_string.substr(0, drive.length()) == drive || line_string.substr(0, drivetrav.length()) == drivetrav) && (commentpos >= xpos) && (commentpos >= ypos) && (commentpos >= zpos)) {
			if (( xpos != std::string::npos || ypos != std::string::npos || zpos != std::string::npos ) && extrudepos == std::string::npos) {
				isatravel = true;
				++travelcount;
			}
			else {
				isatravel = false;
			}
		}
		// method to use with orlas slm compatibility, if above will always evaluate to false
		if (this->slmcompat) {
			// sections with laser on are started with "NXY M45" and ended with M46
			// since if statements aren't reached in lines withouth M45/6 current extrusion state is kept
			size_t extrudestoppos = line_string.find("M46");
			if (extrudepos != std::string::npos) {
				isatravel = false;
			}
			else if (extrudestoppos != std::string::npos) {
				isatravel = true;
			}

		}
		//********************************************************************

		//********** update printing position before start of object *********
		if (firstlines && (commentpos >= xpos) && (commentpos >= ypos) && (commentpos >= zpos)) {
			// start of object has not been reached, continue until found  && (travelcount <= (skirtnum + 6))
			if (xpos != std::string::npos) { //found x coord
				size_t endposx = line_string.substr(xpos + 1, line_string.length()).find_first_not_of(coordchars); //length of coord
				origin[0] = stof(line_string.substr(xpos + 1, endposx));
			}
			if (ypos != std::string::npos) { //found y coord
				size_t endposy = line_string.substr(ypos + 1, line_string.length()).find_first_not_of(coordchars); //length of coord
				origin[1] = stof(line_string.substr(ypos + 1, endposy));
			}
			if (zpos != std::string::npos) { //found z coord
				size_t endposz = line_string.substr(zpos + 1, line_string.length()).find_first_not_of(coordchars); //length of coord
				origin[2] = stof(line_string.substr(zpos + 1, endposz));
			}
			if (feedratepos != std::string::npos) { //found feedrate
				size_t endfeedratepos = line_string.substr(feedratepos + 1, line_string.length()).find_first_not_of(coordchars); //length of feedrate
				origin[7] = stof(line_string.substr(feedratepos + 1, endfeedratepos));
			}
			if (this->slmcompat && feedratepos != std::string::npos) { //found feedrate -- SLM
				size_t endfeedratepos = line_string.substr(feedratepos + 2, line_string.length()).find_first_not_of(coordchars); //length of feedrate
				origin[7] = stof(line_string.substr(feedratepos + 2, endfeedratepos));
			}
			if (this->slmcompat && beamexpanderpos != std::string::npos) { //found beam expander value
				size_t endbeamexpanderpos = line_string.substr(beamexpanderpos + 5, line_string.length()).find_first_not_of(coordchars); //length of beam expander value
				origin[8] = stof(line_string.substr(beamexpanderpos + 5, endbeamexpanderpos));
			}
			if (this->slmcompat && powerpos != std::string::npos) { //found laser power value
				size_t endpowerpos = line_string.substr(powerpos + 5, line_string.length()).find_first_not_of(coordchars); //length of laser power value
				origin[9] = stof(line_string.substr(powerpos + 5, endpowerpos));
			}
		}
		//********** start of object detection *******************************
		// check which start detection method has been set: 
		if (firstlines && !startbystring) {
			// this is the "offsetfromskirt" technique; techniques have to switch off the bool firstlines and care to push_back an origin
			// boundariesbyskirtoffset(bool &firstlines, const size_t &travelcount, const size_t &skirtoffset, const std::array &origin)
			if (travelcount == (this->skirtnum + this->skirtoffset)) {
				//set traveldestination as start coordinates
				std::cout << "pushed back origin [startbyskirtoffset]: " << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;
				this->PrintVector.push_back({ origin }); //"origin" contains starts coords + start feedrate
				firstlines = false;
				feedrate = origin[7]; beamexpandervalue = origin[8]; laserpowervalue = origin[9];
				continue;
			}
			else {
				continue; // otherwise crash because printvec gets accessed to early, or last travel is not saved as origin
			}
		}
		else if (firstlines) { // use startbystring
			if ( (line_string.substr(0, this->startstring.length()) == this->startstring) ) {
				//set traveldestination as start coordinates
				std::cout << "pushed back origin [startbystring]: " << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;
				this->PrintVector.push_back({ origin }); //"origin"
				firstlines = false;
				feedrate = origin[7]; beamexpandervalue = origin[8]; laserpowervalue = origin[9];
				continue;
			}
			else {
				continue; // otherwise crash because printvec gets accessed to early, or last travel is not saved as origin
			}
		}
		//********************************************************************

		//*********** end of object detection ********************************
		// ignore commands after end of printing code ----> needs improvement
		if (line_string.substr(0, this->endstring.length()) == this->endstring) {  //|| line_string.substr(0, 4) == ";End"
			break;
		}
		//********************************************************************

		//********* update feedrate/BED from commands without printing ***********
		if ( (line_string.substr(0, drive.length()) == drive || line_string.substr(0, drivetrav.length()) == drivetrav) && (commentpos >= feedratepos) && (feedratepos != std::string::npos) ) {
			size_t endfeedratepos = line_string.substr(feedratepos + 1, line_string.length()).find_first_not_of(coordchars); //length of feedrate
			feedrate = stof(line_string.substr(feedratepos + 1, endfeedratepos));
		}
		// SLM - feed: method to use with orlas slm compatibility, if above will always evaluate to false
		if (this->slmcompat && (commentpos >= feedratepos) && (feedratepos != std::string::npos)) {
			size_t endfeedratepos = line_string.substr(feedratepos + 2, line_string.length()).find_first_not_of(coordchars); //length of feedrate
			feedrate = stof(line_string.substr(feedratepos + 2, endfeedratepos));
		}
		// SLM - BED: method to use with orlas slm compatibility, if above will always evaluate to false
		if (this->slmcompat && (commentpos >= beamexpanderpos) && (beamexpanderpos != std::string::npos)) {
			size_t endbeamexpanderpos = line_string.substr(beamexpanderpos + 5, line_string.length()).find_first_not_of(coordchars); //length of BED
			beamexpandervalue = stof(line_string.substr(beamexpanderpos + 5, endbeamexpanderpos));
		}
		// SLM - laser power: method to use with orlas slm compatibility, if above will always evaluate to false
		if (this->slmcompat && (commentpos >= powerpos) && (powerpos != std::string::npos)) {
			size_t endpowerpos = line_string.substr(powerpos + 5, line_string.length()).find_first_not_of(coordchars); //length of BED
			laserpowervalue = stof(line_string.substr(powerpos + 5, endpowerpos));
		}
		//********************************************************************
		// helper position for slmcompat mode
		slmnpos = line_string.find("N");
		endslmnpos = line_string.substr(slmnpos + 1, line_string.length()).find_first_not_of(coordchars);

		//*************** construct coordinate vector ************************
		// filter lines with G1 command containing coordinates, containing at least one coordinate
		if ( (  (line_string.substr(0, drive.length()) == drive || line_string.substr(0, drivetrav.length()) == drivetrav) && (xpos != std::string::npos || ypos != std::string::npos || zpos != std::string::npos)  ) 
			|| (this->slmcompat && slmnpos != std::string::npos && (line_string.substr(endslmnpos+2, drive.length()) == drive || line_string.substr(endslmnpos+2, drivetrav.length()) == drivetrav) && (xpos != std::string::npos || ypos != std::string::npos || zpos != std::string::npos) )) {
			// extract coordinate values
			if (xpos != std::string::npos) { //found x coord
				size_t endposx = line_string.substr(xpos + 1, line_string.length()).find_first_not_of(coordchars); //length of coord
				xcoord = stof(line_string.substr(xpos + 1, endposx));
			}
			else {
				xcoord = this->PrintVector.back().back()[0]; //did not find x coord; use last one
			}
			if (ypos != std::string::npos) {
				size_t endposy = line_string.substr(ypos + 1, line_string.length()).find_first_not_of(coordchars);
				ycoord = stof(line_string.substr(ypos + 1, endposy));
			}
			else {
				ycoord = this->PrintVector.back().back()[1];
			}
			if (zpos != std::string::npos) {
				size_t endposz = line_string.substr(zpos + 1, line_string.length()).find_first_not_of(coordchars);
				zcoord = stof(line_string.substr(zpos + 1, endposz));
			}
			else {
				zcoord = this->PrintVector.back().back()[2];
			}
			// feedrate is determined above
			if (isatravel) { // isatravel so start new vec and and push back the travel end coords
				//if ((abs(xcoord - PrintVector.back().back()[0]) < tol) && (abs(ycoord - PrintVector.back().back()[1]) < tol) && (abs(zcoord - PrintVector.back().back()[2]) < tol)) { // coordinates need to differ from last coords, else just skip	
				if ( arefloatsequal(xcoord, this->PrintVector.back().back()[0]) && arefloatsequal(ycoord, this->PrintVector.back().back()[1]) && arefloatsequal(zcoord, this->PrintVector.back().back()[2]) ) {
					// do nothing when coords are not changing from last one
				}
				else if (startednewchunk) {
					this->PrintVector.back().back() = { xcoord,ycoord,zcoord,0,0,0,0, feedrate, beamexpandervalue, laserpowervalue } ;
				}
				else { // if printer has just travelled, this point will show wrong movement speed from travel move this->PrintVector.back().back()[7]
					this->PrintVector.push_back({{ xcoord,ycoord,zcoord,0,0,0,0, feedrate, beamexpandervalue, laserpowervalue } });
					startednewchunk = true;
				}
				//linenum = 1; // reset linenum with starting a new Vector
			}
			else {
				startednewchunk = false;
				if (arefloatsequal(xcoord, this->PrintVector.back().back()[0]) && arefloatsequal(ycoord, this->PrintVector.back().back()[1]) && arefloatsequal(zcoord, this->PrintVector.back().back()[2])) { 
					// coordinates need to differ from last coords, else just skip	
				}
				else {
					this->PrintVector.back().push_back({ xcoord,ycoord,zcoord,0,0,0,0,feedrate, beamexpandervalue, laserpowervalue });
				}
				
				//++linenum; //added entry in this path
			}
		}

	}

	GCodeFile.close();
	
	//************** calculate directions + disttonext *********************************
	for (auto &pathitem : this->PrintVector) {
		for (std::size_t i = 0; i < pathitem.size() - 1; ++i) {
			dx = pathitem[i + 1][0] - pathitem[i][0];
			dy = pathitem[i + 1][1] - pathitem[i][1];
			dz = pathitem[i + 1][2] - pathitem[i][2];
			
			disttonext = sqrt(dx * dx + dy * dy + dz * dz);
			dx = dx / disttonext; dy = dy / disttonext; dz = dz / disttonext; //norm direction
			pathitem[i][3] = dx; pathitem[i][4] = dy; pathitem[i][5] = dz; pathitem[i][6] = disttonext;
			// calulate printpathlength
			this->printpathlength += disttonext;
		}
	}
	for (std::size_t i = 0; i < this->PrintVector.size(); ++i) { // update last elements
		/*dx = this->PrintVector[i + 1][0][0] - this->PrintVector[i].back()[0];
		dy = this->PrintVector[i + 1][0][1] - this->PrintVector[i].back()[1];
		dz = this->PrintVector[i + 1][0][2] - this->PrintVector[i].back()[2];*/
		//if ((dx == 0.0 && dy == 0.0 && dz == 0.0) && (i != 0)) { // superfluous coord definition in GCode file
		//	dx = this->PrintVector[i-1].back()[3]; dy = this->PrintVector[i - 1].back()[4]; dz = this->PrintVector[i - 1].back()[5]; disttonext = this->PrintVector[i - 1].back()[6];
		//}
		if (this->PrintVector[i].size() >= 2) {
			dx = this->PrintVector[i][this->PrintVector[i].size() - 2][3];
			dy = this->PrintVector[i][this->PrintVector[i].size() - 2][4]; //caution with double travels !!!!!!
			dz = this->PrintVector[i][this->PrintVector[i].size() - 2][5];
		}
		else { // take last direction of -1 -1 -1 as fallback ----> Todo!
			dx = -1;
			dy = -1; 
			dz = -1;
		}
		// disttonext is 0 by definition
		this->PrintVector[i].back()[3] = dx; this->PrintVector[i].back()[4] = dy; this->PrintVector[i].back()[5] = dz; this->PrintVector[i].back()[6] = 0;
	}	
}

PathBase::PathVector GCodeParser::GetOutput() const
{
	if (this->PrintVector.size() == 0) {
		throw std::runtime_error("Error after parsing: No path was found!");
	}
	return this->PrintVector;
}

void GCodeParser::getGCodeparams(float & extr_width, float & lay_height) const
{
	extr_width = this->extrusion_width;
	lay_height = this->layer_height;
	return;
}

