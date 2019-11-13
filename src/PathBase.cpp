#include "PathBase.h"

//declare static class vars for contour fields
bool PathBase::contourfieldwritten = false;
size_t PathBase::numberofpathpoints = 0; 
bool PathBase::csvfilewritten = false;
std::string PathBase::csvfilename = "";


void PathBase::outputvtkfile(const PathBase::PathVector &coordinateVec, const std::string &filename, const std::vector<std::string> &argList)
{
	// preparation:
	// count number of points  --  PathBase::numberofpathpoints
	
	std::vector<size_t> chunklength; // stores length of chunk i
	for (auto const &chunkitem : coordinateVec) {
		chunklength.push_back(chunkitem.size());
	}
	 
	PathBase::numberofpathpoints = std::accumulate(chunklength.begin(), chunklength.end(), static_cast<size_t>(0)); // WTF!, needed to tell template to work with a zero of type size_t...;

	std::ofstream vtkfile;
	// set buffer
	//const size_t buffersize = 16384 * 16384;
	//std::unique_ptr<char[]> buf(new char[buffersize]);
	//char* buffer = new char[buffersize];
	//setbuf(vtkfile., buffer);

	//vtkfile.rdbuf()->pubsetbuf(buffer, buffersize);
	
	vtkfile.open(filename);
	
	// print header
	vtkfile << "# vtk DataFile Version 4.2 # "; // the command line comment has to be in the same line as the version string to preserve compatibility
	for (size_t j = 0; j < argList.size(); ++j) {
		vtkfile << argList[j] << " ";
	}
	vtkfile << "\n";
	vtkfile << "GCodeCoordinates\n";
	vtkfile << "ASCII\n";
	vtkfile << "DATASET POLYDATA\n";
	vtkfile << "POINTS " << PathBase::numberofpathpoints << " float\n";
	for (auto const &chunkitem : coordinateVec) {
		for (std::size_t i = 0; i < chunkitem.size(); ++i) {
			vtkfile << chunkitem[i][0] << " " << chunkitem[i][1] << " " << chunkitem[i][2] << "\n";
		}
	}
	vtkfile << "\n";
	vtkfile << "LINES " << chunklength.size() << " " << (PathBase::numberofpathpoints + chunklength.size()) << "\n";
	size_t countcoords = 0;
	for (auto const &linelength : chunklength) {
		vtkfile << linelength << " ";
		for (size_t j = 0; j < linelength; ++j) {
			vtkfile << countcoords << " ";
			++countcoords;
		}
		vtkfile << "\n";
	}
	vtkfile.flush();
	vtkfile.close();
	std::cout << "successfully wrote vtk file! \n";
}

void PathBase::setCSVfilename(const std::string & filename)
{
	if (filename.empty()) {
		std::cerr << "provide filename for CSV file!" << std::endl;
		throw std::runtime_error("csv filename empty!");
	}
	if (filename.substr(filename.length() - 4, filename.length()) != ".csv") {
		PathBase::csvfilename = filename + ".csv";
	}
	else {
		PathBase::csvfilename = filename;
	}
}


// DTOR
PathBase::~PathBase() {}
