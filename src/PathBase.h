#pragma once
#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <fstream>
#include <numeric> // accumulate
#include <algorithm> // count

// base class to provide basics for GCode and path processing

class PathBase
{
public:

	//X;Y;Z,dx,dy,dz,disttonext,feedrate
	static constexpr uint32_t dim_coords = 8;


	typedef std::vector< std::vector< std::array< float, dim_coords > > > PathVector;
	typedef std::vector< std::array< float, dim_coords > > ChunkVector;
	

	// this writes coordinates to a .vtk file
	static void outputvtkfile(const PathBase::PathVector &coordinateVec, const std::string &filename, const std::vector<std::string> &argList = std::vector<std::string> () );

	// this writes coordinates to a .vtk file
	template <typename T>
	static void outputcontourvtkfile(const std::vector<T> &fieldvec, const std::string &filename, const std::string &fielddenominator) {
	
		std::ofstream vtkfile;
		vtkfile.open(filename, std::ios_base::app);
		// print color scalars
		vtkfile << "\n";
		//point data can only be specified once
		if (!PathBase::contourfieldwritten) { // this is omitted if there were other fields written before
			PathBase::contourfieldwritten = true;
			vtkfile << "POINT_DATA " << fieldvec.size() << " \n";
		}
		// just a consistency check
		if (PathBase::numberofpathpoints != fieldvec.size()) {
			std::cerr << "WARNING! contourfield size doesn't match path size in field: " << fielddenominator << "!" << std::endl;
		}
		// typeid().name may give mangled names on not MSVC compilers 
		vtkfile << "SCALARS" << " " << fielddenominator << " " << typeid(T).name() <<" 1 " << "\n";
		vtkfile << "LOOKUP_TABLE default \n";
		for (size_t i = 0; i < fieldvec.size(); ++i) {
			vtkfile << fieldvec[i] << " \n";
		}

		vtkfile.flush();
		vtkfile.close();
	}

	// set name of csv file
	static void setCSVfilename(const std::string &filename);

	// write a csv file with PathCoords and a number of user defined fields (variadic function)

	// writeCSV with empty argument list to stop recursion
	static void writeCSV() {}

	// provide fieldvecs + fielname string
	template<typename T, typename... Args>
	static void writeCSV(const std::vector<T> &DataVec, const std::string &fieldname, Args... args) {
		// write DataVec to CSV and forward args to next run
		// first check if DataVec is nonempty
		if (!DataVec.empty()) {
			const std::string filename = PathBase::csvfilename;
			const std::string tempfilename = PathBase::csvfilename.substr(0, PathBase::csvfilename.length() - 4) + "_temp.csv"; // write changes to this file
			if (!PathBase::csvfilewritten) {
				PathBase::csvfilewritten = true;
				// 1st run: create file, write vec to it --> static private bool to check if file written
				std::ofstream CSVFile;
				CSVFile.open(filename);
				// write header for this vec
				CSVFile << fieldname << " ; \n";
				for (size_t i = 0; i < DataVec.size(); ++i) {
					CSVFile << DataVec[i] << " ; \n";
				}
			}
			else {
				// read line from old file, append data, write to new file
				std::ifstream CSVFile;
				CSVFile.open(filename);
				std::string line_string;
				std::ofstream CSVtempFile;
				CSVtempFile.open(tempfilename);
				std::size_t line_num = 0; // data starts in line 2
										  // check if field lengths match...
										  // append header
				getline(CSVFile, line_string);
				std::string header = line_string;
				CSVtempFile << line_string << fieldname << " ; \n";
				while (getline(CSVFile, line_string)) {
					// append datavec to string and write to new file
					if (line_num < DataVec.size()) {
						CSVtempFile << line_string << DataVec[line_num] << " ; \n";
						++line_num;
					}
					else { // only copy old one and insert field delimiter
						CSVtempFile << line_string <<  " ; \n";
					}
				}
				// if DataVec is bigger than previous fields file has to be extended
				// insert enough delimiters
				size_t fieldnum = std::count(header.begin(), header.end(), ';');
				std::string emptyfields(fieldnum, ';');
				while(line_num < DataVec.size()) {
					CSVtempFile << emptyfields << DataVec[line_num] << " ; \n";
					++line_num;
				}

				CSVtempFile.flush();
				CSVtempFile.close();
				CSVFile.close();
				// rename new file to old file to overwrite, delete beforehand
				std::remove(filename.c_str());
				if (std::rename(tempfilename.c_str(), filename.c_str())) {
					std::cerr << "Encountered error while saving .csv file!\n";
				}
			}
		}
		// appending inline needs whole rewrite of file -> new file, rewrite every line + new data. overwrite old file with this
		writeCSV(args...); // how does this terminate when arglist is empty?
	}

	
	// worker process that outputs data
	static void writerWorker(const PathBase::PathVector &pathVec, const std::string &vtkfilename, const std::vector<std::string> &arglist, const std::vector<int> &pathclassification, const std::vector<float> &chunklength) {
		std::cout << "writing vtk file!\n";
		PathBase::outputvtkfile(pathVec, vtkfilename, arglist);
		std::cout << "writing pathclassifications to vtk file!\n";
		PathBase::outputcontourvtkfile(pathclassification, vtkfilename, "pathclassification");
		std::cout << "writing pathlengths to vtk file!\n";
		PathBase::outputcontourvtkfile(chunklength, vtkfilename, "pathlengths");
	}




	// ensure PathBase stays virtual...
	virtual ~PathBase() = 0;

private:

	static bool contourfieldwritten;
	static size_t numberofpathpoints;

	static bool csvfilewritten;
	static std::string csvfilename;



};

