/************************************************************************************
*************************************************************************************
 *  GCodePores
 *  analyse µCT porosity data in AM samples with respect to toolpath data
 *
 *  Copyright 2019 Lukas Englert
 *
 *
 *  Institute of Materials Science and Engineering, <http://www.iam.kit.edu/wk/english/index.php>
 *  Karlsruhe Institute of Technology
 *
 *  If you intend to use this work for your scientific publication  please cite the
 *  appropriate publications listed on
 *  <https://github.com/IAM-WK/GCodePores>
 *
 *************************************************************************************
 *************************************************************************************/

#include <limits>
#include <future>
#include <omp.h>

// own classes and functions
#include "GCodeParser.h"
#include "GCodeTransform.h"
#include "GCodeInterpolation.h"
#include "GCodeAnalysis.h"
#include "helperfunctions.h"
#include "PoreParser.h"
#include "mhdParser.h"
#include "TimeMeasurement.h"

// boost libs
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
	std::ios_base::sync_with_stdio(false);// said to improve IO perf...

	//***********************************************************************
	// INPUT 
	//***********************************************************************

	std::string GCodeFilename, PoreFilename, ImageFilename, OutputFilename, imagetranslationstring, startstring, endstring, vtkfilename;
	float interpolationdist, anglealpha, anglebeta, anglegamma, classwidth, sphericitythreshold, aspectratiothreshold, pathlengthclasswidth;
	std::vector<float> imagetranslation, imageangles;
	unsigned int nextneighbours, layersconsidered, volumethreshold;
	bool vtkfileonly = false, interpolationoff = false, vtkoutputoff = false, freeformcompat = false, slmcompat = false;
	std::size_t skirtoffset;

	try {
		// Declare the supported options
		po::options_description general("general");
		general.add_options()
			("help,h", "produce help message")
			("ImageFilename,i", po::value<std::string>(&ImageFilename)->required(), "specify Image filename - .mhd ONLY!")
			("GCodeFilename,g", po::value<std::string>(&GCodeFilename)->required(), "specify GCode filename")
			("PoreFilename,p", po::value<std::string>(&PoreFilename)->required(), "specify filename for file generated from MorphoLibJ")
			("OutputFilename,o", po::value<std::string>(&OutputFilename), "specify Output filename")
			("imagetranslation,u", po::value< std::string >(&imagetranslationstring)->required(), "Vector to translate GCode to fit in CT image (as string): \"x y z\" [mm]")
			("imagerotation,a", po::value< std::vector<float> >(&imageangles)->multitoken()->required(), "angle about z-axis between object in GCode and in CT-image [optional: provide angle about y- and x-axis; order: z y x] \nCAUTION: Object should still be approximately in the z-plane!")
			("interpolationdistance,d", po::value<float>(&interpolationdist), "specify interpolation distance, not required if --interpolationoff is specified")
			("neighbours,k", po::value<unsigned int>(&nextneighbours)->required(), "specify number of next neighbours to be considered")
			;
		po::options_description optional("optional");
		optional.add_options()
			("interpolationoff", "do not interpolate coordinates from GCode, just use points from GCode directly")
			("vtkfileonly", "only write the .vtk file containing the extracted printing coordinates as a line")
			("vtkoutputoff", "do not output a .vtk file")
			("startbystring,s", po::value<std::string>(&startstring), "if specified GCode will be extracted starting from this string (has to be unique!)")
			("startbyskirtoffset,j", po::value<size_t>(&skirtoffset)->default_value(7), "object starts numskirts+offset travels from start [default = 7]")
			("endstring,q", po::value<std::string>(&endstring)->default_value(";End"), "string until which GCode will be extracted [default = ;End]")
			("vtkfilename", po::value<std::string>(&vtkfilename)->default_value("coordinatefile.vtk"), "specify a filename for .vtk file")
			("freeform", "process GCode file in freeformer compatibility mode")
			("slm", "process GCode file in ORLAS SLM compatibility mode")
			("classwidth,c", po::value<float>(&classwidth)->default_value(10.f), "specify class width for output histogram (degrees)")
			("layers,l", po::value<unsigned int>(&layersconsidered)->default_value(2), "specify layers to be considered as neighbours")
			("sphthres,t", po::value<float>(&sphericitythreshold)->default_value(0.85f), "specify threshold; only pores with sphericity lower than t will be considered; use value >>1.0 to ignore")
			("arthres,r", po::value<float>(&aspectratiothreshold)->default_value(1.3f), "specify threshold; only pores with aspectratio R1/R2;;R1/R3 higher than r will be considered; use value <1.0 to ignore")
			("volumethres,v", po::value<unsigned int>(&volumethreshold)->default_value(27), "specify threshold; only pores with volume higher than v will be considered; use value of 0 to ignore")
			("pathlengthclasswidth,w", po::value<float>(&pathlengthclasswidth)->default_value(1.0f), "specify class width for pathlengths (mm)")
			; //end add options
		po::options_description all("Allowed options");
		all.add(general).add(optional);
		po::variables_map vm; // declare vm object
		po::store(po::parse_command_line(argc, argv, all), vm); // let vm contain options from command line

		if (vm.count("help")) {
			std::cout << all << "\n";
			std::cout << "--startbyskirtoffset will be used by default, use this option to override offset. When --startbystring is supplied, this method will be used." << "\n";
			std::cout << "sample invocation: \n";
			std::cout << "GCodePores.exe -i pacman_bw.mhd -g Program.xml -d 0.1 -u \"5.2 5.2 12.75\" -a 1.57079 -o porefile -k 4 -p pacman_pores_morpho.csv --slm --startbystring \"N1 G90\" -c 10 -t 0.7 -l 0 \n";
			return 1;
		}
		if (!(vm.count("interpolationoff")) && !(vm.count("interpolationdistance"))) {
			std::cerr << "please provide either interpolation distance or switch off interpolation by using --interpolationoff! \n";
			throw po::required_option("interpolationdistance");
			// throw exception that opt is missing
		}
		if ( !(vm.count("vtkfileonly")) && !(vm.count("OutputFilename")) ) {
			std::cerr << "please provide either output .csv filename or supply --vtkfileonly! \n";
			throw po::required_option("OutputFilename");
			// throw exception that opt is missing
		}
		po::notify(vm); // run notify after dealing with help, otherwise required options will throw exception
		
		if (!vm["imagerotation"].empty() && (imageangles = vm["imagerotation"].as<std::vector<float> >()).size() == 1) {
			anglegamma = imageangles[0];
			anglebeta = 0;
			anglealpha = 0;
		}
		else if (!vm["imagerotation"].empty() && (imageangles = vm["imagerotation"].as<std::vector<float> >()).size() == 2) {
			anglegamma = imageangles[0];
			anglebeta = imageangles[1];
			anglealpha = 0;
		}
		else if (!vm["imagerotation"].empty() && (imageangles = vm["imagerotation"].as<std::vector<float> >()).size() == 3) {
			anglegamma = imageangles[0];
			anglebeta = imageangles[1];
			anglealpha = imageangles[2];
		}
		else {
			std::cerr << "invalid value(s) for imagerotation! \n";
			return EXIT_FAILURE;
		}

		size_t start = 0, end = 0, tokens = 0; // token counting is necessary to catch last value
		while ((end = imagetranslationstring.find(' ', start)) != std::string::npos || (tokens < 3)) {
			imagetranslation.push_back(stof(imagetranslationstring.substr(start, end - start)));
			start = end + 1;
			++tokens;
		}
		if (vm.count("vtkfileonly")) { vtkfileonly = true; }
		if (vm.count("interpolationoff")) { interpolationoff = true; }
		if (vm.count("vtkoutputoff")) { vtkoutputoff = true; }
		if (vm.count("freeform")) { freeformcompat = true; }
		if (vm.count("slm")) { slmcompat = true; }
		if (vm.count("freeform") && !vm.count("startbystring")) { std::cerr << "freeform compatibility mode only works with startbystring! \n"; return EXIT_FAILURE; }
		

		if (!(ImageFilename.substr(ImageFilename.length() - 4, ImageFilename.length()) == ".mhd")) {
			ImageFilename = ImageFilename + ".mhd";
		}
		if (vtkoutputoff && (vtkfilename != "coordinatefile.vtk")) {
			std::cerr << "WARNING: You specified a vtk filename but disabled vtk output!" << std::endl;
		}

	}
	catch (std::exception& e)
	{
		std::cerr << "Error: " << e.what() << "\n";
		return EXIT_FAILURE;
	}

	// store command line arguments for documentation
	std::vector<std::string> argList(argv + 1, argv + argc);

	//***********************************************************************
	// END OF INPUT 
	//***********************************************************************

	//***********************************************************************
	// GCODE PREPROCESSING : PARSING - TRANSFORMATION - INTERPOLATION - ANALYSIS
	//***********************************************************************

	//*************** read GCode and generate object ************************
	const double wall_start = get_wall_time();
	GCodeParser PathParser;
	PathParser.setFilename(GCodeFilename);
	PathParser.setSkirtoffset(skirtoffset);
	PathParser.setBoundarystrings(startstring, endstring);
	PathParser.setCompatibility(freeformcompat, slmcompat);
	PathParser.Execute();// read GCode, catch read error 
	const double wall_parser = get_wall_time();
	std::cout << "GCodeParser finished! Time needed = " << std::to_string(wall_parser - wall_start) << "\n";
	//************************* end read GCode *******************************



	//*************** parse .mhd file for dimensions ************************
	// image dimensions are needed to transform GCode coordinates to image coordinate system
	mhdParser ImageParser;
	ImageParser.setFilename(ImageFilename);
	ImageParser.Execute();
	//********************* end mhd parser **********************************



	//********************** align GCode to Image ***************************
	GCodeTransform GCodeTransformer;
	GCodeTransformer.setOrigin(imagetranslation);
	GCodeTransformer.setAngles(anglegamma, anglebeta, anglealpha);
	GCodeTransformer.setImagedimensions(ImageParser.imagedimensionmm[1], ImageParser.imagedimensionmm[2]);
	GCodeTransformer.setInput(PathParser.GetOutput());
	GCodeTransformer.Execute(); // start of printing in image; {X,Y,Z}, Phi (radians), length of image in y' direction, z' direction; 5.32f,8.10f,3.0f, 2.73
								// select GCode representation
	const double wall_transform = get_wall_time();
	std::cout << "Path coordinate transformation finished! Time needed = " << std::to_string(wall_transform - wall_parser) << "\n";
	//********************* end coordinate transformation *******************



	//************************* interpolate GCode ***************************
	GCodeInterpolation GCodeInterpolator;
	GCodeInterpolator.setDistance(interpolationdist);
	GCodeInterpolator.setInput(GCodeTransformer.GetOutput());
	PathBase::PathVector pathVec; // to store Vector with points, PathVector is typedefd in the virtual base class PathBase
	if (!interpolationoff) {
		GCodeInterpolator.Execute(); //
		pathVec = GCodeInterpolator.GetOutput(); // use enum 
	}
	else {
		pathVec = GCodeTransformer.GetOutput();
	}
	const double wall_interpol = get_wall_time();
	std::cout << "Path interpolation finished! Time needed = " << std::to_string(wall_interpol - wall_transform) << "\n";
	//************************* end interpolate GCode ***********************




	//************************* analyse GCode *********************************
	//****************** write VTK file for PARAVIEW **************************
	if (!(vtkfilename.substr(vtkfilename.length() - 4, vtkfilename.length()) == ".vtk")) {
		vtkfilename = vtkfilename + ".vtk";
	}


	GCodeAnalysis GCodeAnalyser;
	GCodeAnalyser.setvtkfilename(vtkfilename);
	//GCodeAnalyser.calcdeltaangle(pathVec, true, true, 0.4f); // calc angles at points and create vtk colorfield, degree, filter, filterlength
	float extr_width = 0.f, l_hgt = 0.f;
	PathParser.getGCodeparams(extr_width, l_hgt);
	if (slmcompat) { extr_width = 0.125f; } //todo: estimate from GCode
	else if (arefloatsequal(extr_width, 0.f)) { extr_width = 0.45f; };
	GCodeAnalyser.classifypoints(pathVec, 0.5f*extr_width); // maybe a bit too conservative to use width as tolerance...
	GCodeAnalyser.calcpathlength(pathVec);
	const float layer_height_tol = 0.001f; // 1 mue tolerance
	GCodeAnalyser.setlayerheighttol(layer_height_tol);
	GCodeAnalyser.calclayerheights(pathVec);
	GCodeAnalyser.findfeedrate(pathVec);

	const double wall_analysis = get_wall_time();
	std::cout << "Path analysis finished! Time needed = " << std::to_string(wall_analysis - wall_interpol) << "\n";

	std::cout << "\nstarting to write data... \n";
	
	//************************* end analyse GCode ****************************

	// print normed coords -- just for debugging
	// for (auto &pathitem : pathVec) {
	//	std::cout << "begin chunk" << std::endl;
	//	for (std::size_t i = 0; i < pathitem.size(); ++i) { // " -- dX: " << pathitem[i][3] << " -- dY: " << pathitem[i][4] << " -- dZ: " << pathitem[i][5] << " -- disttonext " << pathitem[i][6] << std::endl
	//		std::cout << "X: " << pathitem[i][0] << "Y: " << pathitem[i][1] << "Z: " << pathitem[i][2] << "feedrate: " << pathitem[i][7] << std::endl;
	//	}
	//}
	//************************* write data ****************************
	//** data writing consumes much time, so writing will be done asynchronous to the rest of the work **
	// async launches the worker which writes the fields to the vtk file
	std::future<void> writework;
	if (!vtkoutputoff) {
		writework = std::async(std::launch::async, PathBase::writerWorker, pathVec, vtkfilename, argList, GCodeAnalyser.pathclassification, GCodeAnalyser.lengthofchunks_points, GCodeAnalyser.feedrate);
	}
	if (vtkfileonly) {
		if (!vtkoutputoff) {
			writework.get();// this makes sure, writer finishes before exit
		}
		return EXIT_SUCCESS; // do not process further
	}

	//***********************************************************************
	// END OF GCODE PREPROCESSING
	//***********************************************************************
	
	//***********************************************************************
	// MLJ POREFILE PREPROCESSING: PARSING - ANALYSIS - TRANSFORMATION
	//***********************************************************************

	//************************* parse pore file from MLJ ********************

	PoreParser poreparser;
	poreparser.setFilename(PoreFilename);
	poreparser.setSpherthres(sphericitythreshold);
	poreparser.setVolthres(volumethreshold);
	poreparser.setARthres(aspectratiothreshold);
	//std::cout << PoreFilename << std::endl;
	poreparser.Execute();
	// get columns of needed values
	unsigned int volumecol = poreparser.volumecol, Xcol = poreparser.Xcol, Ycol = poreparser.Ycol, Zcol = poreparser.Zcol, azimuthcol = poreparser.azimuthcol;
	//  sphericitycol = poreparser.sphericitycol, TODO
	//************************* end parse pore file  ************************


	//*************** output processed information from MLJ *****************
	std::cout << "volume of all pores = " << poreparser.porevolume << " voxels; => " << (poreparser.porevolume*ImageParser.pixelwidth[0] * ImageParser.pixelwidth[1] * ImageParser.pixelwidth[2]) << "mm^3 \n";
	std::cout << "volume of considered pores = " << poreparser.filteredporesvolume << " voxels; => " << (poreparser.filteredporesvolume* ImageParser.pixelwidth[0] * ImageParser.pixelwidth[1] * ImageParser.pixelwidth[2]) << "mm^3 \n";
	std::cout << "this corresponds to a fraction of " << (std::roundf(static_cast<float>(poreparser.filteredporesvolume) / static_cast<float>(poreparser.porevolume) * 1000)/10) << "% \n";
	std::cout << "fraction of pores (by number) that were classified as non gas pores: " << (std::roundf(static_cast<float>(poreparser.khpores) / static_cast<float>(poreparser.khpores+poreparser.gaspores) * 1000) / 10) << "% \n";
	// save this informations to a vec to output it to the csv file
	std::vector<std::string> poreinformations;
	poreinformations.push_back("volume of all pores = " + std::to_string(poreparser.porevolume) + " voxels; => " + std::to_string(poreparser.porevolume* ImageParser.pixelwidth[0] * ImageParser.pixelwidth[1] * ImageParser.pixelwidth[2]) + "mm^3 ");
	poreinformations.push_back("volume of considered pores = " + std::to_string(poreparser.filteredporesvolume) + " voxels; => " + std::to_string(poreparser.filteredporesvolume* ImageParser.pixelwidth[0] * ImageParser.pixelwidth[1] * ImageParser.pixelwidth[2]) + "mm^3 ");
	poreinformations.push_back("this corresponds to a fraction of " + std::to_string(std::roundf(static_cast<float>(poreparser.filteredporesvolume) / static_cast<float>(poreparser.porevolume) * 1000) / 10) + "% ");
	poreinformations.push_back("fraction of pores (by number) that were classified as non gas pores: " + std::to_string(std::roundf(static_cast<float>(poreparser.khpores) / static_cast<float>(poreparser.khpores + poreparser.gaspores) * 1000) / 10) + "% ");
	//*************** end output processed information **********************

	std::vector < std::array<float, 26>> porevec;
	porevec = poreparser.GetOutput();
	// entry 5,6,7 is center x,y,z // 8,9,10 is rad1,2,3 // 
	// post process pore coords: pixels to mm
	for (unsigned int porenum = 0; porenum < porevec.size(); ++porenum) {
		porevec[porenum][Xcol] = porevec[porenum][Xcol] * ImageParser.pixelwidth[0];
		porevec[porenum][Ycol] = porevec[porenum][Ycol] * ImageParser.pixelwidth[1];
		porevec[porenum][Zcol] = porevec[porenum][Zcol] * ImageParser.pixelwidth[2];
	}

	//***********************************************************************
	// END MLJ POREFILE PREPROCESSING
	//***********************************************************************

	//********************************************************************************************************************************************
	//********************************************************************************************************************************************


	// float: deltaangle, int: class, int: volume
	std::vector<std::tuple<float, int, int>> deltaangles(porevec.size(), std::make_tuple(150.f, 0, 0));
	// vec to store how many pores a chunk contains; size is number of chunks, start with 0
	std::vector<size_t> poreinchunkfrequencies(pathVec.size(), 0);
	// vec to store how many defect voxels a chunk contains; size is number of chunks, start with 0
	std::vector<size_t> poreinchunkfrequencies_voxel_weighted(pathVec.size(), 0);
	// what is the highest class found?
	size_t maxclass = *max_element(GCodeAnalyser.pathclassificationwithchunks.begin(), GCodeAnalyser.pathclassificationwithchunks.end());
	// vec to store number of pores in each class; eg class 5 -> highest entry = 5
	std::vector<size_t> frequenciesinclasses(maxclass+1, 0);
	// vec to store number of pores in each class; eg class 5 -> highest entry = 5
	//TODO: unsigned int maybe to small
	std::vector<size_t> frequenciesinclasses_volumeweighted(maxclass + 1, 0);
	// vec to store how many pores a layer contains
	std::vector<size_t> poresinlayer(GCodeAnalyser.heights_avail.size(), 0);
	// same, but summed volume instead of number
	std::vector<size_t> porevolumeinlayer(GCodeAnalyser.heights_avail.size(), 0);

	//**************************************************************************

	float current_pore_height, minheightdiff, GCodeAzimuth, angledelta, PoreAzimuth;

	unsigned int chunknum_next;
	int deltaindex, height_index;

	// loop through each pore

	size_t progresscounter = 0; // for progressbar
	const size_t sizeporevec = porevec.size();


#pragma omp parallel for private(current_pore_height, minheightdiff, GCodeAzimuth, angledelta, PoreAzimuth, chunknum_next, deltaindex, height_index)
	for (int porenum = 0; porenum < porevec.size(); ++porenum) {

		// unfortunately this is not so easy in parallel ... :(
#pragma omp atomic
		++progresscounter;

#pragma omp critical (outpinfo)
		{
			// show and update progressbar
			progressbar(progresscounter, sizeporevec);
		}
		//********************************************************************************************************************************************
		// reduce points to be considered by sorting out relevant layers
		// first find next layer coord for porecenter.z
		current_pore_height = porevec[porenum][Zcol];
		chunknum_next = std::numeric_limits<unsigned int>::max();
		minheightdiff = std::numeric_limits<float>::max();

		// iterate through height of chunks (layer_height)
		// if height pore - height chunk -> minimal ---> this is the next layer resp. next chunk, but other chunks in same layer are also possible
		for (unsigned int j = 0; j < GCodeAnalyser.layer_height.size(); ++j) {
			if (abs(current_pore_height - GCodeAnalyser.layer_height[j]) < minheightdiff) {
				chunknum_next = j;
				minheightdiff = abs(current_pore_height - GCodeAnalyser.layer_height[j]);
			}
		}

		//***////***////***////***////***////***////***////***////***////***////***////***////***////***//
		// we also can now classify the pore as being in Perimeter or Hatching
		// increment frequency for class of chunk
		//

		// this pore is in layer with height of layer_height[chunknum_next]
		// search entry in heights_avail with this layer height -> index in heights_avail is index for poresinlayer
		for (unsigned int layernum = 0; layernum < GCodeAnalyser.heights_avail.size(); ++layernum) {
			if (arefloatsequal(GCodeAnalyser.layer_height[chunknum_next], GCodeAnalyser.heights_avail[layernum], layer_height_tol)) {
#pragma omp atomic
				poresinlayer[layernum] += 1;
#pragma omp critical (calcvolume) // atomic only supports normal incrementation
				{
					porevolumeinlayer[layernum] += static_cast<size_t>(porevec[porenum][volumecol]);
				}
			}
		}
		// next step is to reduce heights_avail to only the heights we need
		// formula to alternatingly select height +1/-1/+2/-2 ...:
		// deltaindex = static_cast<int>(round((static_cast<float>(deltaindex) * (-2.f) + 1.f) / 2.f));
		deltaindex = 0, height_index = 0;

		std::vector< float> heights_considered;
		for (unsigned int j = 0; j < GCodeAnalyser.heights_avail.size(); ++j) {
			if (arefloatsequal(GCodeAnalyser.heights_avail[j], GCodeAnalyser.layer_height[chunknum_next], layer_height_tol)) {
				heights_considered.push_back(GCodeAnalyser.heights_avail[j]); // nearest layer itself
				// loop with deltaindex to consider neighbouring layers
				for (unsigned int numdeltalayer = 0; numdeltalayer < layersconsidered; ++numdeltalayer) {
					deltaindex = static_cast<int>(round((static_cast<float>(deltaindex) * (-2.f) + 1.f) / 2.f));
					height_index = j + deltaindex;
					if ((height_index >= 0) && (height_index < GCodeAnalyser.heights_avail.size())) {
						heights_considered.push_back(GCodeAnalyser.heights_avail[height_index]);
					}
				}
				// if height was found, then there should not be anything else to be found...
				break;
			}
		}
		/*std::cout << " size avail heights =  " << heights_considered.size() << std::endl;
		for (unsigned int j = 0; j < heights_considered.size(); ++j) {
			std::cout << heights_considered[j] << std::endl;
		}*/
		// next we need to translate considered heights to chunknums
		// so we generate a vector with all chunknumbers to be considered
		// make a comparison with 
		std::vector<unsigned int> chunks_considered;
		for (unsigned int j = 0; j < GCodeAnalyser.layer_height.size(); ++j) {
			for (unsigned int m = 0; m < heights_considered.size(); ++m) {
				if (arefloatsequal(GCodeAnalyser.layer_height[j], heights_considered[m], layer_height_tol)) {
					chunks_considered.push_back(j);
				}
			}
		}
		//std::cout << "considchunk chunkssize =  " << chunks_considered.size() << std::endl;
		/*for (unsigned int j = 0; j < chunks_considered.size(); ++j) {
			std::cout << chunks_considered[j] << std::endl;
		}
*/
// end reduce neighbouring points
// chunks considered contains every chunk to be considered
//********************************************************************************************************************************************

//********************************************************************************************************************************************
// find nearest neighbours of pore
//**************************************************************
// stores X,Y,Z,dX,dY,dZ,distance to pore, pathclass, chunknumber -- of neighbourhood path points 
		std::vector<std::array<float, 9> > neighbourhood;

		// loop through remaining path points, calc distance and store to neighbourhood
		for (unsigned int chunknum_index = 0; chunknum_index < chunks_considered.size(); ++chunknum_index) {
			std::array<float, 9> current_neighbour;
			for (unsigned int point_index = 0; point_index < pathVec[chunks_considered[chunknum_index]].size(); ++point_index) {
				// transfer location and direction of point
				for (unsigned int a = 0; a < 6; ++a) {
					current_neighbour[a] = pathVec[chunks_considered[chunknum_index]][point_index][a];
				}
				// calculate distance pore to point
				current_neighbour[6] = sqrt(((porevec[porenum][Xcol] - pathVec[chunks_considered[chunknum_index]][point_index][0]) * (porevec[porenum][Xcol] - pathVec[chunks_considered[chunknum_index]][point_index][0])) + ((porevec[porenum][Ycol] - pathVec[chunks_considered[chunknum_index]][point_index][1]) * (porevec[porenum][Ycol] - pathVec[chunks_considered[chunknum_index]][point_index][1])) + ((porevec[porenum][Zcol] - pathVec[chunks_considered[chunknum_index]][point_index][2]) * (porevec[porenum][Zcol] - pathVec[chunks_considered[chunknum_index]][point_index][2])));
				// find out pathclassification for this neighbouring point
				current_neighbour[7] = static_cast<float>(GCodeAnalyser.pathclassificationwithchunks[chunks_considered[chunknum_index]]);
				// store chunknumber of this very point --> needed for porespermm
				current_neighbour[8] = static_cast<float>(chunks_considered[chunknum_index]);
				neighbourhood.push_back(current_neighbour);

			}
		}
		// sort neighbourhood to pore distance (lambda function)
		std::sort(neighbourhood.begin(), neighbourhood.end(), [](const std::array<float, 9> &a, const std::array<float, 9> &b) {
			return (a[6] < b[6]);
		});


		/*	for (int i = 0; i < 5; ++i) {
				std::cout << "distance to nearest point number " << i << " is = " << neighbourhood[i][6];
			}*/


			// this is erroneus, we need further classification before we know which one is right
			//pathclassificationwithchunks[chunknum_next] = wrong
			// nearest neighbours chunknumber is needed to find out pathclassification
			// this classification is stored to freqinclasses
			// neighbourhood[0][7] is nearest neighbours' pathclassification
#pragma omp atomic
		frequenciesinclasses[static_cast<size_t>(neighbourhood[0][7])] += 1;
#pragma omp critical (pathclassification_volweight) // atomic only supports normal incrementation
		{
			frequenciesinclasses_volumeweighted[static_cast<int>(neighbourhood[0][7])] += static_cast<size_t>(porevec[porenum][volumecol]);
		}
		// now we know the nearest chunk to the pore
		//***////***// this is enough to classify for pathlengths correlation, so do that now //***////***//
		// increment number of pores the nearest chunk contains
		// this vector contains the number of pores for each chunk

#pragma omp atomic
		//poreinchunkfrequencies[chunknum_next] += 1;
		poreinchunkfrequencies[static_cast<size_t>(neighbourhood[0][8])] += 1;
#pragma omp critical (poreinchunkfreq_volweighted) // atomic only supports normal incrementation
		{
			poreinchunkfrequencies_voxel_weighted[static_cast<size_t>(neighbourhood[0][8])] += static_cast<size_t>(porevec[porenum][volumecol]);
		}

		// calculate average direction of GCode in pore neighbourhood
		// to do: weight direction with distance
		std::array<float, 3 > neighbour_direction; neighbour_direction.fill(0);
		for (unsigned int neighbour_index = 0; neighbour_index < nextneighbours; ++neighbour_index) {
			for (size_t h = 0; h < 3; ++h) {
				neighbour_direction[h] = neighbour_direction[h] + neighbourhood[neighbour_index][h+3];
			}
		}
		for (int h = 0; h < 3; ++h) {
			neighbour_direction[h] = neighbour_direction[h] / static_cast<float>(nextneighbours);
		}

		// convert GCode direction to azimuth/elevation system
		// GCode will usually have elevation = 0°
		// so elevation of pores cannot be correlated with GCode
		// azimuth = 180 / pi * atan2(dy, dx);
		// azimuth of pore is entry 11;
		GCodeAzimuth = (180.f / static_cast<float>(PI)) * atan2(neighbour_direction[1], neighbour_direction[0]);

		// azimuth should be between 0 - 180, only one sense of direction
		if (GCodeAzimuth < 0.f) {
			GCodeAzimuth = GCodeAzimuth + 180.f;
		}
		
		// PoreAziumuth 0 - 180
		PoreAzimuth = porevec[porenum][azimuthcol];
		if (PoreAzimuth < 0.f) {
			PoreAzimuth = PoreAzimuth + 180.f;
		}



		// calculate delta angle between pore and gcode
		angledelta = abs( GCodeAzimuth - PoreAzimuth);
		// delta should be between 0 - 90 : angles between 90 - 180 should be converted to smaller angle
		if (angledelta > 90.f) {
			angledelta = 180.f - angledelta;
		}

		// store delta phi (classified or do that later on?)
		deltaangles[porenum] = std::make_tuple(angledelta, static_cast<int>(neighbourhood[0][7]), static_cast<int>(porevec[porenum][volumecol]) );
		//std::cout << "calced delta " << angledelta << std::endl;

	}
	std::cout << "\n";
	// now do a classification on deltaangles
	// 1. sort vec
	// 2. specify class width
	// loop
	//	3. count objects below i*class width
	//	4. if object exceeds i*class width-> ++i (next class)
	//	5. store in a array with 180/class width entries the number of objects for each i
	// 6. write this to a output file as: i*class width/2 || num angles ----> histogram
	std::cout << "klassenbreite = " << classwidth << std::endl;

	//std::sort(deltaangles.begin(), deltaangles.end());
	// sort deltaangles by angle value;
	std::sort(deltaangles.begin(), deltaangles.end(), [](const std::tuple<float, int, int> &left, const std::tuple<float, int, int> &right) {
		return std::get<0>(left) < std::get<0>(right);
	});

	// prepare classwidth vectors for deltaangles and pathlengthclassification
	std::vector<std::string> classwidthvecangles;
	std::vector<std::string> classwidthvecpathlength;

	/*
	for (int i = 0; i < 100; ++i) {
		std::cout << deltaangles[i] << std::endl;
	}
*/
	// histogram of all deltaangles = histogram
	// histogram of deltaangles in all perims = perimhistogram
	// histogram of deltaangles in all hatchings = hatchhistogram
	unsigned int classnumber = 1, currentclassfrequency = 0, currentclassvolume = 0, perimcurrentclassvolume = 0, hatchingcurrentclassvolume = 0, perimclassnumber = 1, perimcurrentclassfrequency = 0, hatchclassnumber = 1, hatchcurrentclassfrequency = 0;
	std::vector <unsigned int> histogram;
	std::vector <unsigned int> perimhistogram;
	std::vector <unsigned int> hatchhistogram;
	std::vector <unsigned int> histogram_vol_in_angle;
	std::vector <unsigned int> perimhistogram_vol_in_angle;
	std::vector <unsigned int> hatchinghistogram_vol_in_angle;



	for (unsigned int anglenum = 0; anglenum < deltaangles.size(); ++anglenum) {
		// overall 
		if (std::get<0>(deltaangles[anglenum]) < (static_cast<float>(classnumber)*classwidth)) {
			++currentclassfrequency;
			currentclassvolume += std::get<2>(deltaangles[anglenum]);
		}
		else if (std::get<0>(deltaangles[anglenum]) >= (static_cast<float>(classnumber)*classwidth)) {
			histogram.push_back(currentclassfrequency);
			histogram_vol_in_angle.push_back(currentclassvolume);
			++classnumber;
			currentclassfrequency = 1; // 1 was found just now
			currentclassvolume = std::get<2>(deltaangles[anglenum]);
		}
		// perim
		if ( (std::get<1>(deltaangles[anglenum]) != 0) && (std::get<0>(deltaangles[anglenum]) < (static_cast<float>(perimclassnumber)*classwidth)) ) {
			++perimcurrentclassfrequency;
			perimcurrentclassvolume += std::get<2>(deltaangles[anglenum]);
		}
		else if ( (std::get<1>(deltaangles[anglenum]) != 0) && (std::get<0>(deltaangles[anglenum]) >= (static_cast<float>(perimclassnumber)*classwidth)) ) {
			perimhistogram.push_back(perimcurrentclassfrequency);
			perimhistogram_vol_in_angle.push_back(perimcurrentclassvolume);
			++perimclassnumber;
			perimcurrentclassfrequency = 1; // 1 was found just now
			perimcurrentclassvolume = std::get<2>(deltaangles[anglenum]);
		}
		// hatching
		if ((std::get<1>(deltaangles[anglenum]) == 0) && (std::get<0>(deltaangles[anglenum]) < (static_cast<float>(hatchclassnumber)*classwidth))) {
			++hatchcurrentclassfrequency;
			hatchingcurrentclassvolume += std::get<2>(deltaangles[anglenum]);
		}
		else if ((std::get<1>(deltaangles[anglenum]) == 0) && (std::get<0>(deltaangles[anglenum]) >= (static_cast<float>(hatchclassnumber)*classwidth))) {
			hatchhistogram.push_back(hatchcurrentclassfrequency);
			hatchinghistogram_vol_in_angle.push_back(hatchingcurrentclassvolume);
			++hatchclassnumber;
			hatchcurrentclassfrequency = 1; // 1 was found just now
			hatchingcurrentclassvolume = std::get<2>(deltaangles[anglenum]);
		}
	}
	histogram.push_back(currentclassfrequency);
	perimhistogram.push_back(perimcurrentclassfrequency);
	hatchhistogram.push_back(hatchcurrentclassfrequency);
	histogram_vol_in_angle.push_back(currentclassvolume);
	hatchinghistogram_vol_in_angle.push_back(hatchingcurrentclassvolume);
	perimhistogram_vol_in_angle.push_back(perimcurrentclassvolume);
	std::vector<float> histogram_vol_in_mm_angle;
	histogram_vol_in_mm_angle.resize(histogram_vol_in_angle.size());
	
	//std::transform(histogram_vol_in_angle.begin(), histogram_vol_in_angle.end(), histogram_vol_in_mm_angle.begin(),std::bind(std::multiplies<float>(), std::placeholders::_1, (pixelwidth[0]*pixelwidth[1]*pixelwidth[2])));
	//for_each(histogram_vol_in_mm_angle.begin(), histogram_vol_in_mm_angle.end(), [ &pixelwidth[0]] ( float &el) {});

	//auto it_vox = histogram_vol_in_angle.begin();
	//for_each(histogram_vol_in_mm_angle.begin(), histogram_vol_in_mm_angle.end(), [&it_vox, &pixelwidth](auto &it_mm) { 
	//	it_mm = *it_vox++;
	//	*it_mm = it_vox*pixelwidth[0] * pixelwidth[1] * pixelwidth[2]; 
	//});
	
	// normal for loop...
	for (unsigned int index = 0; index < histogram_vol_in_mm_angle.size(); ++index) {
		histogram_vol_in_mm_angle[index] = histogram_vol_in_angle[index] * ImageParser.pixelwidth[0] * ImageParser.pixelwidth[1] * ImageParser.pixelwidth[2];
	}

	/*
	std::cout << "anzahl winkel " << deltaangles.size() << std::endl;
	std::cout << "histogram size" << histogram.size() << std::endl;*/
	std::cout << "klassenzahl = " << classnumber << std::endl;
	
	for (int i = 0; i < histogram.size(); ++i) {
		std::cout << classwidth*i << " - " << classwidth*(i+1) << " : " << histogram[i] << std::endl;
		classwidthvecangles.push_back(std::to_string(classwidth*i) + " - " + std::to_string(classwidth * (i + 1)));
	}

	//***// //***// //***// 
	// we have a vector with the frequencies of pores for each chunk
	// classify lengths of chunks, so that similar chunklengths go in the same class
	//  chunk i has lengthofchunks[i]; poreinchunkfrequencies gives num of pores in them
	// first: which chunks go in which class?
	// number of classes = max lengthofchunks / classwidth +1
	float pathlengthnumofclasses = (*(std::max_element(GCodeAnalyser.lengthofchunks.begin(), GCodeAnalyser.lengthofchunks.end())) / pathlengthclasswidth) + 1;
	// stores pores per mm for each pathlengthclass
	std::vector<float> pathlengthshisto (static_cast<size_t>(pathlengthnumofclasses), 0); 
	// stores pores per mm for each pathlengthclass - voxel weighted
	std::vector<float> pathlengthshisto_volumeweighted(static_cast<size_t>(pathlengthnumofclasses), 0);
	// stores length of all chunks in the same class
	std::vector<float> pathlengthsinclass(static_cast<size_t>(pathlengthnumofclasses), 0);



	// classify chunk in a pathlength; iterate through classes and take every chunk that fits in class
	for (unsigned int pathlengthclass = 0; pathlengthclass < pathlengthshisto.size(); ++pathlengthclass) {
		for (unsigned int chunknumber = 0; chunknumber < GCodeAnalyser.lengthofchunks.size(); ++chunknumber) {
			// check if this chunk fits in current class
			if ((GCodeAnalyser.lengthofchunks[chunknumber] >= static_cast<float>(pathlengthclass)*pathlengthclasswidth) && (GCodeAnalyser.lengthofchunks[chunknumber] < static_cast<float>(pathlengthclass+1)*pathlengthclasswidth) ) {
				// add pores to pathlengthshisto[pathlenghtsclass] and add length to pathlengthsinclass for normalising at the end
				pathlengthshisto[pathlengthclass] += static_cast<float>(poreinchunkfrequencies[chunknumber]);
				pathlengthsinclass[pathlengthclass] += GCodeAnalyser.lengthofchunks[chunknumber];
				pathlengthshisto_volumeweighted[pathlengthclass] += poreinchunkfrequencies_voxel_weighted[chunknumber];
			}
		}
		// finished iterating over all chunks, now normalise pores in chunk to pores per mm
		// make sure there was at least one pore in chunk class
		if (!arefloatsequal(pathlengthshisto[pathlengthclass], 0.f)) {
			pathlengthshisto[pathlengthclass] = pathlengthshisto[pathlengthclass] / pathlengthsinclass[pathlengthclass];
		}
		// make sure there was at least one pore in chunk class
		if (pathlengthshisto_volumeweighted[pathlengthclass] != 0) {
			pathlengthshisto_volumeweighted[pathlengthclass] = static_cast<float>(static_cast<float>(pathlengthshisto_volumeweighted[pathlengthclass]) / pathlengthsinclass[pathlengthclass]);
		}
	}

	// output the pores per mm as a histogram (sorted)
	std::cout << "Histogram: pores per mm over length of the paths" << std::endl;
	for (int i = 0; i < pathlengthshisto.size(); ++i) {
		std::cout << pathlengthclasswidth * i << " - " << pathlengthclasswidth * (i + 1) << " : " << pathlengthshisto[i] << std::endl;
		classwidthvecpathlength.push_back(std::to_string(pathlengthclasswidth*i) + " - " + std::to_string(pathlengthclasswidth * (i + 1)));

	}


	//***// //***// //***// 
	PathBase::setCSVfilename(OutputFilename);
	PathBase::writeCSV(argList, "list of cmd arguments", classwidthvecangles, "angle classes", histogram, "number of angles (overall)",
		perimhistogram, "number of angles (perim only)", hatchhistogram, "number of angles (hatching only)",
		histogram_vol_in_angle, "volume of pores per angle class (overall) (voxels)", histogram_vol_in_mm_angle, "volume of pores per angle class (overall) (mm^3)",
		perimhistogram_vol_in_angle, "volume of pores per angle class (perim) (voxels)", hatchinghistogram_vol_in_angle, "volume of pores per angle class (hatching) (voxels)",
		classwidthvecpathlength, "path length classes", pathlengthshisto, "pores per mm", pathlengthshisto_volumeweighted, "defect volume of pores per mm", frequenciesinclasses, "number of pores in path class", frequenciesinclasses_volumeweighted, "defect volume in path class",
		poresinlayer, "number of pores per layer", porevolumeinlayer, "overall volume of pores in layer", poreparser.poreorientationhisto_azim_overall, "Azimuth histogram of all pores from file -as is",
		poreparser.poreorientationhisto_azim_filtered, "Azimuth histogram of pores from file filtered with rules", poreparser.poreorientationhisto_elev_overall, "Elevation histogram of all pores from file -as is",
		poreparser.poreorientationhisto_elev_filtered, "Elevation histogram of pores from file filtered with rules", poreparser.poreazimuth_histo_overall_volume, "Azimuth histogram of all pores from file -as is - volume weighted",
		poreparser.poreazimuth_histo_filtered_volume, "Azimuth histogram of pores from file filtered with rules - volume weighted", poreparser.poreelevation_histo_overall_volume, "Elevation histogram of all pores from file -as is - volume weighted",
		poreparser.poreelevation_histo_filtered_volume, "Elevation histogram of pores from file filtered with rules - volume weighted",
		poreinformations, "information about pores");

	// output histos as csv file
	if (!vtkoutputoff) {
		writework.get(); // writing finished until here?
	}
	std::cout << "GCodePores finished with success!" << std::endl;
	return EXIT_SUCCESS;
}