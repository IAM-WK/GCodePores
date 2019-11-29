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
#include "PorePathCorrelate.h"

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
	constexpr float layer_height_tol = 0.001f; // 1 mue tolerance
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
	const double pixelvolume = static_cast<double>(ImageParser.pixelwidth[0]) * static_cast<double>(ImageParser.pixelwidth[1]) * static_cast<double>(ImageParser.pixelwidth[2]);
	std::cout << "volume of all pores = " << poreparser.porevolume << " voxels; => " << (poreparser.porevolume* pixelvolume) << "mm^3 \n";
	std::cout << "volume of considered pores = " << poreparser.filteredporesvolume << " voxels; => " << (poreparser.filteredporesvolume* pixelvolume) << "mm^3 \n";
	std::cout << "this corresponds to a fraction of " << (std::roundf(static_cast<float>(poreparser.filteredporesvolume) / static_cast<float>(poreparser.porevolume) * 1000)/10) << "% \n";
	std::cout << "fraction of pores (by number) that were classified as non gas pores: " << (std::roundf(static_cast<float>(poreparser.khpores) / static_cast<float>(poreparser.khpores+poreparser.gaspores) * 1000) / 10) << "% \n";
	// save this informations to a vec to output it to the csv file
	std::vector<std::string> poreinformations;
	poreinformations.push_back("volume of all pores = " + std::to_string(poreparser.porevolume) + " voxels; => " + std::to_string(poreparser.porevolume* pixelvolume) + "mm^3 ");
	poreinformations.push_back("volume of considered pores = " + std::to_string(poreparser.filteredporesvolume) + " voxels; => " + std::to_string(poreparser.filteredporesvolume* pixelvolume) + "mm^3 ");
	poreinformations.push_back("this corresponds to a fraction of " + std::to_string(std::roundf(static_cast<float>(poreparser.filteredporesvolume) / static_cast<float>(poreparser.porevolume) * 1000) / 10) + "% ");
	poreinformations.push_back("fraction of pores (by number) that were classified as non gas pores: " + std::to_string(std::roundf(static_cast<float>(poreparser.khpores) / static_cast<float>(poreparser.khpores + poreparser.gaspores) * 1000) / 10) + "% ");
	//*************** end output processed information **********************

	PoreParser::PoreVecType porevec;
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

	//***********************************************************************
	// PORES TO PATH CORRELATION
	//***********************************************************************
	PorePathCorrelate Correlator(&pathVec, &GCodeAnalyser.pathclassificationwithchunks, &GCodeAnalyser.layer_height, &GCodeAnalyser.heights_avail, &GCodeAnalyser.lengthofchunks, &porevec);
	Correlator.setInputSettings(layer_height_tol, layersconsidered, nextneighbours, classwidth, pathlengthclasswidth, pixelvolume);
	Correlator.setPoreCols(Xcol, Ycol, Zcol, azimuthcol, volumecol);
	Correlator.Execute();
	//***********************************************************************
	// END PORES TO PATH CORRELATION
	//***********************************************************************

	//***********************************************************************
	// OUTPUT CSV FILE
	//***********************************************************************
	PathBase::setCSVfilename(OutputFilename);
	PathBase::writeCSV(argList, "list of cmd arguments", Correlator.classwidthvecangles, "angle classes", Correlator.PoretoPathAngle_Histo, "number of angles (overall)",
		Correlator.PoretoPathAngle_Histo_Perim, "number of angles (perim only)", Correlator.PoretoPathAngle_Histo_Hatching, "number of angles (hatching only)",
		Correlator.PoretoPathAngle_Histo_VW, "volume of pores per angle class (overall) (voxels)", Correlator.PoretoPathAngle_Histo_VW_cmm, "volume of pores per angle class (overall) (mm^3)",
		Correlator.PoretoPathAngle_Histo_Perim_VW, "volume of pores per angle class (perim) (voxels)", Correlator.PoretoPathAngle_Histo_Hatching_VW, "volume of pores per angle class (hatching) (voxels)",
		Correlator.classwidthvecpathlength, "path length classes", Correlator.PoresperPathlength_Histo, "pores per mm", Correlator.PoresperPathlength_Histo_VW, "defect volume of pores per mm", Correlator.PoreinClassFrequencies, "number of pores in path class", Correlator.PoreinClassFrequencies_VW, "defect volume in path class",
		Correlator.PoreinLayerFrequencies, "number of pores per layer", Correlator.PoreinLayerFrequencies_VW, "overall volume of pores in layer", poreparser.poreorientationhisto_azim_overall, "Azimuth histogram of all pores from file -as is",
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