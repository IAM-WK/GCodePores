#include "mhdParser.h"



void mhdParser::setFilename(const std::string& Filename)
{
	mhdfilename = Filename;
}


void mhdParser::Execute()
{

	std::ifstream ImageFile;
	ImageFile.open(this->mhdfilename);
	std::string line_string;
	const std::string dimsize = "DimSize", pixsize = "ElementSize";
	size_t pixsizeposx = std::string::npos, pixsizeposy = std::string::npos, pixsizeposz = std::string::npos, sizeposx = std::string::npos, sizeposy = std::string::npos, sizeposz = std::string::npos;
	size_t endpixposx = std::string::npos, endpixposy = std::string::npos, endpixposz = std::string::npos, endsizeposx = std::string::npos, endsizeposy = std::string::npos, endsizeposz = std::string::npos;
	const std::string coordchars = "0123456789.+-"; //allowed characters in a coordinate specification // +- for robustness

	
	imagedimensionvoxels.fill(0); imagedimensionmm.fill(0);

	while (getline(ImageFile, line_string))
	{
		// dimension
		if (line_string.substr(0, dimsize.length()) == dimsize) {
			sizeposx = line_string.find("= "); // x is first after "="
			endsizeposx = line_string.substr(sizeposx + 2, line_string.length()).find_first_not_of(coordchars); //length of coord
			imagedimensionvoxels[0] = stof(line_string.substr(sizeposx + 1, endsizeposx + 1));

			sizeposy = line_string.find(" ", sizeposx + endsizeposx); // y after next space
			endsizeposy = line_string.substr(sizeposy + 2, line_string.length()).find_first_not_of(coordchars); //length of coord
			imagedimensionvoxels[1] = stof(line_string.substr(sizeposy + 1, endsizeposy + 1));

			sizeposz = line_string.find(" ", sizeposy + endsizeposy); // z after next space
			endsizeposz = line_string.substr(sizeposz + 2, line_string.length()).find_first_not_of(coordchars); //length of coord
			imagedimensionvoxels[2] = stof(line_string.substr(sizeposz + 1, endsizeposz));

		}
		// pixel size
		if (line_string.substr(0, pixsize.length()) == pixsize) {
			pixsizeposx = line_string.find("= "); // x is first after "="
			endpixposx = line_string.substr(pixsizeposx + 2, line_string.length()).find_first_not_of(coordchars); //length of coord
			this->pixelwidth[0] = stof(line_string.substr(pixsizeposx + 1, endpixposx + 1));

			pixsizeposy = line_string.find(" ", pixsizeposx + endpixposx); // y after next space
			endpixposy = line_string.substr(pixsizeposy + 2, line_string.length()).find_first_not_of(coordchars); //length of coord
			this->pixelwidth[1] = stof(line_string.substr(pixsizeposy + 1, endpixposy + 1));

			pixsizeposz = line_string.find(" ", pixsizeposy + endpixposy); // z after next space
			endpixposz = line_string.substr(pixsizeposz + 2, line_string.length()).find_first_not_of(coordchars); //length of coord
			this->pixelwidth[2] = stof(line_string.substr(pixsizeposz + 1, endpixposz));
		}

	}

	this->imagedimensionmm[0] = (imagedimensionvoxels[0]) * pixelwidth[0];
	this->imagedimensionmm[1] = (imagedimensionvoxels[1]) * pixelwidth[1];
	this->imagedimensionmm[2] = (imagedimensionvoxels[2]) * pixelwidth[2];

}