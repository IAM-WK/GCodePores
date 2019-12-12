#include "GCodeTransform.h"


GCodeTransform::GCodeTransform()
{
}

void GCodeTransform::setOrigin(const std::vector<float>& ImageOriginVal) noexcept
{
	this->ImageOrigin = ImageOriginVal;
}

void GCodeTransform::setAngles(const float & ImageAngleGammaVal, const float & ImageAngleBetaVal, const float & ImageAngleAlphaVal) noexcept
{
	this->ImageAngleGamma = ImageAngleGammaVal;
	this->ImageAngleBeta = ImageAngleBetaVal;
	this->ImageAngleAlpha = ImageAngleAlphaVal;
}

void GCodeTransform::setImagedimensions(const float & y_image_lengthVal, const float & z_image_lengthVal) noexcept
{
	this->y_image_length = y_image_lengthVal;
	this->z_image_length = z_image_lengthVal;
}

void GCodeTransform::setInput(const PathBase::PathVector &GCodeVector) noexcept
{
	this->ImageCoordinateGCode = GCodeVector; // copy values to new vector
}



void GCodeTransform::Execute()
{
	
	// make vector relative: start with 0,0,0 and subtract first entry from every entry
	std::array<float, 3> printpathbegin = { this->ImageCoordinateGCode.front().front()[0], this->ImageCoordinateGCode.front().front()[1] ,this->ImageCoordinateGCode.front().front()[2] };
	for (auto &chunkitem : this->ImageCoordinateGCode) {
		for (auto &coorditem : chunkitem) {
			for (std::size_t k = 0; k < 3; ++k) {
				coorditem[k] = coorditem[k] - printpathbegin[k] ; // make vector relative to first printing point
			}
		}
	}

	// first do an transformation of KOS : Roation and translation (Alias transform)
	// GCodeKOS to ImageKOS --> 180° rotation about x-axis ----------> always like this?
	// --> invert y-axis and z-axis
	for (auto &chunkitem : this->ImageCoordinateGCode) {
		for (auto &coorditem : chunkitem) {
			coorditem[1] = (-1)*coorditem[1]; // invert y
			coorditem[2] = (-1)*coorditem[2]; // invert z
			}
	}

	// GCodeKOS to ImageKOS --> translate z about h_image; y about length_image
	for (auto &chunkitem : this->ImageCoordinateGCode) {
		for (auto &coorditem : chunkitem) {
			coorditem[1] = coorditem[1] + this->y_image_length; // translate y
			coorditem[2] = coorditem[2] + this->z_image_length; // translate z
		}
	}

	// then rotate+translate object accordingly
	std::array<float, 3> DifferenceVector;
	
	// translate GCode to imageorigin
	// GCodeOrigin = PrintVector[0][0]
	// DifferenceVector = ImageOrigin - GCodeOrigin
	// loop: PrintVector[i][j] + DifferenceVector
	DifferenceVector[0] = this->ImageOrigin[0] - this->ImageCoordinateGCode.front().front()[0];
	DifferenceVector[1] = this->ImageOrigin[1] - this->ImageCoordinateGCode.front().front()[1];
	DifferenceVector[2] = this->ImageOrigin[2] - this->ImageCoordinateGCode.front().front()[2];

	for (auto &pathitem : this->ImageCoordinateGCode) {
		for (auto &coorditem : pathitem) {
			for (std::size_t k = 0; k < 3; ++k) {
				coorditem[k] =  coorditem[k] + DifferenceVector[k]; //   use ImageCoordinateGCode
			}
		}
	}

	// rotate GCode with angle
	// rotate XYZ
	// rotate dxdydz or recalculate
	// x' = xursp + cos * (x - x1) - sin * (y - y1)
	// y' = yursp + sin * (x - x1) + cos * (y - y1)   
	// rotation about z-axis
	float tempx, tempy, tempz;
	for (auto &pathitem : this->ImageCoordinateGCode) {
		for (auto &coorditem : pathitem) {
			tempx = this->ImageOrigin[0] + (coorditem[0] - this->ImageOrigin[0]) * cos(this->ImageAngleGamma) - (coorditem[1] - this->ImageOrigin[1]) * sin(this->ImageAngleGamma);
			tempy = this->ImageOrigin[1] + (coorditem[0] - this->ImageOrigin[0]) * sin(this->ImageAngleGamma) + (coorditem[1] - this->ImageOrigin[1]) * cos(this->ImageAngleGamma);
			coorditem[0] = tempx;
			coorditem[1] = tempy;
		}
	}
	// rotation about y-axis
	for (auto &pathitem : this->ImageCoordinateGCode) {
		for (auto &coorditem : pathitem) {
			tempx = this->ImageOrigin[0] + (coorditem[0] - this->ImageOrigin[0]) * cos(this->ImageAngleBeta) + (coorditem[2] - this->ImageOrigin[2]) * sin(this->ImageAngleBeta);
			tempz = this->ImageOrigin[2] - (coorditem[0] - this->ImageOrigin[0]) * sin(this->ImageAngleBeta) + (coorditem[2] - this->ImageOrigin[2]) * cos(this->ImageAngleBeta);
			coorditem[0] = tempx;
			coorditem[2] = tempz;
		}
	}
	// rotation about x-axis
	for (auto &pathitem : this->ImageCoordinateGCode) {
		for (auto &coorditem : pathitem) {
			tempy = this->ImageOrigin[1] + (coorditem[1] - this->ImageOrigin[1]) * cos(this->ImageAngleAlpha) - (coorditem[2] - this->ImageOrigin[2]) * sin(this->ImageAngleAlpha);
			tempz = this->ImageOrigin[2] + (coorditem[1] - this->ImageOrigin[1]) * sin(this->ImageAngleAlpha) + (coorditem[2] - this->ImageOrigin[2]) * cos(this->ImageAngleAlpha);
			coorditem[1] = tempy;
			coorditem[2] = tempz;
		}
	}


	// now we need to regenerate the directions
	float dx, dy, dz, disttonext;
	for (auto &chunkitem : this->ImageCoordinateGCode) {
		for (std::size_t i = 0; i < chunkitem.size() - 1; ++i) {
			dx = chunkitem[i + 1][0] - chunkitem[i][0];
			dy = chunkitem[i + 1][1] - chunkitem[i][1];
			dz = chunkitem[i + 1][2] - chunkitem[i][2];
			disttonext = sqrt(dx * dx + dy * dy + dz * dz);
			dx = dx / disttonext; dy = dy / disttonext; dz = dz / disttonext; //norm direction
			chunkitem[i][3] = dx; chunkitem[i][4] = dy; chunkitem[i][5] = dz; chunkitem[i][6] = disttonext;
			// prints out extracted object printing coords
			//std::cout << " -- X: " << pathitem[i][0] << " -- Y: " << pathitem[i][1] << " -- Z: " << pathitem[i][2] << " -- disttonext " << pathitem[i][6] << std::endl;
		}
	}
	for (std::size_t i = 0; i < this->ImageCoordinateGCode.size(); ++i) { // update last elements
		if (this->ImageCoordinateGCode[i].size() >= 2) {
			dx = this->ImageCoordinateGCode[i][this->ImageCoordinateGCode[i].size() - 2][3];
			dy = this->ImageCoordinateGCode[i][this->ImageCoordinateGCode[i].size() - 2][4]; //caution with double travels !!!!!!
			dz = this->ImageCoordinateGCode[i][this->ImageCoordinateGCode[i].size() - 2][5];
		}
		else { // take last direction of -1 -1 -1 as fallback ----> Todo!
			dx = -1;
			dy = -1;
			dz = -1;
		}
		// disttonext is 0 by definition
		this->ImageCoordinateGCode[i].back()[3] = dx; this->ImageCoordinateGCode[i].back()[4] = dy; this->ImageCoordinateGCode[i].back()[5] = dz; this->ImageCoordinateGCode[i].back()[6] = 0;
	}



	return;
}

PathBase::PathVector GCodeTransform::GetOutput() const
{
	if (this->ImageCoordinateGCode.size() == 0) {
		throw std::runtime_error("Error after path transformation: Transformed coordinate vector is empty!");
	}
	return this->ImageCoordinateGCode;
}

