#include "Grid.h"
#include <math.h>
#include <iomanip>
#include <stdexcept>
#include <stdio.h>

Grid::Grid (double height, double width, double resolution,
		double upperLeftLatitude, double U, double initHeight,
		double simulationSeconds)
{
	this->height = height;
	this->width = width;
	this->resolution = resolution;
	this->upperLeftLatitude = upperLeftLatitude;
	this->simulationSeconds = simulationSeconds;
	this->U = U;
	this->initHeight = initHeight * 1000;

	// The timestep is determined by the Courant condition number.
	this->dt = 0.9 * resolution * 1000.0 / U;
	this->timestepNum = 0;

	// Compute, how many computation cells are needed and add 2 because of the boundary cells.
	int numRows = ceil(height / resolution) + 2;
	int numCols = ceil(width / resolution) + 2;

	// Reserve space for grid matrix.
	gridCells.resize(numRows);
	for (int i = 0; i < numRows; i++)
		gridCells.at(i).reserve(numCols);

	// Create the output file
	char const* filename = "shallowWaterOut.nc";
	ncfile = new NcFile(filename, NcFile::Replace);

	const NcDim* t1 = ncfile->add_dim("time");
	t = ncfile->add_var("time", ncDouble, 1, &t1);
	dx = ncfile->add_var("dx", ncDouble, 1, &t1);
	dx->put(&resolution, 1);

	const NcDim* dims[2] = { ncfile->add_dim("x", numCols), ncfile->add_dim(
			"y", numRows) };
	u = ncfile->add_var("u", ncDouble, t1, dims[0], dims[1]);

	v = ncfile->add_var("v", ncDouble, t1, dims[0], dims[1]);

	h = ncfile->add_var("h", ncDouble, t1, dims[0], dims[1]);

	htot = ncfile->add_var("htot", ncDouble, t1, dims[0], dims[1]);

	// Initialize grid values
	setInitialConditions();
}

void Grid::setInitialConditions ()
{
	int numRows = (int) gridCells.size();
	int numCols = (int) gridCells[1].capacity();

	for (int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numCols; j++)
		{
			GridCell g(j, 0.0, initHeight, 0.0, 0.0);
			gridCells.at(i).push_back(g);
		}
	}

	writeVariablesToFile();
}

void Grid::writeVariablesToFile ()
{
	double simTime = dt * timestepNum;

	t->put(&simTime, 1);

	int numRows = (int) gridCells.size();
	int numCols = (int) gridCells[1].size();

	std::vector<double > uVec, vVec, hVec, htotVec;
	for (int j = 0; j < numCols; j++)
	{
		for (int i = 0; i < numRows; i++)
		{
			uVec.push_back(gridCells.at(i).at(j).getU());
			vVec.push_back(gridCells.at(i).at(j).getV());
			hVec.push_back(gridCells.at(i).at(j).getH());
			htotVec.push_back(gridCells.at(i).at(j).getHtot());
		}
	}

	printf("%12.3g\n", uVec[0]);

	u->put(&uVec[0], 1, numCols, numRows);
	v->put(&vVec[0], 1, numCols, numRows);
	h->put(&hVec[0], 1, numCols, numRows);
	htot->put(&htotVec[0], 1, numCols, numRows);
}

Grid::~Grid()
{
	delete ncfile; // Flush the i/o and close the file.
}

