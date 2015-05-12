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

	const NcDim* dims[2] = { ncfile->add_dim("x", numCols), ncfile->add_dim("y",
			numRows) };
	u = ncfile->add_var("u", ncDouble, t1, dims[0], dims[1]);

	v = ncfile->add_var("v", ncDouble, t1, dims[0], dims[1]);

	h = ncfile->add_var("h", ncDouble, t1, dims[0], dims[1]);

	htot = ncfile->add_var("htot", ncDouble, t1, dims[0], dims[1]);

	// Initialize grid values
	setInitialConditions();
}

void Grid::setInitialConditions ()
{
	initializeCenterCells();
	initializeBoundaryCells();
	writeVariablesToFile();
}

void Grid::initializeCenterCells ()
{
	// Get dimensions.
	int numRows = (int) gridCells.size();
	int numCols = (int) gridCells[1].capacity();

	// Initialize ALL grid cells with given initial conditions.
	for (int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numCols; j++)
		{
			GridCell g(j, 0.0, initHeight, 0.0, 0.0, resolution, dt);
			gridCells.at(i).push_back(g);
		}
	}

	// Set neighbors for center cells.
	for (int i = 1; i < numRows - 2; i++)
	{
		for (int j = 1; j < numCols - 2; j++)
		{
			GridCell* g = &gridCells.at(i).at(j);
			g->setNorthernNeighbour(&gridCells.at(i - 1).at(j));
			g->setSouthernNeighbour(&gridCells.at(i + 1).at(j));
			g->setEasternNeighbour(&gridCells.at(i).at(j + 1));
			g->setWesternNeighbour(&gridCells.at(i).at(j - 1));
			g->setBoundaryType('0');
		}
	}
}

void Grid::initializeBoundaryCells ()
{
	// Get dimensions.
	int numRows = (int) gridCells.size();
	int numCols = (int) gridCells[1].capacity();

	// Set Western and Eastern boundary cells' neighbors, excluding corners.
	for (int i = 1; i < numRows - 2; i++)
	{
		GridCell* g = &gridCells.at(i).at(0);
		g->setNorthernNeighbour(&gridCells.at(i - 1).at(0));
		g->setSouthernNeighbour(&gridCells.at(i + 1).at(0));
		g->setEasternNeighbour(&gridCells.at(i).at(1));
		g->setWesternNeighbour(0);
		g->setBoundaryType('W');

		GridCell* g = &gridCells.at(i).at(numCols - 1);
		g->setNorthernNeighbour(&gridCells.at(i - 1).at(numCols - 1));
		g->setSouthernNeighbour(&gridCells.at(i + 1).at(numCols - 1));
		g->setEasternNeighbour(0);
		g->setWesternNeighbour(&gridCells.at(i).at(numCols - 2));
		g->setBoundaryType('E');
	}

	// Set Northern and Southern boundary cells' neighbors, excluding corners.
	for (int j = 1; j < numCols - 2; j++)
	{
		GridCell* g = &gridCells.at(0).at(j);
		g->setNorthernNeighbour(0);
		g->setSouthernNeighbour(&gridCells.at(1).at(j));
		g->setEasternNeighbour(&gridCells.at(0).at(j + 1));
		g->setWesternNeighbour(&gridCells.at(0).at(j - 1));
		g->setBoundaryType('N');

		GridCell* g = &gridCells.at(numRows - 1).at(j);
		g->setNorthernNeighbour(&gridCells.at(numRows - 2).at(j));
		g->setSouthernNeighbour(0);
		g->setEasternNeighbour(&gridCells.at(numRows - 1).at(j + 1));
		g->setWesternNeighbour(&gridCells.at(numRows - 1).at(j - 1));
		g->setBoundaryType('S');
	}

	// Set neighbors for corner cells. 1 = NW, 2 = NE, 3 = SW, 4 = SE corner.

	GridCell* g = &gridCells.at(0).at(0);
	g->setNorthernNeighbour(0);
	g->setSouthernNeighbour(&gridCells.at(1).at(0));
	g->setEasternNeighbour(&gridCells.at(0).at(1));
	g->setWesternNeighbour(0);
	g->setBoundaryType('1');

	GridCell* g = &gridCells.at(0).at(numCols-1);
	g->setNorthernNeighbour(0);
	g->setSouthernNeighbour(&gridCells.at(1).at(numCols-1));
	g->setEasternNeighbour(0);
	g->setWesternNeighbour(&gridCells.at(0).at(numCols-2));
	g->setBoundaryType('2');

	GridCell* g = &gridCells.at(numRows-1).at(0);
	g->setNorthernNeighbour(&gridCells.at(numRows-2).at(0));
	g->setSouthernNeighbour(0);
	g->setEasternNeighbour(&gridCells.at(numRows-1).at(1));
	g->setWesternNeighbour(0);
	g->setBoundaryType('3');

	GridCell* g = &gridCells.at(numRows-1).at(numCols-1);
	g->setNorthernNeighbour(&gridCells.at(numRows-2).at(numCols-1));
	g->setSouthernNeighbour(0);
	g->setEasternNeighbour(0);
	g->setWesternNeighbour(&gridCells.at(numRows-1).at(numCols-2));
	g->setBoundaryType('4');
}

void Grid::writeVariablesToFile ()
{
	double simTime = dt * timestepNum;

	t->put(&simTime, 1);

	int numRows = (int) gridCells.size();
	int numCols = (int) gridCells[1].size();

	std::vector<double> uVec, vVec, hVec, htotVec;
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

Grid::~Grid ()
{
	delete ncfile; // Flush the i/o and close the file.
}

