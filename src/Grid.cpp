#include "Grid.h"
#include <math.h>
#include <iomanip>
#include <stdexcept>
#include <stdio.h>

Grid::Grid (double height, double width, double resolution,
		double upperLeftLatitude, double U, double initHeight,
		double simulationSeconds, double writeInterval)
{
	this->height = height;
	this->width = width;
	this->resolution = resolution;
	this->upperLeftLatitude = upperLeftLatitude;
	this->simulationLength = simulationSeconds;
	this->U = U;
	this->initHeight = initHeight * 1000;
	this->writeInterval = writeInterval;

	// The timestep is determined by the Courant condition number.
	this->courantNumber = 0.01;
	this->dt = 0.01; // Start out with a VERY small dt. This will be adapted later.

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

void Grid::simulate ()
{
	long timestepNum = 1, writeCount = 1;
	double elapsedSeconds = 0.0;
	double previousWriteTime = 0.0;

	while (elapsedSeconds <= simulationLength)
	{
		elapsedSeconds += dt;
		advanceGrid(timestepNum);
		if (elapsedSeconds - previousWriteTime >= writeInterval)
		{
			writeVariablesToFile(writeCount, elapsedSeconds);
			previousWriteTime = elapsedSeconds;
			writeCount++;
		}
		timestepNum++;
	}
}

void Grid::advanceGrid (long stepNum)
{
	// Get grid dimensions.
	int numRows = (int) gridCells.size();
	int numCols = (int) gridCells[1].capacity();

	// First compute predictors for all cells.
#pragma omp parallel for
	for (int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numCols; j++)
		{
			gridCells.at(i).at(j).computePredictor(stepNum);
		}
	}

	// Because all predictors must be known before computing any corrector,
	// traverse the grid again, now computing the correctors and updating the values.
#pragma omp for
	for (int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numCols; j++)
		{
			gridCells.at(i).at(j).computeCorrector(stepNum);
		}
	}


	double maxVelocity = 0.0;
	for (int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numCols; j++)
		{
			double u = gridCells.at(i).at(j).getU();
			double v = gridCells.at(i).at(j).getV();
			// The Courant number requires component addition in linear fashion.
			double vel = u + v;
			if (vel > maxVelocity)
			{
				maxVelocity = vel;
			}
		}
	}

//	if (stepNum > 1)
//	{
//		for (int i = 0; i < numRows; i++)
//		{
//			for (int j = 0; j < numCols; j++)
//			{
//				gridCells.at(i).at(j).clampSolution(stepNum);
//			}
//		}
//	}

	dt = courantNumber * resolution * 1000.0 / maxVelocity;

	for (int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numCols; j++)
		{
			gridCells.at(i).at(j).setTimestep(dt);
		}
	}
}

void Grid::setInitialConditions ()
{
	initializeCenterCells();
	initializeBoundaryCells();
	writeVariablesToFile(0, 0.0);
}

void Grid::initializeCenterCells ()
{
	// Get dimensions.
	int numRows = (int) gridCells.size();
	int numCols = (int) gridCells[1].capacity();
	int mountPos = numCols / 4;
	int ic = numRows / 2;
	int jc = numCols / 2;

	// Initialize ALL grid cells with given initial conditions.
	for (int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numCols; j++)
		{
			double drop = 15.0
					* exp(0.01 * (pow(i - ic, 2.0) + pow(j - jc, 2.0)));
			GridCell g(U, 0.0, initHeight, 00.0, resolution, dt);
			gridCells.at(i).push_back(g);
		}
	}

	// Set neighbors for center cells.
	for (int i = 1; i < numRows - 1; i++)
	{
		for (int j = 1; j < numCols - 1; j++)
		{
			GridCell* g = &gridCells.at(i).at(j);
			g->setNorthernNeighbour(&gridCells.at(i - 1).at(j));
			g->setSouthernNeighbour(&gridCells.at(i + 1).at(j));
			g->setEasternNeighbour(&gridCells.at(i).at(j + 1));
			g->setWesternNeighbour(&gridCells.at(i).at(j - 1));
			g->setBoundaryType('0');
			double surfHeight = 2000.0 * exp(-0.25 * pow(j - mountPos, 2.0));
			g->setSurfaceHeight(surfHeight);
		}
	}
}

void Grid::initializeBoundaryCells ()
{
	// Get dimensions.
	int numRows = (int) gridCells.size();
	int numCols = (int) gridCells[1].capacity();

	GridCell* g;

	// Set Western and Eastern boundary cells' neighbors, excluding corners.
	for (int i = 1; i < numRows - 1; i++)
	{
		g = &gridCells.at(i).at(0);
		g->setNorthernNeighbour(&gridCells.at(i - 1).at(0));
		g->setSouthernNeighbour(&gridCells.at(i + 1).at(0));
		g->setEasternNeighbour(&gridCells.at(i).at(1));
		g->setWesternNeighbour(0);
		g->setBoundaryType('W');

		g = &gridCells.at(i).at(numCols - 1);
		g->setNorthernNeighbour(&gridCells.at(i - 1).at(numCols - 1));
		g->setSouthernNeighbour(&gridCells.at(i + 1).at(numCols - 1));
		g->setEasternNeighbour(0);
		g->setWesternNeighbour(&gridCells.at(i).at(numCols - 2));
		g->setBoundaryType('E');
	}

	// Set Northern and Southern boundary cells' neighbors, excluding corners.
	for (int j = 1; j < numCols - 1; j++)
	{
		g = &gridCells.at(0).at(j);
		g->setNorthernNeighbour(0);
		g->setSouthernNeighbour(&gridCells.at(1).at(j));
		g->setEasternNeighbour(&gridCells.at(0).at(j + 1));
		g->setWesternNeighbour(&gridCells.at(0).at(j - 1));
		g->setBoundaryType('N');

		g = &gridCells.at(numRows - 1).at(j);
		g->setNorthernNeighbour(&gridCells.at(numRows - 2).at(j));
		g->setSouthernNeighbour(0);
		g->setEasternNeighbour(&gridCells.at(numRows - 1).at(j + 1));
		g->setWesternNeighbour(&gridCells.at(numRows - 1).at(j - 1));
		g->setBoundaryType('S');
	}

	// Set neighbors for corner cells. 1 = NW, 2 = NE, 3 = SW, 4 = SE corner.

	g = &gridCells.at(0).at(0);
	g->setNorthernNeighbour(0);
	g->setSouthernNeighbour(&gridCells.at(1).at(0));
	g->setEasternNeighbour(&gridCells.at(0).at(1));
	g->setWesternNeighbour(0);
	g->setBoundaryType('1');

	g = &gridCells.at(0).at(numCols - 1);
	g->setNorthernNeighbour(0);
	g->setSouthernNeighbour(&gridCells.at(1).at(numCols - 1));
	g->setEasternNeighbour(0);
	g->setWesternNeighbour(&gridCells.at(0).at(numCols - 2));
	g->setBoundaryType('2');

	g = &gridCells.at(numRows - 1).at(0);
	g->setNorthernNeighbour(&gridCells.at(numRows - 2).at(0));
	g->setSouthernNeighbour(0);
	g->setEasternNeighbour(&gridCells.at(numRows - 1).at(1));
	g->setWesternNeighbour(0);
	g->setBoundaryType('3');

	g = &gridCells.at(numRows - 1).at(numCols - 1);
	g->setNorthernNeighbour(&gridCells.at(numRows - 2).at(numCols - 1));
	g->setSouthernNeighbour(0);
	g->setEasternNeighbour(0);
	g->setWesternNeighbour(&gridCells.at(numRows - 1).at(numCols - 2));
	g->setBoundaryType('4');
}

void Grid::writeVariablesToFile (long timestepNum, double elapsedSeconds)
{

	t->set_cur(timestepNum);
	t->put(&elapsedSeconds, 1);

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

	printf("%ld\n", timestepNum);

	u->set_cur(timestepNum, 0, 0);
	u->put(&uVec[0], 1, numCols, numRows);

	v->set_cur(timestepNum, 0, 0);
	v->put(&vVec[0], 1, numCols, numRows);

	h->set_cur(timestepNum, 0, 0);
	h->put(&hVec[0], 1, numCols, numRows);

	htot->set_cur(timestepNum, 0, 0);
	htot->put(&htotVec[0], 1, numCols, numRows);

}

Grid::~Grid ()
{
	delete ncfile; // Flush the output and close the file.
}

