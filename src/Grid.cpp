#include "Grid.h"
#include <math.h>
#include <iomanip>
#include <stdexcept>
#include <stdio.h>

// PUBLIC MEMBER FUNCTIONS &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// CONSTRUCTOR ****************************************************
Grid::Grid (double height, double width, double resolution,
		double upperLeftLatitude, double U, double depth,
		double simulationSeconds, double writeInterval, double mountHeight)
{
	// Set the member variables to given values.
	this->height = height;
	this->width = width;
	this->resolution = resolution;
	this->upperLeftLatitude = upperLeftLatitude;
	this->simulationLength = simulationSeconds;
	this->U = U;
	this->depth = depth;
	this->writeInterval = writeInterval;
	this->mountainHeight = mountHeight;

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

	// Create an allocatable dimension for time.
	const NcDim* t1 = ncfile->add_dim("time");
	t = ncfile->add_var("time", ncDouble, 1, &t1);
	// Store resolution.
	dx = ncfile->add_var("dx", ncDouble, 1, &t1);
	dx->put(&resolution, 1);
	// Store the latitude.
	lat = ncfile->add_var("lat", ncDouble, 1, &t1);
	lat->put(&upperLeftLatitude, 1);

	// Create the remaining dimensions and variables.
	const NcDim* dims[2] = { ncfile->add_dim("x", numCols), ncfile->add_dim("y",
			numRows) };
	u = ncfile->add_var("u", ncDouble, t1, dims[0], dims[1]);
	v = ncfile->add_var("v", ncDouble, t1, dims[0], dims[1]);
	h = ncfile->add_var("h", ncDouble, t1, dims[0], dims[1]);
	htot = ncfile->add_var("htot", ncDouble, t1, dims[0], dims[1]);

	// Initialize grid values
	setInitialConditions();
}
// ****************************************************************

// MAIN SIMULATION FUNCTION ***************************************
void Grid::simulate ()
{
	// Start counting from 1, because initial condition is defined as iteration 0.
	long timestepNum = 1, writeCount = 1;

	// Set timers to zero.
	double elapsedSeconds = 0.0; // Elapsed seconds since the simulation beginning.
	double previousWriteTime = 0.0; // Elapsed seconds since last file write time.

	while (elapsedSeconds <= simulationLength)
	{
		elapsedSeconds += dt; // Timestep dt is NOT constant.
		advanceGrid(timestepNum); // Compute the MacCormack scheme.

		// Write simulation fields to file, if enough seconds has passed.
		if (elapsedSeconds - previousWriteTime >= writeInterval)
		{
			writeVariablesToFile(writeCount, elapsedSeconds);
			previousWriteTime = elapsedSeconds;
			writeCount++;
		}
		timestepNum++;
	}

	// Write final results to file.
	if (previousWriteTime != elapsedSeconds)
	{
		writeVariablesToFile(writeCount, elapsedSeconds);
	}
}
// ****************************************************************

// DESTRUCTOR *****************************************************
Grid::~Grid ()
{
	delete ncfile; // Flush the output and close the file.
}
// ****************************************************************
// END OF PUBLIC FUNCTIONS &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

// PRIVATE FUNCTIONS &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
void Grid::advanceGrid (long stepNum)
{
	// Get grid dimensions.
	int numRows = (int) gridCells.size();
	int numCols = (int) gridCells[1].capacity();

	#pragma omp parallel
	{
		// First compute all the predictors.
		#pragma omp for schedule(static)
		for (int i = 0; i < numRows; i++)
		{
			for (int j = 0; j < numCols; j++)
			{
				gridCells[i][j].computePredictor(stepNum);
			}
		}

		// Because all the predictors must be known before computing any of the correctors,
		// traverse the grid again, now computing the correctors and updating the values.
		#pragma omp for schedule(static)
		for (int i = 0; i < numRows; i++)
		{
			for (int j = 0; j < numCols; j++)
			{
				gridCells[i][j].computeCorrector(stepNum);
			}
		}
	}

	// Find the maximum velocity for the computation of the new timestep.
	double maxVelocity = 0.0;
	for (int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numCols; j++)
		{
			double u = gridCells[i][j].getU();
			double v = gridCells[i][j].getV();
			// The Courant number requires component addition in linear fashion.
			double vel = u + v;
			if (vel > maxVelocity)
			{
				maxVelocity = vel;
			}
		}
	}

	// Limit the new timestep between 0.1 and 60 seconds.
	dt = fmax(0.1,
			fmin(courantNumber * resolution * 1000.0 / maxVelocity, 60.0));

	// Set the new timestep for all of the grid cells.
	// (This could also be implemented as a static member variable).
	for (int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numCols; j++)
		{
			gridCells[i][j].setTimestep(dt);
		}
	}
}

void Grid::setInitialConditions ()
{
	// Initialize the grid.
	initializeCenterCells();
	initializeBoundaryCells();

	// Write initial conditions to file.
	writeVariablesToFile(0, 0.0);
}

void Grid::initializeCenterCells ()
{
	// Get dimensions.
	int numRows = (int) gridCells.size();
	int numCols = (int) gridCells[1].capacity();

	int mountPos = numCols / 4; // In longitudinal sense.
	int mountMin = numRows / 4; // Northern boundary for the mountain range.
	int mountMax = 3 * numRows / 4; // Southern boundary.

	// Some variables to prevent sharp mountain edges near the sides.
	int slopeBegin = numRows / 8;
	int slopeEnd = 7 * numRows / 8;
	double mountHeight = mountainHeight;
	double slope = mountHeight / (mountMin - slopeBegin);
	double mountWidth = 33.5 / numCols; // = 0.05 in the default configuration.

	// Correction to the initial height, which is necessary for the
	// pressure gradient.
	double htotCorrection = 0.0;

	// Initialize grid cells with given initial conditions.
	for (int i = 0; i < numRows; i++)
	{
		double lat = computeCellLatitude(i);
		htotCorrection = computeHtotCorrection(lat);
		for (int j = 0; j < numCols; j++)
		{
			// Call grid cell constructor.
			GridCell g(0.01, 0.0, depth + htotCorrection, lat, resolution, dt);
			// Push the grid cell into the grid.
			gridCells.at(i).push_back(g);
		}
	}

	// Set neighbors for center cells.
	for (int i = 1; i < numRows - 1; i++)
	{
		for (int j = 1; j < numCols - 1; j++)
		{
			GridCell* g = &gridCells[i][j];
			g->setNorthernNeighbour(&gridCells.at(i - 1).at(j));
			g->setSouthernNeighbour(&gridCells.at(i + 1).at(j));
			g->setEasternNeighbour(&gridCells.at(i).at(j + 1));
			g->setWesternNeighbour(&gridCells.at(i).at(j - 1));
			g->setBoundaryType('0'); // 0 means not on boundary.
		}
	}

	// Create the mountain.
	for (int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numCols; j++)
		{
			double surfHeight = 0;
			if (i >= mountMin && i <= mountMax)
			{
				surfHeight = mountHeight
						* exp(-mountWidth * pow(j - mountPos, 2.0));
			}
			else if (i >= slopeBegin && i < mountMin)
			{
				surfHeight = slope * (i - slopeBegin)
						* exp(-mountWidth * pow(j - mountPos, 2.0));
			}
			else if (i > mountMax && i <= slopeEnd)
			{
				surfHeight = slope * (slopeEnd - i)
						* exp(-mountWidth * pow(j - mountPos, 2.0));
			}
			gridCells[i][j].setSurfaceHeight(surfHeight);
		}
	}
}

void Grid::initializeBoundaryCells ()
{
	// Get dimensions.
	int numRows = (int) gridCells.size();
	int numCols = (int) gridCells[1].capacity();

	// Set boundary conditions.
	GridCell::northH = depth;
	GridCell::northV = 0.0;

	GridCell::southV = 0.0;

	GridCell::westU = U;
	GridCell::westV = 0.0;

	GridCell* g; // Temporary gridCell.

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

double Grid::computeCellLatitude (int i)
{
	return upperLeftLatitude - 360.0 * resolution * i / (2.0 * M_PI * 6371.0);
}

double Grid::computeHtotCorrection (double lat)
{
	// 2 * omega * U * a * (cos(lat) - cos(lat_0)) / g;
	return 2.0 * 7.2921150E-5 * U * 6371E3
			* (cos(lat * M_PI / 180.0) - cos(upperLeftLatitude * M_PI / 180.0))
			/ 9.80665;
}

void Grid::writeVariablesToFile (long timestepNum, double elapsedSeconds)
{
	t->set_cur(timestepNum); // Set current time index.
	t->put(&elapsedSeconds, 1); // Put elapsed time.

	int numRows = (int) gridCells.size();
	int numCols = (int) gridCells[1].size();

	// Collect the two-dimensional output fields.
	std::vector<double> uVec, vVec, hVec, htotVec;
	for (int j = 0; j < numCols; j++)
	{
		for (int i = 0; i < numRows; i++)
		{
			uVec.push_back(gridCells[i][j].getU());
			vVec.push_back(gridCells[i][j].getV());
			hVec.push_back(gridCells[i][j].getH());
			htotVec.push_back(gridCells[i][j].getHtot());
		}
	}

	// Print progress.
	if (timestepNum == 0)
	{
		printf("Progress:\n");
	}
	printf("%4.1f %%\n",
			timestepNum / (simulationLength / writeInterval) * 100);

	// Write the variables.
	u->set_cur(timestepNum, 0, 0);
	u->put(&uVec[0], 1, numCols, numRows);

	v->set_cur(timestepNum, 0, 0);
	v->put(&vVec[0], 1, numCols, numRows);

	h->set_cur(timestepNum, 0, 0);
	h->put(&hVec[0], 1, numCols, numRows);

	htot->set_cur(timestepNum, 0, 0);
	htot->put(&htotVec[0], 1, numCols, numRows);
}
// END OF PRIVATE FUNCTIONS &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

