#ifndef GRID_H_
#define GRID_H_

#include "GridCell.h"
#include <netcdfcpp.h>
#include <netcdf.h>
#include <vector>
// Netcdf Install: sudo apt-get install netcdf-bin libnetcdfc++4 libnetcdfc7 libnetcdf-dev

class Grid
{
public:
	// Three public member functions.

	// Constructor sets the member variables to their given values
	// and calls the private member setInitialConditions().
	Grid (double height, double width, double resolution,
			double upperLeftLatitude, double U, double depth,
			double simulationSeconds, double writeInterval, double mountHeight);

	// The main simulation function, called once in the main program.
	void simulate ();

	// Destructor closes the nc-file properly.
	~Grid ();

private:
	double height;  // Grid size in North-South direction (km).
	double width;  // Grid size in East-West direction (km).
	double resolution;  // Grid cell width (km).
	double upperLeftLatitude;  // Latitude of North-West corner of the grid.
	double simulationLength; // Length of simulation in seconds.
	double U; // Speed of inflow at the western boundary (m/s).
	double depth; // Initial total depth of the atmosphere at the northern boundary (m)
	double dt; // Timestep defined such that Courant condition (U*dt/dx) = 0.01
	double courantNumber; // Courant number for this simulation.
	double writeInterval; // When to write the nc-file.
	double mountainHeight; // Height of the mountain.

	NcFile* ncfile; // Output file.
	NcVar *t, *u, *v, *h, *htot, *dx, *lat; // Output variables.

	// 2-dimensional vector, where the Grid Cells are stored.
	std::vector<std::vector<GridCell> > gridCells;

	// Initial condition setter.
	void setInitialConditions ();

	// Initialize non-boundary cells
	void initializeCenterCells ();

	// Take special care with initializing the boundaries.
	void initializeBoundaryCells ();

	// A short function that computes the latitude of each cell.
	double computeCellLatitude (int i);

	// Compute the height of geopotential isosurface as a function of latitude.
	double computeHtotCorrection (double lat);

	// Write the netcdf output file.
	void writeVariablesToFile (long timestepNum, double elapsedSeconds);

	// Step one timestep forward using the MacCormack method.
	void advanceGrid (long stepNum);
};

#endif /* GRID_H_ */
