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
	Grid (double height, double width, double resolution,
			double upperLeftLatitude, double U, double initHeight, double simulationSeconds);
	~Grid();

private:
	double height;  // Grid size in North-South direction (km).
	double width;  // Grid size in East-West direction (km).
	double resolution;  // Grid cell width (km).
	double upperLeftLatitude;  // Latitude of North-West corner of the grid.
	double simulationSeconds; // Length of simulation in seconds.
	double U; // Speed of inflow at the western boundaty (m/s).
	double initHeight; // Initial total height (km)
	double dt; // Timestep defined such that Courant condition (U*dt/dx) = 0.9
	long timestepNum;
	NcFile* ncfile; // Output file.
	NcVar *t, *u, *v, *h, *htot, *dx;

	// 2-dimensional vector, where the Grid Cells are stored.
	std::vector<std::vector<GridCell> > gridCells;

	// Initial condition setter.
	void setInitialConditions();

	// Write the netcdf
	void writeVariablesToFile();
};

#endif /* GRID_H_ */
