#ifndef GRIDCELL_H_
#define GRIDCELL_H_

class GridCell
{
public:
	/* Constructor GridCell sets the initial conditions except the neighboring cells.
	 * All parameters in units described below and latitude in DEGREES!
	 */
	GridCell (double u, double v, double h, double hs, double latitude,
			double width, double dt);

	// Getters
	double getV () const;
	double getU () const;
	double getH () const;
	double getHtot () const;

	// Setters
	void setNorthernNeighbour (GridCell* c);
	void setEasternNeighbour (GridCell* c);
	void setSouthernNeighbour (GridCell* c);
	void setWesternNeighbour (GridCell* c);
	void setBoundaryType (char boundaryType);

	void computePredictor (long stepNum);

private:
	double u;  // Eastern wind component (m/s).
	double v;  // Northern wind component (m/s).
	double h;  // Depth of the "atmosphere" (m).
	double hs;  // Height of surface (m);
	double htot;  // Total height = h + hs (m).
	double f;  // Coriolis parameter = 2 * omega * sin(latitude) (1/s).
	double g;  // Gravitational acceleration (m/s^2).
	double width; // Cell width
	double dt; // Timestep
	double duPredictor, dvPredictor, dhPredictor;
	double uPredictor, vPredictor, hPredictor, hTotPredictor;
	char boundaryType;

	struct stencil
	{
		double current;
		double east;
		double west;
		double north;
		double south;
	};

	void giveInitialStencil (stencil& us, stencil& vs, stencil& hs,
			stencil& hts, long stepNum);

	// Pointers to neighboring cells.
	GridCell* northernNeighbour;
	GridCell* easternNeighbour;
	GridCell* southernNeighbour;
	GridCell* westernNeighbour;
};

#endif /* GRIDCELL_H_ */
