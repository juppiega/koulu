#ifndef GRIDCELL_H_
#define GRIDCELL_H_

class GridCell
{
public:
	/* Constructor GridCell sets the initial conditions except the neighboring cells.
	 * All parameters in units described below and latitude in DEGREES!
	 */
	GridCell (double u, double v, double htot, double latitude, double width,
			double dt);

	// Getters
	double getV () const;
	double getU () const;
	double getH () const;
	double getHtot () const;
	double getUPredictor () const;
	double getVPredictor () const;
	double getHPredictor () const;
	double getHtotPredictor () const;

	// Setters
	void setNorthernNeighbour (GridCell* c);
	void setEasternNeighbour (GridCell* c);
	void setSouthernNeighbour (GridCell* c);
	void setWesternNeighbour (GridCell* c);
	void setBoundaryType (char boundaryType);
	void setSurfaceHeight (double hs);
	void setTimestep (double dt);

	// Time marching (advancing) methods.
	void computePredictor (long stepNum);
	void computeCorrector (long stepNum);

	// Boundary condition declarations.
	static double westU, westV;
	static double northH, northV;
	static double southV;

private:
	double u;  // Eastern wind component (m/s).
	double v;  // Northern wind component (m/s).
	double h;  // Depth of the "atmosphere" (m).
	double hsurf;  // Height of surface (m);
	double htot;  // Total height = h + hs (m).
	double f;  // Coriolis parameter = 2 * omega * sin(latitude) (1/s).
	double g;  // Gravitational acceleration (m/s^2).
	double width; // Cell width
	double dt; // Timestep
	double duPredictor, dvPredictor, dhPredictor; // Predictor time derivatives
	double uPredictor, vPredictor, hPredictor, htotPredictor; // Predictor values
	char boundaryType;

	// Pointers to neighboring cells.
	GridCell* northernNeighbour;
	GridCell* easternNeighbour;
	GridCell* southernNeighbour;
	GridCell* westernNeighbour;

	// Struct to store neighboring values for finite difference.
	struct stencil
	{
		double current;
		double east;
		double west;
		double north;
		double south;
	};

	// Struct to store the derivative information for either forward or backward Euler step.
	struct oneSidedDifferences
	{
		double du, dv, dh;
	};

	// Get neighboring cell values.
	void giveStencil (stencil& us, stencil& vs, stencil& hs,
			stencil& hts, char stepType);

	// Computes one sided derivatives using differences to directions given by xdir and ydir.
	oneSidedDifferences computeDerivatives (stencil& us, stencil& vs,
			stencil& hs, stencil& hts, char xdir, char ydir);
};

#endif /* GRIDCELL_H_ */
