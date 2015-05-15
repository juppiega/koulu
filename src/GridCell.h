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
	double getVPrevious () const;
	double getUPrevious () const;
	double getHPrevious () const;

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

	// Clamper
	void clampSolution (long stepNum);

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
	double duPredictor, dvPredictor, dhPredictor;
	double uPredictor, vPredictor, hPredictor, htotPredictor;
	double uPrevious, vPrevious, hPrevious;
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

	void giveInitialStencil (stencil& us, stencil& vs, stencil& hs,
			stencil& hts);
	void givePredictorStencil (stencil& us, stencil& vs, stencil& hs,
			stencil& hts);
	void givePreviousStencil (stencil& us, stencil& vs, stencil& hs);
	oneSidedDifferences stepForward (stencil& us, stencil& vs, stencil& hs,
			stencil& hts);
	oneSidedDifferences stepBackward (stencil& us, stencil& vs, stencil& hs,
			stencil& hts);
	oneSidedDifferences computeDerivatives (stencil& us, stencil& vs,
			stencil& hs, stencil& hts, char xdir, char ydir);

	void clampIfNeeded (double uXPredDir, double uYPredDir, double vXPredDir,
			double vYPredDir, double hXPredDir, double hYPredDir);
};

#endif /* GRIDCELL_H_ */
