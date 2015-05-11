#ifndef GRIDCELL_H_
#define GRIDCELL_H_

class GridCell
{
public:
	/* Constructor GridCell sets the initial conditions except the neighboring cells.
	 * All parameters in units described below and latitude in DEGREES!
	 */
	GridCell (double u, double v, double h, double hs, double latitude);

	// Getters
	double getV () const;
	double getU () const;
	double getH () const;
	double getHtot () const;

private:
	double u;  // Eastern wind component (m/s).
	double v;  // Northern wind component (m/s).
	double h;  // Depth of the "atmosphere" (m).
	double hs;  // Height of surface (m);
	double htot;  // Total height = h + hs (m).
	double f;  // Coriolis parameter = 2 * omega * sin(latitude) (1/s).
	double g;  // Gravitational acceleration (m/s^2).
};

#endif /* GRIDCELL_H_ */
