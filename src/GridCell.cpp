#include "GridCell.h"
#include <math.h>
#include <stdio.h>
#include <algorithm>

// Boundary condition definitions.
double GridCell::northH;
double GridCell::northV;
double GridCell::southV;
double GridCell::westU;
double GridCell::westV;

// CONSTRUCTOR ****************************************************
GridCell::GridCell (double u, double v, double htot, double latitude,
		double width, double dt)
{
	// Variables set in the constructor.
	this->u = u;
	this->v = v;
	this->htot = htot;
	this->f = 2 * 7.2921150E-5 * sin(latitude * M_PI / 180.0);
	this->g = 9.80665;
	this->width = width * 1000;
	this->dt = dt;

	// Variables modified using setters. These are set later.
	northernNeighbour = 0;
	easternNeighbour = 0;
	southernNeighbour = 0;
	westernNeighbour = 0;

	boundaryType = '0';

	hsurf = 0.0;
	this->h = htot - hsurf;

	duPredictor = dvPredictor = dhPredictor = 0.0;
	uPredictor = u;
	vPredictor = v;
	hPredictor = h;
	htotPredictor = htot;
}
// ****************************************************************

// GETTERS ********************************************************
double GridCell::getU () const
{
	return u;
}

double GridCell::getV () const
{
	return v;
}

double GridCell::getH () const
{
	return h;
}

double GridCell::getHtot () const
{
	return htot;
}

double GridCell::getUPredictor () const
{
	return uPredictor;
}

double GridCell::getVPredictor () const
{
	return vPredictor;
}

double GridCell::getHPredictor () const
{
	return hPredictor;
}

double GridCell::getHtotPredictor () const
{
	return htotPredictor;
}
// ****************************************************************

// SETTERS ********************************************************
void GridCell::setNorthernNeighbour (GridCell* c)
{
	northernNeighbour = c;
}

void GridCell::setEasternNeighbour (GridCell* c)
{
	easternNeighbour = c;
}

void GridCell::setSouthernNeighbour (GridCell* c)
{
	southernNeighbour = c;
}

void GridCell::setWesternNeighbour (GridCell* c)
{
	westernNeighbour = c;
}

void GridCell::setBoundaryType (char boundaryType)
{
	this->boundaryType = boundaryType;
}

void GridCell::setSurfaceHeight (double hs)
{
	this->hsurf = hs;
	h = htot - hsurf;
}

void GridCell::setTimestep (double dt)
{
	this->dt = dt;
}
// ****************************************************************

// PUBLIC TIME MARCHING METHODS ***********************************
void GridCell::computePredictor (long stepNum)
{
	// Values of neighboring cells. us = "u stencil", vs = "v stencil" etc.
	stencil us, vs, hs, hts;
	// Call the giveStencil function with predictor setting 'p'.
	giveStencil(us, vs, hs, hts, 'p');

	// Struct to store the predictor derivatives.
	oneSidedDifferences predictorStep;

	// Compute predictor in a cyclic manner, alternating the directions of
	// the x and y differences between forward and backward in a cycle of four.
	switch (stepNum % 4)
	{
	case 0:
		// Now, both derivatives are approximated using forward differences.
		predictorStep = computeDerivatives(us, vs, hs, hts, 'f', 'f');
		break;
	case 1:
		// Now, use forward difference for x and backward for y.
		predictorStep = computeDerivatives(us, vs, hs, hts, 'f', 'b');
		break;
	case 2:
		// Backward differences for both.
		predictorStep = computeDerivatives(us, vs, hs, hts, 'b', 'b');
		break;
	case 3:
		// Backward for x and forward for y.
		predictorStep = computeDerivatives(us, vs, hs, hts, 'b', 'f');
		break;
	}

	// Store the predictor derivatives.
	duPredictor = predictorStep.du;
	dvPredictor = predictorStep.dv;
	dhPredictor = predictorStep.dh;

	// Compute predictor values.
	uPredictor = u + duPredictor * dt;
	vPredictor = v + dvPredictor * dt;
	hPredictor = h + dhPredictor * dt;
	htotPredictor = hsurf + hPredictor;
}

void GridCell::computeCorrector (long stepNum)
{
	// Values of neighboring cells. us = "u stencil", vs = "v stencil" etc.
	stencil us, vs, hs, hts;
	// Call the giveStencil function with corrector setting 'c'.
	giveStencil(us, vs, hs, hts, 'c');

	// Struct to store the corrector derivatives.
	oneSidedDifferences correctorStep;

	// Compute corrector in a cyclic manner, alternating the directions of
	// the x and y differences between forward and backward in a cycle of four.
	// The directions are always exactly the opposite as in the predictor computation.
	// This allows for second order accuracy in time and space.
	switch (stepNum % 4)
	{
	case 0:
		correctorStep = computeDerivatives(us, vs, hs, hts, 'b', 'b');
		break;
	case 1:
		correctorStep = computeDerivatives(us, vs, hs, hts, 'b', 'f');
		break;
	case 2:
		correctorStep = computeDerivatives(us, vs, hs, hts, 'f', 'f');
		break;
	case 3:
		correctorStep = computeDerivatives(us, vs, hs, hts, 'f', 'b');
		break;
	}

	double duCorrector = correctorStep.du;
	double dvCorrector = correctorStep.dv;
	double dhCorrector = correctorStep.dh;

	// Compute the new values.
	u = u + 0.5 * (duCorrector + duPredictor) * dt;
	v = v + 0.5 * (dvCorrector + dvPredictor) * dt;
	h = h + 0.5 * (dhCorrector + dhPredictor) * dt;
	if (u < 0) // Clamp slightly such that u is never negative.
	{
		u = 0;
	}
	htot = h + hsurf;
}
// ****************************************************************

// PRIVATE UTILITY FUNCTIONS **************************************
GridCell::oneSidedDifferences GridCell::computeDerivatives (stencil& us,
		stencil& vs, stencil& hs, stencil& hts, char xdir, char ydir)
{
	oneSidedDifferences results; // Data to be returned.

	double dx = width;
	double dy = dx;
	double dhtotDx, dudx, dvdx, dhdx; // x derivatives.
	double dhtotDy, dudy, dvdy, dhdy; // y derivatives.

	if (xdir == 'f')
	{
		// Use forward difference (downwind to east).
		dhtotDx = (hts.east - hts.current) / dx;
		dudx = (us.east - us.current) / dx;
		dvdx = (vs.east - vs.current) / dx;
		dhdx = (hs.east - hs.current) / dx;
	}
	else
	{
		// Use backward difference (upwind to west).
		dhtotDx = (hts.current - hts.west) / dx;
		dudx = (us.current - us.west) / dx;
		dvdx = (vs.current - vs.west) / dx;
		dhdx = (hs.current - hs.west) / dx;
	}

	if (ydir == 'f')
	{
		// Forward difference.
		dhtotDy = (hts.north - hts.current) / dy;
		dudy = (us.north - us.current) / dy;
		dvdy = (vs.north - vs.current) / dy;
		dhdy = (hs.north - hs.current) / dy;
	}
	else
	{
		// Backward difference.
		dhtotDy = (hts.current - hts.south) / dy;
		dudy = (us.current - us.south) / dy;
		dvdy = (vs.current - vs.south) / dy;
		dhdy = (hs.current - hs.south) / dy;
	}

	// THE SHALLOW WATER EQUATIONS &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	results.du = -g * dhtotDx + f * vs.current - us.current * dudx //&
			- vs.current * dudy;								   //&
	results.dv = -g * dhtotDy - f * us.current - us.current * dvdx //&
			- vs.current * dvdy;                                   //&
	results.dh = -hs.current * (dudx + dvdy) - us.current * dhdx   //&
			- vs.current * dhdy;                                   //&
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	return results;
}

void GridCell::giveStencil (stencil& us, stencil& vs, stencil& hs, stencil& hts,
		char stepType)
{
	typedef double (GridCell::*ptrToGetFunction) () const;
	// Pointers used to call the getter functions.
	ptrToGetFunction getUFunc, getVFunc, getHFunc, getHtotFunc;

	if (stepType == 'p')
	{
		// Give the predictor variables.
		getUFunc = &GridCell::getU;
		getVFunc = &GridCell::getV;
		getHFunc = &GridCell::getH;
		getHtotFunc = &GridCell::getHtot;
	}
	else
	{
		// Give the corrector variables.
		getUFunc = &GridCell::getUPredictor;
		getVFunc = &GridCell::getVPredictor;
		getHFunc = &GridCell::getHPredictor;
		getHtotFunc = &GridCell::getHtotPredictor;
	}

	// Get values of this cell.
	us.current = (this->*getUFunc)();
	vs.current = (this->*getVFunc)();
	hs.current = (this->*getHFunc)();
	hts.current = (this->*getHtotFunc)();

	// Northern boundary conditions
	if (boundaryType != 'N' && boundaryType != '1' && boundaryType != '2')
	{
		// If NOT on the Northern border:
		us.north = (northernNeighbour->*getUFunc)();
		vs.north = (northernNeighbour->*getVFunc)();
		hs.north = (northernNeighbour->*getHFunc)();
		hts.north = (northernNeighbour->*getHtotFunc)();
	}
	else
	{
		// If on the Northern border (northernNeighbour == NULL), use boundary conditions:
		us.north = us.current;
		vs.north = northV;
		hs.north = northH;
		hts.north = hsurf + northH;
	}

	// Eastern BCs.
	if (boundaryType != 'E' && boundaryType != '2' && boundaryType != '4')
	{
		us.east = (easternNeighbour->*getUFunc)();
		vs.east = (easternNeighbour->*getVFunc)();
		hs.east = (easternNeighbour->*getHFunc)();
		hts.east = (easternNeighbour->*getHtotFunc)();
	}
	else
	{
		us.east = us.current;
		vs.east = vs.current;
		hs.east = hs.current;
		hts.east = hts.current;
	}

	// Southern BCs.
	if (boundaryType != 'S' && boundaryType != '3' && boundaryType != '4')
	{
		us.south = (southernNeighbour->*getUFunc)();
		vs.south = (southernNeighbour->*getVFunc)();
		hs.south = (southernNeighbour->*getHFunc)();
		hts.south = (southernNeighbour->*getHtotFunc)();
	}
	else
	{
		us.south = us.current;
		vs.south = southV;
		hs.south = hs.current;
		hts.south = hts.current;
	}

	// Western BCs.
	if (boundaryType != 'W' && boundaryType != '1' && boundaryType != '3')
	{
		us.west = (westernNeighbour->*getUFunc)();
		vs.west = (westernNeighbour->*getVFunc)();
		hs.west = (westernNeighbour->*getHFunc)();
		hts.west = (westernNeighbour->*getHtotFunc)();
	}
	else
	{
		us.west = westU;
		vs.west = westV;
		hs.west = hs.current;
		hts.west = hts.current;
	}
}
// END OF PRIVATE FUNCTIONS ***************************************
