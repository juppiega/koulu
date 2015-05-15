#include "GridCell.h"
#include <math.h>
#include <stdio.h>
#include <algorithm>

using std::min;
using std::max;

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

	// Variables modified using setters.
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

	uPrevious = u;
	vPrevious = v;
	hPrevious = h;
}

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

double GridCell::getUPrevious () const
{
	return uPrevious;
}

double GridCell::getVPrevious () const
{
	return vPrevious;
}

double GridCell::getHPrevious () const
{
	return hPrevious;
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
	// The Western boundary won't change.
	if (boundaryType == 'W' || boundaryType == '1' || boundaryType == '3')
	{
		return;
	}

	stencil us, vs, hs, hts;
	giveInitialStencil(us, vs, hs, hts);

	oneSidedDifferences predictorStep;

	switch (stepNum % 4)
	{
	case 0:
		predictorStep = computeDerivatives(us, vs, hs, hts, 'f', 'f');
		break;
	case 1:
		predictorStep = computeDerivatives(us, vs, hs, hts, 'f', 'b');
		break;
	case 2:
		predictorStep = computeDerivatives(us, vs, hs, hts, 'b', 'b');
		break;
	case 3:
		predictorStep = computeDerivatives(us, vs, hs, hts, 'b', 'f');
		break;
	}

	duPredictor = predictorStep.du;
	dvPredictor = predictorStep.dv;
	dhPredictor = predictorStep.dh;

	uPredictor = u + duPredictor * dt;
	vPredictor = v + dvPredictor * dt;
	hPredictor = h + dhPredictor * dt;
	htotPredictor = hsurf + hPredictor;
}

void GridCell::computeCorrector (long stepNum)
{
	// The Western boundary won't change.
	if (boundaryType == 'W' || boundaryType == '1' || boundaryType == '3')
	{
		return;
	}

	stencil us, vs, hs, hts;
	givePredictorStencil(us, vs, hs, hts);

	oneSidedDifferences correctorStep;

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

	uPrevious = u;
	vPrevious = v;
	hPrevious = h;
	u = u + 0.5 * (duCorrector + duPredictor) * dt;
	v = v + 0.5 * (dvCorrector + dvPredictor) * dt;
	h = h + 0.5 * (dhCorrector + dhPredictor) * dt;
	htot = h + hsurf;
}
// ****************************************************************

// CLAMPER ********************************************************
void GridCell::clampSolution (long stepNum)
{
	// The Western boundary won't change.
	if (boundaryType == 'W' || boundaryType == '1' || boundaryType == '3')
	{
		return;
	}

	stencil uPrev, vPrev, hPrev;
	givePreviousStencil(uPrev, vPrev, hPrev);

	switch (stepNum % 4)
	{
	case 0:
		clampIfNeeded(uPrev.east, uPrev.north, vPrev.east, vPrev.north,
				hPrev.east, hPrev.north);
		break;
	case 1:
		clampIfNeeded(uPrev.east, uPrev.south, vPrev.east, vPrev.south,
				hPrev.east, hPrev.south);
		break;
	case 2:
		clampIfNeeded(uPrev.west, uPrev.south, vPrev.west, vPrev.south,
				hPrev.west, hPrev.south);
		break;
	case 3:
		clampIfNeeded(uPrev.west, uPrev.north, vPrev.west, vPrev.north,
				hPrev.west, hPrev.north);
		break;
	}
}
// ****************************************************************

// PRIVATE UTILITY FUNCTIONS **************************************
GridCell::oneSidedDifferences GridCell::computeDerivatives (stencil& us,
		stencil& vs, stencil& hs, stencil& hts, char xdir, char ydir)
{
	oneSidedDifferences results;

	double dx = width;
	double dy = dx;
	double dhtotDx, dudx, dvdx, dhdx;
	double dhtotDy, dudy, dvdy, dhdy;

	if (xdir == 'f')
	{
		dhtotDx = (hts.east - hts.current) / dx;
		dudx = (us.east - us.current) / dx;
		dvdx = (vs.east - vs.current) / dx;
		dhdx = (hs.east - hs.current) / dx;
	}
	else
	{
		dhtotDx = (hts.current - hts.west) / dx;
		dudx = (us.current - us.west) / dx;
		dvdx = (vs.current - vs.west) / dx;
		dhdx = (hs.current - hs.west) / dx;
	}

	if (ydir == 'f')
	{
		dhtotDy = (hts.north - hts.current) / dy;
		dudy = (us.north - us.current) / dy;
		dvdy = (vs.north - vs.current) / dy;
		dhdy = (hs.north - hs.current) / dy;
	}
	else
	{
		dhtotDy = (hts.current - hts.south) / dy;
		dudy = (us.current - us.south) / dy;
		dvdy = (vs.current - vs.south) / dy;
		dhdy = (hs.current - hs.south) / dy;
	}

	results.du = -g * dhtotDx + f * vs.current - us.current * dudx
			- vs.current * dudy;
	results.dv = -g * dhtotDy - f * us.current - us.current * dvdx
			- vs.current * dvdy;
	results.dh = -hs.current * (dudx + dvdy) - us.current * dhdx
			- vs.current * dhdy;

//	if (results.du > 0.1)
//	{
//		printf("YOU FAILED!\n");
//	}

	return results;
}

void GridCell::clampIfNeeded (double uXPredDir, double uYPredDir,
		double vXPredDir, double vYPredDir, double hXPredDir, double hYPredDir)
{
	if (u < min( { uPrevious, uXPredDir, uYPredDir }) || u > max( { uPrevious,
			uXPredDir, uYPredDir }))
	{
		u = uPredictor;
		//printf("CLAMPED!\n");
	}
	if (v < min( { vPrevious, vXPredDir, vYPredDir }) || v > max( { vPrevious,
			vXPredDir, vYPredDir }))
	{
		v = vPredictor;
		//printf("CLAMPED!\n");
	}
	if (h < min( { hPrevious, hXPredDir, hYPredDir }) || h > max( { hPrevious,
			hXPredDir, hYPredDir }))
	{
		h = hPredictor;
		htot = h + hsurf;
		//printf("CLAMPED!\n");
	}

//	double uMin = min( { uPrevious, uXPredDir, uYPredDir });
//	double uMax = max( { uPrevious, uXPredDir, uYPredDir });
//	double vMin = min( { vPrevious, vXPredDir, vYPredDir });
//	double vMax = max( { vPrevious, vXPredDir, vYPredDir });
//	double hMin = min( { hPrevious, hXPredDir, hYPredDir });
//	double hMax = max( { hPrevious, hXPredDir, hYPredDir });
//
//	if (u < uMin)
//		u = uMin;
//	else if (u > uMax)
//		u = uMax;
//
//	if (v < vMin)
//		v = vMin;
//	else if (v > vMax)
//		v = vMax;
//
//	if (h < hMin)
//		h = hMin;
//	else if (h > hMax)
//	{
//		h = hMax;
//		htot = h + hsurf;
//	}

}

void GridCell::giveInitialStencil (stencil& us, stencil& vs, stencil& hs,
		stencil& hts)
{
	us.current = u;
	vs.current = v;
	hs.current = h;
	hts.current = htot;

	us.west = westernNeighbour->getU();
	vs.west = westernNeighbour->getV();
	hs.west = westernNeighbour->getH();
	hts.west = westernNeighbour->getHtot();

	if (boundaryType == '0')
	{
		us.north = northernNeighbour->getU();
		us.east = easternNeighbour->getU();
		vs.north = northernNeighbour->getV();
		vs.east = easternNeighbour->getV();
		hs.north = northernNeighbour->getH();
		hs.east = easternNeighbour->getH();
		hts.north = northernNeighbour->getHtot();
		hts.east = easternNeighbour->getHtot();

		us.south = southernNeighbour->getU();
		vs.south = southernNeighbour->getV();
		hs.south = southernNeighbour->getH();
		hts.south = southernNeighbour->getHtot();
	}
	else if (boundaryType == 'N')
	{
		us.south = southernNeighbour->getU();
		vs.south = southernNeighbour->getV();
		hs.south = southernNeighbour->getH();
		hts.south = southernNeighbour->getHtot();

		us.east = easternNeighbour->getU();
		vs.east = easternNeighbour->getV();
		hs.east = easternNeighbour->getH();
		hts.east = easternNeighbour->getHtot();

		us.north = us.current; //+ (us.current - us.south);
		vs.north = vs.current; // + (vs.current - vs.south);
		hs.north = hs.current; // + (hs.current - hs.south);
		hts.north = hts.current; // + (hts.current - hts.south);
	}
	else if (boundaryType == 'E')
	{
		us.north = northernNeighbour->getU();
		vs.north = northernNeighbour->getV();
		hs.north = northernNeighbour->getH();
		hts.north = northernNeighbour->getHtot();

		us.east = us.current; // + (us.current - us.west);
		vs.east = vs.current; // + (vs.current - vs.west);
		hs.east = hs.current; // + (hs.current - hs.west);
		hts.east = hts.current; // + (hts.current - hts.west);

		us.south = southernNeighbour->getU();
		vs.south = southernNeighbour->getV();
		hs.south = southernNeighbour->getH();
		hts.south = southernNeighbour->getHtot();
	}
	else if (boundaryType == 'S')
	{
		us.north = northernNeighbour->getU();
		vs.north = northernNeighbour->getV();
		hs.north = northernNeighbour->getH();
		hts.north = northernNeighbour->getHtot();

		us.east = easternNeighbour->getU();
		vs.east = easternNeighbour->getV();
		hs.east = easternNeighbour->getH();
		hts.east = easternNeighbour->getHtot();

		us.south = us.current; // - (us.north - us.current);
		vs.south = vs.current; // - (vs.north - vs.current);
		hs.south = hs.current; // - (hs.north - hs.current);
		hts.south = hts.current; // - (hts.north - hts.current);

	}
	else if (boundaryType == '2')
	{
		us.south = southernNeighbour->getU();
		vs.south = southernNeighbour->getV();
		hs.south = southernNeighbour->getH();
		hts.south = southernNeighbour->getHtot();

		us.north = us.current; // + (us.current - us.south);
		vs.north = vs.current; // + (vs.current - vs.south);
		hs.north = hs.current; // + (hs.current - hs.south);
		hts.north = hts.current; // + (hts.current - hts.south);

		us.east = us.current; // + (us.current - us.west);
		vs.east = vs.current; // + (vs.current - vs.west);
		hs.east = hs.current; // + (hs.current - hs.west);
		hts.east = hts.current; // + (hts.current - hts.west);
	}
	else if (boundaryType == '4')
	{
		us.north = northernNeighbour->getU();
		vs.north = northernNeighbour->getV();
		hs.north = northernNeighbour->getH();
		hts.north = northernNeighbour->getHtot();

		us.east = us.current; // + (us.current - us.west);
		vs.east = vs.current; // + (vs.current - vs.west);
		hs.east = hs.current; // + (hs.current - hs.west);
		hts.east = hts.current; // + (hts.current - hts.west);

		us.south = us.current; // - (us.north - us.current);
		vs.south = vs.current; // - (vs.north - vs.current);
		hs.south = hs.current; // - (hs.north - hs.current);
		hts.south = hts.current; // - (hts.north - hts.current);
	}
}

void GridCell::givePredictorStencil (stencil& us, stencil& vs, stencil& hs,
		stencil& hts)
{
	us.current = uPredictor;
	vs.current = vPredictor;
	hs.current = hPredictor;
	hts.current = htotPredictor;

	us.west = westernNeighbour->getUPredictor();
	vs.west = westernNeighbour->getVPredictor();
	hs.west = westernNeighbour->getHPredictor();
	hts.west = westernNeighbour->getHtotPredictor();

	if (boundaryType == '0')
	{
		us.north = northernNeighbour->getUPredictor();
		us.east = easternNeighbour->getUPredictor();
		vs.north = northernNeighbour->getVPredictor();
		vs.east = easternNeighbour->getVPredictor();
		hs.north = northernNeighbour->getHPredictor();
		hs.east = easternNeighbour->getHPredictor();
		hts.north = northernNeighbour->getHtotPredictor();
		hts.east = easternNeighbour->getHtotPredictor();

		us.south = southernNeighbour->getUPredictor();
		vs.south = southernNeighbour->getVPredictor();
		hs.south = southernNeighbour->getHPredictor();
		hts.south = southernNeighbour->getHtotPredictor();
	}
	else if (boundaryType == 'N')
	{
		us.south = southernNeighbour->getUPredictor();
		vs.south = southernNeighbour->getVPredictor();
		hs.south = southernNeighbour->getHPredictor();
		hts.south = southernNeighbour->getHtotPredictor();

		us.east = easternNeighbour->getUPredictor();
		vs.east = easternNeighbour->getVPredictor();
		hs.east = easternNeighbour->getHPredictor();
		hts.east = easternNeighbour->getHtotPredictor();

		us.north = us.current; // + (us.current - us.south);
		vs.north = vs.current; // + (vs.current - vs.south);
		hs.north = hs.current; // + (hs.current - hs.south);
		hts.north = hts.current; // + (hts.current - hts.south);
	}
	else if (boundaryType == 'E')
	{
		us.north = northernNeighbour->getUPredictor();
		vs.north = northernNeighbour->getVPredictor();
		hs.north = northernNeighbour->getHPredictor();
		hts.north = northernNeighbour->getHtotPredictor();

		us.east = us.current; // + (us.current - us.west);
		vs.east = vs.current; // + (vs.current - vs.west);
		hs.east = hs.current; // + (hs.current - hs.west);
		hts.east = hts.current; // + (hts.current - hts.west);

		us.south = southernNeighbour->getUPredictor();
		vs.south = southernNeighbour->getVPredictor();
		hs.south = southernNeighbour->getHPredictor();
		hts.south = southernNeighbour->getHtotPredictor();
	}
	else if (boundaryType == 'S')
	{
		us.north = northernNeighbour->getUPredictor();
		vs.north = northernNeighbour->getVPredictor();
		hs.north = northernNeighbour->getHPredictor();
		hts.north = northernNeighbour->getHtotPredictor();

		us.east = easternNeighbour->getUPredictor();
		vs.east = easternNeighbour->getVPredictor();
		hs.east = easternNeighbour->getHPredictor();
		hts.east = easternNeighbour->getHtotPredictor();

		us.south = us.current; // - (us.north - us.current);
		vs.south = vs.current; // - (vs.north - vs.current);
		hs.south = hs.current; // - (hs.north - hs.current);
		hts.south = hts.current; // - (hts.north - hts.current);
	}
	else if (boundaryType == '2')
	{
		us.south = southernNeighbour->getUPredictor();
		vs.south = southernNeighbour->getVPredictor();
		hs.south = southernNeighbour->getHPredictor();
		hts.south = southernNeighbour->getHtotPredictor();

		us.north = us.current; // + (us.current - us.south);
		vs.north = vs.current; // + (vs.current - vs.south);
		hs.north = hs.current; // + (hs.current - hs.south);
		hts.north = hts.current; // + (hts.current - hts.south);

		us.east = us.current; // + (us.current - us.west);
		vs.east = vs.current; // + (vs.current - vs.west);
		hs.east = hs.current; // + (hs.current - hs.west);
		hts.east = hts.current; // + (hts.current - hts.west);
	}
	else if (boundaryType == '4')
	{
		us.north = northernNeighbour->getUPredictor();
		vs.north = northernNeighbour->getVPredictor();
		hs.north = northernNeighbour->getHPredictor();
		hts.north = northernNeighbour->getHtotPredictor();

		us.east = us.current; // + (us.current - us.west);
		vs.east = vs.current; // + (vs.current - vs.west);
		hs.east = hs.current; // + (hs.current - hs.west);
		hts.east = hts.current; // + (hts.current - hts.west);

		us.south = us.current; // - (us.north - us.current);
		vs.south = vs.current; // - (vs.north - vs.current);
		hs.south = hs.current; // - (hs.north - hs.current);
		hts.south = hts.current; // - (hts.north - hts.current);
	}
}

void GridCell::givePreviousStencil (stencil& us, stencil& vs, stencil& hs)
{
	us.west = westernNeighbour->getUPrevious();
	vs.west = westernNeighbour->getVPrevious();
	hs.west = westernNeighbour->getHPrevious();

	if (boundaryType == '0')
	{
		us.north = northernNeighbour->getUPrevious();
		us.east = easternNeighbour->getUPrevious();
		vs.north = northernNeighbour->getVPrevious();
		vs.east = easternNeighbour->getVPrevious();
		hs.north = northernNeighbour->getHPrevious();
		hs.east = easternNeighbour->getHPrevious();

		us.south = southernNeighbour->getUPrevious();
		vs.south = southernNeighbour->getVPrevious();
		hs.south = southernNeighbour->getHPrevious();
	}
	else if (boundaryType == 'N')
	{
		us.south = southernNeighbour->getUPrevious();
		vs.south = southernNeighbour->getVPrevious();
		hs.south = southernNeighbour->getHPrevious();

		us.east = easternNeighbour->getUPrevious();
		vs.east = easternNeighbour->getVPrevious();
		hs.east = easternNeighbour->getHPrevious();

		us.north = us.current; // + (us.current - us.south);
		vs.north = vs.current; // + (vs.current - vs.south);
		hs.north = hs.current; // + (hs.current - hs.south);
	}
	else if (boundaryType == 'E')
	{
		us.north = northernNeighbour->getUPrevious();
		vs.north = northernNeighbour->getVPrevious();
		hs.north = northernNeighbour->getHPrevious();

		us.east = us.current; // + (us.current - us.west);
		vs.east = vs.current; // + (vs.current - vs.west);
		hs.east = hs.current; // + (hs.current - hs.west);

		us.south = southernNeighbour->getUPrevious();
		vs.south = southernNeighbour->getVPrevious();
		hs.south = southernNeighbour->getHPrevious();
	}
	else if (boundaryType == 'S')
	{
		us.north = northernNeighbour->getUPrevious();
		vs.north = northernNeighbour->getVPrevious();
		hs.north = northernNeighbour->getHPrevious();

		us.east = easternNeighbour->getUPrevious();
		vs.east = easternNeighbour->getVPrevious();
		hs.east = easternNeighbour->getHPrevious();

		us.south = us.current; // - (us.north - us.current);
		vs.south = vs.current; // - (vs.north - vs.current);
		hs.south = hs.current; // - (hs.north - hs.current);
	}
	else if (boundaryType == '2')
	{
		us.south = southernNeighbour->getUPrevious();
		vs.south = southernNeighbour->getVPrevious();
		hs.south = southernNeighbour->getHPrevious();

		us.north = us.current; // + (us.current - us.south);
		vs.north = vs.current; // + (vs.current - vs.south);
		hs.north = hs.current; // + (hs.current - hs.south);

		us.east = us.current; // + (us.current - us.west);
		vs.east = vs.current; // + (vs.current - vs.west);
		hs.east = hs.current; // + (hs.current - hs.west);
	}
	else if (boundaryType == '4')
	{
		us.north = northernNeighbour->getUPrevious();
		vs.north = northernNeighbour->getVPrevious();
		hs.north = northernNeighbour->getHPrevious();

		us.east = us.current; // + (us.current - us.west);
		vs.east = vs.current; // + (vs.current - vs.west);
		hs.east = hs.current; // + (hs.current - hs.west);

		us.south = us.current; // - (us.north - us.current);
		vs.south = vs.current; // - (vs.north - vs.current);
		hs.south = hs.current; // - (hs.north - hs.current);
	}
}
// ****************************************************************
