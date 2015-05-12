#include "GridCell.h"
#include <math.h>

GridCell::GridCell (double u, double v, double h, double hs, double latitude,
		double width, double dt)
{
	// Variables set in the constructor.
	this->u = u;
	this->v = v;
	this->h = h;
	this->hs = hs;
	this->htot = h + hs;
	this->f = 2 * 7.2921150E-5 * sin(latitude * M_PI / 180.0);
	this->g = 9.80665;
	this->width = width;
	this->dt = dt;

	duPredictor = dvPredictor = dhPredictor = 0.0;
	uPredictor = u;
	vPredictor = v;
	hPredictor = h;
	hTotPredictor = htot;

	// Variables modified using setters.
	northernNeighbour = 0;
	easternNeighbour = 0;
	southernNeighbour = 0;
	westernNeighbour = 0;

	boundaryType = '0';
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
// ****************************************************************

void GridCell::computePredictor (long stepNum)
{
	// The Western boundary won't change.
	if (boundaryType == 'W' || boundaryType == '1' || boundaryType == '3')
	{
		return;
	}

	stencil us, vs, hs, hts;

}

void GridCell::giveInitialStencil (stencil& us, stencil& vs, stencil& hs,
		stencil& hts, long stepNum)
{
	us.current = u;
	vs.current = v;
	hs.current = h;
	hts.current = htot;

	if (boundaryType == '0')
	{
		if (stepNum % 2 == 0)
		{
			us.north = northernNeighbour->getU();
			us.east = easternNeighbour->getU();
			vs.north = northernNeighbour->getV();
			vs.east = easternNeighbour->getV();
			hs.north = northernNeighbour->getH();
			hs.east = easternNeighbour->getH();
			hts.north = northernNeighbour->getHtot();
			hts.east = easternNeighbour->getHtot();
		}
		else
		{
			us.south = southernNeighbour->getU();
			us.west = westernNeighbour->getU();
			vs.south = southernNeighbour->getV();
			vs.west = westernNeighbour->getV();
			hs.south = southernNeighbour->getH();
			hs.west = westernNeighbour->getH();
			hts.south = southernNeighbour->getHtot();
			hts.west = westernNeighbour->getHtot();
		}
	}
	else if (boundaryType == 'N')
	{
		us.south = southernNeighbour->getU();
		vs.south = southernNeighbour->getV();
		hs.south = southernNeighbour->getH();
		hts.south = southernNeighbour->getHtot();

		if (stepNum % 2 == 0)
		{
			us.east = easternNeighbour->getU();
			vs.east = easternNeighbour->getV();
			hs.east = easternNeighbour->getH();
			hts.east = easternNeighbour->getHtot();

			us.north = us.current + (us.current - us.south);
			vs.north = vs.current + (vs.current - vs.south);
			hs.north = hs.current + (hs.current - hs.south);
			hts.north = hts.current + (hts.current - hts.south);
		}
		else
		{
			us.west = westernNeighbour->getU();
			vs.west = westernNeighbour->getV();
			hs.west = westernNeighbour->getH();
			hts.west = westernNeighbour->getHtot();
		}
	}
	else if (boundaryType == 'E')
	{
		us.west = westernNeighbour->getU();
		vs.west = westernNeighbour->getV();
		hs.west = westernNeighbour->getH();
		hts.west = westernNeighbour->getHtot();

		if (stepNum % 2 == 0)
		{
			us.north = northernNeighbour->getU();
			vs.north = northernNeighbour->getV();
			hs.north = northernNeighbour->getH();
			hts.north = northernNeighbour->getHtot();

			us.east = us.current + (us.current - us.west);
			vs.east = vs.current + (vs.current - vs.west);
			hs.east = hs.current + (hs.current - hs.west);
			hts.east = hts.current + (hts.current - hts.west);
		}
		else
		{
			us.south = southernNeighbour->getU();
			vs.south = southernNeighbour->getV();
			hs.south = southernNeighbour->getH();
			hts.south = southernNeighbour->getHtot();
		}
	}
	else if (boundaryType == 'S')
	{
		us.north = northernNeighbour->getU();
		vs.north = northernNeighbour->getV();
		hs.north = northernNeighbour->getH();
		hts.north = northernNeighbour->getHtot();

		if (stepNum % 2 == 0)
		{
			us.east = easternNeighbour->getU();
			vs.east = easternNeighbour->getV();
			hs.east = easternNeighbour->getH();
			hts.east = easternNeighbour->getHtot();
		}
		else
		{
			us.west = westernNeighbour->getU();
			vs.west = westernNeighbour->getV();
			hs.west = westernNeighbour->getH();
			hts.west = westernNeighbour->getHtot();

			us.south = us.current - (us.north - us.current);
			vs.south = vs.current - (vs.north - vs.current);
			hs.south = hs.current - (hs.north - hs.current);
			hts.south = hts.current - (hts.north - hts.current);
		}
	}
	else if (boundaryType == '2')
	{
		us.west = westernNeighbour->getU();
		vs.west = westernNeighbour->getV();
		hs.west = westernNeighbour->getH();
		hts.west = westernNeighbour->getHtot();
		us.south = southernNeighbour->getU();
		vs.south = southernNeighbour->getV();
		hs.south = southernNeighbour->getH();
		hts.south = southernNeighbour->getHtot();

		if (stepNum % 2 == 0)
		{
			us.north = us.current + (us.current - us.south);
			vs.north = vs.current + (vs.current - vs.south);
			hs.north = hs.current + (hs.current - hs.south);
			hts.north = hts.current + (hts.current - hts.south);

			us.east = us.current + (us.current - us.west);
			vs.east = vs.current + (vs.current - vs.west);
			hs.east = hs.current + (hs.current - hs.west);
			hts.east = hts.current + (hts.current - hts.west);
		}
	}
	else if (boundaryType == '3')
	{
		us.west = westernNeighbour->getU();
		vs.west = westernNeighbour->getV();
		hs.west = westernNeighbour->getH();
		hts.west = westernNeighbour->getHtot();
		us.north = northernNeighbour->getU();
		vs.north = northernNeighbour->getV();
		hs.north = northernNeighbour->getH();
		hts.north = northernNeighbour->getHtot();

		if (stepNum % 2 == 0)
		{
			us.east = us.current + (us.current - us.west);
			vs.east = vs.current + (vs.current - vs.west);
			hs.east = hs.current + (hs.current - hs.west);
			hts.east = hts.current + (hts.current - hts.west);
		}
		else
		{
			us.south = us.current - (us.north - us.current);
			vs.south = vs.current - (vs.north - vs.current);
			hs.south = hs.current - (hs.north - hs.current);
			hts.south = hts.current - (hts.north - hts.current);
		}
	}
}

