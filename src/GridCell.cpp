#include "GridCell.h"
#include <math.h>

GridCell::GridCell (double u, double v, double h, double hs, double latitude)
{
	this->u = u;
	this->v = v;
	this->h = h;
	this->hs = hs;
	this->htot = h + hs;
	this->f = 2 * 7.2921150E-5 * sin(latitude * M_PI / 180.0);
	this->g = 9.80665;
}

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

