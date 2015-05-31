#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <math.h>
#include "Grid.h"

// Function for parsing command line arguments.
void readCmdArgs (int argc, char** argv, double& NSheight, double& width,
		double& resolution, double& upperLeftLatitude, double& U, double& depth,
		double& simSecs, double& writeInterval, double& mountHeight);

int main (int argc, char** argv)
{
	// Declare input variables.
	double NSheight, width, resolution, upperLeftLatitude, U, depth, simSecs,
			writeInterval, mountHeight;

	// Read command line arguments.
	readCmdArgs(argc, argv, NSheight, width, resolution, upperLeftLatitude, U,
			depth, simSecs, writeInterval, mountHeight);

	// Initialize the simulation grid.
	Grid grid(NSheight, width, resolution, upperLeftLatitude, U, depth, simSecs,
			writeInterval, mountHeight);

	// Run the simulation.
	grid.simulate();

	return 0;
}

void readCmdArgs (int argc, char** argv, double& NSheight, double& width,
		double& resolution, double& upperLeftLatitude, double& U, double& depth,
		double& simSecs, double& writeInterval, double& mountHeight)
{
	NSheight = 5000; // in km.
	width = 10000; // in km.
	resolution = 15; // in km.
	upperLeftLatitude = 70; // in degrees.
	U = 5; // in m/s.
	depth = 10000; // in m.
	simSecs = 691200; // = 8 days
	writeInterval = 7200; // = 2 hours
	mountHeight = 5000; // in m.

	// Check argument count.
	if (argc % 2 == 0)
	{
		throw std::runtime_error(
				std::string("Number of arguments must be even (or zero)!\n")
						+ "Usage: shallowWater [--<argName1> <argValue1> ...]\n"
						+ "See the project pdf for more information.");
	}
	// Read in arguments.
	for (int i = 1; i < argc; i += 2)
	{
		std::string argStr(argv[i]);
		double argVal = atof(argv[i + 1]);

		if (argStr.compare("--n-s-height") == 0) // <str1>.compare(<str2>) returns 0 for equal strings.
		{
			NSheight = argVal;
			if (NSheight < 2000)
			{
				throw std::runtime_error(
						std::string(
								"Height (N-S direction) must be greater than 2000 km."));
			}
		}
		else if (argStr.compare("--resolution") == 0)
		{
			resolution = argVal;
		}
		else if (argStr.compare("--latitude") == 0)
		{
			upperLeftLatitude = argVal;
		}
		else if (argStr.compare("--U") == 0)
		{
			U = argVal;
		}
		else if (argStr.compare("--depth") == 0)
		{
			depth = argVal;
		}
		else if (argStr.compare("--duration") == 0)
		{
			simSecs = argVal * 86400;
		}
		else if (argStr.compare("--write-interval") == 0)
		{
			writeInterval = argVal * 3600;
		}
		else if (argStr.compare("--mountain-height") == 0)
		{
			mountHeight = argVal;
		}
	}
	if (depth < mountHeight)
	{
		throw std::runtime_error(
				std::string("Depth must be greater than the mountain height."));
	}
	if (resolution > NSheight / 8)
	{
		throw std::runtime_error(
				std::string(
						"Resolution must be less than one eighth of grid height."));
	}

	// Compute grid width in West to East direction.
	width = 2.0 * M_PI
			* sqrt(
					U
							/ (2.0 * 7.2921150E-5
									* cos(upperLeftLatitude * M_PI / 180.0)
									/ 6371E3)) / 1000 + 5000;

	printf("Using values:\n");
	printf("N-S Height: %.1f km\n", NSheight);
	printf("W-E width: %.1f km\n", width);
	printf("Resolution: %.1f km\n", resolution);
	printf("Latitude: %.1f\n", upperLeftLatitude);
	printf("U: %.1f m/s\n", U);
	printf("Depth: %.1f m\n", depth);
	printf("Simulation duration: %.2f d\n", simSecs / 86400);
	printf("Write Interval: %.2f h\n", writeInterval / 3600);
	printf("Mountain height: %.1f m\n\n", mountHeight);
}
