#include <stdio.h>
#include <stdlib.h>
#include "Grid.h"

int main (int argc, char** argv)
{
	Grid grid(1000, 2000, 10, 45, 5, 10, 10000, 60);
	grid.simulate();
	return 0;
}
