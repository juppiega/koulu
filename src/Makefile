CC = g++
CFLAGS = -O3 -fopenmp -lnetcdf_c++

TARGET = shallowWater2d

all: $(TARGET)

$(TARGET): GridCell.cpp  GridCell.h  Grid.cpp  Grid.h  main.cpp
	$(CC) -o $(TARGET) main.cpp Grid.cpp GridCell.cpp $(CFLAGS) 

clean:
	rm -f $(TARGET)
