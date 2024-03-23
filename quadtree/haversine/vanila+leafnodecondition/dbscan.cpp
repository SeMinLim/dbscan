#include <sys/resource.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>


// DBSCAN
#define CORE_POINT 1
#define BORDER_POINT 2

#define SUCCESS 0
#define UNCLASSIFIED -1
#define NOISE -2
#define FAILURE -3

#define MINIMUM_POINTS 2
#define PCIECONDITION 4

// Haversine
#define EARTH_RADIUS 6371
#define TO_RADIAN (3.1415926536 / 180)
#define TO_DEGREE (180 / 3.1415926536)


uint64_t numDataPoints = 0;
uint64_t numQuadrants = 0;
uint64_t numHaversine = 0;
uint64_t numPointsCondition = 0;
uint64_t numNormal = 0;
uint64_t numComplete = 0;


typedef struct Point {
	double lat, lon;
}Point;

typedef struct PointQuadTree {
	Point point;
	size_t datasetID;
}PointQuadTree;

typedef struct PointDBSCAN {
	Point point;
	Point northEastern;
	Point northWestern;
	Point southEastern;
	Point southWestern;
	size_t clusterID;
}PointDBSCAN;

typedef struct Quadrant {
	std::vector<Quadrant*> child;
	std::vector<PointQuadTree> cities;
	Point northEastern;
	Point northWestern;
	Point southEastern;
	Point southWestern;
	Point center;
	double diagonal;
	bool done;
	bool same;
}Quadrant;


std::vector<PointDBSCAN> dataset;


// Elapsed time checker
static inline double timeCheckerCPU(void) {
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

// Function for reading benchmark file
// Quadrant
void readBenchmarkDataQuadTree(Quadrant *root, char* filename, size_t length) {
	FILE *f_data = fopen(filename, "rb");
	if( f_data == NULL ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	for ( size_t i = 0; i < length; i ++ ) {
		PointQuadTree temp;
		size_t datasetID = i;
		float lat = 0.00;
		float lon = 0.00;
		fread(&lat, sizeof(float), 1, f_data);
		fread(&lon, sizeof(float), 1, f_data);
		temp.point.lat = (double)lat;
		temp.point.lon = (double)lon;
		temp.datasetID = datasetID;
		root->cities.push_back(temp);
	}

	root->done = 0;
	root->same = 0;
	
	fclose(f_data);
}
// Data points
void readBenchmarkDataDBSCAN(std::vector<PointDBSCAN> &dataset, char* filename, size_t length) {
	FILE *f_data = fopen(filename, "rb");
	if( f_data == NULL ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	for ( size_t i = 0; i < length; i ++ ) {
		PointDBSCAN temp;
		size_t clusterID = UNCLASSIFIED;
		float lat = 0.00;
		float lon = 0.00;
		fread(&lat, sizeof(float), 1, f_data);
		fread(&lon, sizeof(float), 1, f_data);
		temp.point.lat = (double)lat;
		temp.point.lon = (double)lon;
		temp.clusterID = clusterID;
		dataset.push_back(temp);
	}
	
	fclose(f_data);
}

// Function for writing benchmark file
void writeBenchmarkDataDBSCAN(std::vector<PointDBSCAN> &dataset, char* filename, size_t length) {
	FILE* f_data = fopen(filename, "wb");

	for ( size_t i = 0; i < length; i ++ ) {
		fwrite(&dataset[i].point.lat, sizeof(float), 1, f_data);
		fwrite(&dataset[i].point.lon, sizeof(float), 1, f_data);
		fwrite(&dataset[i].clusterID, sizeof(size_t), 1, f_data);
	}

	fclose(f_data);
}

// Haversine
double haversine(const Point pointCore, const Point pointTarget, double startTime) {
	numHaversine++;
	if ( (numHaversine % 7677439641) == 0 ) {
		double middleTime = timeCheckerCPU();
		double time = middleTime - startTime;
		printf( "Haversine    : %ld\n", numHaversine );
		printf( "Elapsed Time : %.8f\n", time );
		printf( "Complete Node: %ld\n", numComplete );
		fflush( stdout );
	}
	// Distance between latitudes and longitudes
	double dlat = (pointTarget.lat-pointCore.lat)*TO_RADIAN;
	double dlon = (pointTarget.lon-pointCore.lon)*TO_RADIAN;

	// Convert to radians
	double rad_lat_core = pointCore.lat*TO_RADIAN;
	double rad_lat_target = pointTarget.lat*TO_RADIAN;

	// Apply formula
	double f = pow(sin(dlat/2),2) + pow(sin(dlon/2),2) * cos(rad_lat_core) * cos(rad_lat_target);
	return asin(sqrt(f)) * 2 * EARTH_RADIUS;
}

// Inverse haversine for latitude
float inverseHaversineLat(const Point pointCore, size_t epsilon) {
	double dlat_1km = 0.008992;
	return dlat_1km * epsilon;
}

// Inverse haversine for longitude
float inverseHaversineLon(const Point pointCore, size_t epsilon) {
	double sinFunc = sin((epsilon * TO_RADIAN) / (2 * EARTH_RADIUS * TO_RADIAN));
	double powFunc = pow(sinFunc, 2);
	double secLat = 1 / cos(pointCore.lat * TO_RADIAN);
	return (2 * asin(sqrt(powFunc * secLat * secLat))) * TO_DEGREE;
}

// Function for four edge points of square
void findEdgePointsEpsilonBox(std::vector<PointDBSCAN> &dataset, size_t epsilon) {
	for ( size_t i = 0; i < dataset.size(); i ++ ) {
		double dlat = inverseHaversineLat(dataset[i].point, epsilon);
		double dlon = inverseHaversineLon(dataset[i].point, epsilon);
		dataset[i].northEastern.lat = dataset[i].point.lat + dlat;
		dataset[i].northEastern.lon = dataset[i].point.lon + dlon;
		dataset[i].northWestern.lat = dataset[i].point.lat - dlat;
		dataset[i].northWestern.lon = dataset[i].point.lon + dlon;
		dataset[i].southEastern.lat = dataset[i].point.lat + dlat;
		dataset[i].southEastern.lon = dataset[i].point.lon - dlon;
		dataset[i].southWestern.lat = dataset[i].point.lat - dlat;
		dataset[i].southWestern.lon = dataset[i].point.lon - dlon;
	}
}

// Function for four edge points of quadrant
void findEdgePointsQuadrant(Quadrant *root) {
	// First child quad
	root->child[0]->northEastern.lat = root->center.lat;
	root->child[0]->northEastern.lon = root->center.lon;
	root->child[0]->northWestern.lat = root->southWestern.lat;
	root->child[0]->northWestern.lon = root->center.lon;
	root->child[0]->southEastern.lat = root->center.lat;
	root->child[0]->southEastern.lon = root->southWestern.lon;
	root->child[0]->southWestern.lat = root->southWestern.lat;
	root->child[0]->southWestern.lon = root->southWestern.lon;

	// Second child quad
	root->child[1]->northEastern.lat = root->northEastern.lat;
	root->child[1]->northEastern.lon = root->center.lon;
	root->child[1]->northWestern.lat = root->center.lat;
	root->child[1]->northWestern.lon = root->center.lon;
	root->child[1]->southEastern.lat = root->northEastern.lat;
	root->child[1]->southEastern.lon = root->southWestern.lon;
	root->child[1]->southWestern.lat = root->center.lat;
	root->child[1]->southWestern.lon = root->southWestern.lon;

	// Third child quad
	root->child[2]->northEastern.lat = root->center.lat;
	root->child[2]->northEastern.lon = root->northEastern.lon;
	root->child[2]->northWestern.lat = root->southWestern.lat;
	root->child[2]->northWestern.lon = root->northEastern.lon;
	root->child[2]->southEastern.lat = root->center.lat;
	root->child[2]->southEastern.lon = root->center.lon;
	root->child[2]->southWestern.lat = root->southWestern.lat;
	root->child[2]->southWestern.lon = root->center.lon;

	// Fourth child quad
	root->child[3]->northEastern.lat = root->northEastern.lat;
	root->child[3]->northEastern.lon = root->northEastern.lon;
	root->child[3]->northWestern.lat = root->center.lat;
	root->child[3]->northWestern.lon = root->northEastern.lon;
	root->child[3]->southEastern.lat = root->northEastern.lat;
	root->child[3]->southEastern.lon = root->center.lon;
	root->child[3]->southWestern.lat = root->center.lat;
	root->child[3]->southWestern.lon = root->center.lon;
}

// Function for finding center mass value
void findCenterMass(Quadrant *root) {
	double totalX = 0.00;
	double totalY = 0.00;
	for ( size_t i = 0; i < root->cities.size(); i ++ ) {
		totalX = totalX + root->cities[i].point.lat;
		totalY = totalY + root->cities[i].point.lon;
	}
	root->center.lat = totalX / (double)root->cities.size();
	root->center.lon = totalY / (double)root->cities.size();
}

// Function for finding a length of diagonal haversine distance
void findDiagonal(Quadrant *root) {
	double useless = 0.00;
	root->diagonal = haversine(root->northEastern, root->southWestern, useless);
}

// Function for initialization
void initialize(Quadrant *root, size_t epsilon) {
	// Highest and lowest
	for ( size_t i = 0; i < root->cities.size(); i ++ ) {
		if ( i == 0 ) {
			root->southWestern.lat = root->cities[i].point.lat;
			root->southWestern.lon = root->cities[i].point.lon;
			root->northEastern.lat = root->cities[i].point.lat;
			root->northEastern.lon = root->cities[i].point.lon;
		} else {
			if ( root->southWestern.lat > root->cities[i].point.lat ) {
				root->southWestern.lat = root->cities[i].point.lat;
			}
			if ( root->southWestern.lon > root->cities[i].point.lon ) {
				root->southWestern.lon = root->cities[i].point.lon;
			}
			if ( root->northEastern.lat < root->cities[i].point.lat ) {
				root->northEastern.lat = root->cities[i].point.lat;
			}
			if ( root->northEastern.lon < root->cities[i].point.lon ) {
				root->northEastern.lon = root->cities[i].point.lon;
			}
		}
	}

	// Other edge points
	root->northWestern.lat = root->southWestern.lat;
	root->northWestern.lon = root->northEastern.lon;
	root->southEastern.lat = root->northEastern.lat;
	root->southEastern.lon = root->southWestern.lon;

	// Center mass of dataset
	findCenterMass(root);

	// Diagonal haversine distance
	findDiagonal(root);

	if ( root->diagonal <= epsilon ) root->done = 1;
}

// Quadtree (Divide parent quadrant to 4 childrent quadrant)
void divideQuad(Quadrant *root) {
	for ( size_t i = 0; i < root->cities.size(); i ++ ) {
		// First child quadrant
		if ( (root->cities[i].point.lat < root->center.lat) && 
		     (root->cities[i].point.lon < root->center.lon) ) {
			PointQuadTree temp;
			temp.point.lat = root->cities[i].point.lat;
			temp.point.lon = root->cities[i].point.lon;
			temp.datasetID = root->cities[i].datasetID;
			root->child[0]->cities.push_back(temp);
		// Second child quadrant
		} else if ( (root->cities[i].point.lat >= root->center.lat) && 
			    (root->cities[i].point.lon < root->center.lon) ) {
			PointQuadTree temp;
			temp.point.lat = root->cities[i].point.lat;
			temp.point.lon = root->cities[i].point.lon;
			temp.datasetID = root->cities[i].datasetID;
			root->child[1]->cities.push_back(temp);
		// Third child quadrant
		} else if ( (root->cities[i].point.lat < root->center.lat) && 
			    (root->cities[i].point.lon >= root->center.lon) ) {
			PointQuadTree temp;
			temp.point.lat = root->cities[i].point.lat;
			temp.point.lon = root->cities[i].point.lon;
			temp.datasetID = root->cities[i].datasetID;
			root->child[2]->cities.push_back(temp);
		// Fourth child quadrant
		} else if ( (root->cities[i].point.lat >= root->center.lat) && 
			    (root->cities[i].point.lon >= root->center.lon) ) {
			PointQuadTree temp;
			temp.point.lat = root->cities[i].point.lat;
			temp.point.lon = root->cities[i].point.lon;
			temp.datasetID = root->cities[i].datasetID;
			root->child[3]->cities.push_back(temp);
		}
	}
}

// Quadtree (Get the needed information for each quadrant)
void getInfoQuad(Quadrant *root, size_t epsilon, size_t pointsCondition, int level) {
	for ( size_t i = 0; i < root->child.size(); ) {
		if ( root->child[i]->cities.size() > (size_t)pointsCondition ) {
			findCenterMass(root->child[i]);
			findDiagonal(root->child[i]);
			if ( root->child[i]->diagonal <= epsilon ) {
				root->child[i]->done = 1;
				numDataPoints = numDataPoints + (uint64_t)root->child[i]->cities.size();
				uint64_t quotient = (uint64_t)root->child[i]->cities.size() / PCIECONDITION;
				uint64_t remainder = (uint64_t)root->child[i]->cities.size() % PCIECONDITION;
				numPointsCondition = numPointsCondition + (quotient * PCIECONDITION);
				numNormal = numNormal + remainder;
			} else {
				size_t same = 1;
				for ( size_t j = 1; j < root->child[i]->cities.size(); j ++ ) {
					if ( (root->child[i]->cities[0].point.lat == root->child[i]->cities[j].point.lat) &&
					     (root->child[i]->cities[0].point.lon == root->child[i]->cities[j].point.lon) ) same++;
				}
				if ( same == root->child[i]->cities.size() ) {
					root->child[i]->done = 1;
					root->child[i]->same = 1;
					numDataPoints = numDataPoints + (uint64_t)root->child[i]->cities.size();
					uint64_t quotient = (uint64_t)root->child[i]->cities.size() / PCIECONDITION;
					uint64_t remainder = (uint64_t)root->child[i]->cities.size() % PCIECONDITION;
					numPointsCondition = numPointsCondition + (quotient * PCIECONDITION);
					numNormal = numNormal + remainder;
				} else root->child[i]->done = 0;
			}
			i++;
		} else if ( root->child[i]->cities.size() == (size_t)pointsCondition ) {
			findCenterMass(root->child[i]);
			findDiagonal(root->child[i]);
			root->child[i]->done = 1;
			size_t same = 1;
			for ( size_t j = 1; j < root->child[i]->cities.size(); j ++ ) {
				if ( (root->child[i]->cities[0].point.lat == root->child[i]->cities[j].point.lat) &&
				     (root->child[i]->cities[0].point.lon == root->child[i]->cities[j].point.lon) ) same++;
			}
			if ( same == root->child[i]->cities.size() ) {
				root->child[i]->same = 1;
			}
			numDataPoints = numDataPoints + (uint64_t)root->child[i]->cities.size();
			uint64_t quotient = (uint64_t)root->child[i]->cities.size() / PCIECONDITION;
			uint64_t remainder = (uint64_t)root->child[i]->cities.size() % PCIECONDITION;
			numPointsCondition = numPointsCondition + (quotient * PCIECONDITION);
			numNormal = numNormal + remainder;
			i++;
		} else if ( root->child[i]->cities.size() < (size_t)pointsCondition ) {
			if ( root->child[i]->cities.size() == 0 ) {
				delete root->child[i];
				root->child.erase(root->child.begin() + i);
			} else {
				findCenterMass(root->child[i]);
				findDiagonal(root->child[i]);
				root->child[i]->done = 1;
				size_t same = 1;
				for ( size_t j = 1; j < root->child[i]->cities.size(); j ++ ) {
					if ( (root->child[i]->cities[0].point.lat == root->child[i]->cities[j].point.lat) &&
					     (root->child[i]->cities[0].point.lon == root->child[i]->cities[j].point.lon) ) same++;
				}
				if ( same == root->child[i]->cities.size() ) {
					root->child[i]->same = 1;
				}
				numDataPoints = numDataPoints + (uint64_t)root->child[i]->cities.size();
				uint64_t quotient = (uint64_t)root->child[i]->cities.size() / PCIECONDITION;
				uint64_t remainder = (uint64_t)root->child[i]->cities.size() % PCIECONDITION;
				numPointsCondition = numPointsCondition + (quotient * PCIECONDITION);
				numNormal = numNormal + remainder;
				i++;
			}
		} 	
	}
	
	// Count the total number of quadrants & delete parent's data
	if ( root->child.size() > 0 ) {
		numQuadrants = numQuadrants + (uint64_t)root->child.size();
		std::vector<PointQuadTree>().swap(root->cities);
	}
}

// Quadtree (Insert new child quadrant to parent quadrant)
void insertQuad(Quadrant *root, size_t epsilon, size_t pointsCondition, int level) {
	if ( root->done == 0 ) {
		// Generate child quadrants first
		root->child.resize(4);
		for ( size_t i = 0; i < 4; i ++ ) {
			root->child[i] = new Quadrant;
			root->child[i]->diagonal = 0.00; 
			root->child[i]->done = 0;
			root->child[i]->same = 0;
		}

		// Divide
		divideQuad(root);

		// Highest and lowest values of each quadrant
		findEdgePointsQuadrant(root);

		// Center mass value and diagonal distance of each quadrant
		getInfoQuad(root, epsilon, pointsCondition, level);
	} else return;
}

// Quadtree (Main)
int quadtree(Quadrant *root, size_t epsilon, size_t pointsCondition) {
	std::vector<Quadrant*> parentsCurr;
	std::vector<Quadrant*> parentsNext;
	parentsCurr.clear();
	parentsCurr.shrink_to_fit();
	parentsNext.clear();
	parentsNext.shrink_to_fit();
	int level = 0;

	// Root quadrant
	insertQuad(root, epsilon, pointsCondition, level);
	if ( root->child.size() > 0 ) {
		for ( size_t i = 0; i < root->child.size(); i ++ ) {
			if ( root->child[i]->done == 0 ) parentsCurr.push_back(root->child[i]);
		}
	}
	level++;
	printf( "DBSCAN Level: %d\n", level );
	fflush( stdout );

	// Iteration until meeting the terminate conditions
	while (1) {
		parentsNext.clear();
		parentsNext.shrink_to_fit();
		for ( size_t i = 0; i < parentsCurr.size(); i ++ ) {
			insertQuad(parentsCurr[i], epsilon, pointsCondition, level);
			if ( parentsCurr[i]->child.size() > 0 ) {
				for ( size_t j = 0; j < parentsCurr[i]->child.size(); j ++ ) {
					if ( parentsCurr[i]->child[j]->done == 0 ) parentsNext.push_back(parentsCurr[i]->child[j]);
				}
			}
		}
		parentsCurr.clear();
		parentsCurr.shrink_to_fit();
		if ( parentsNext.size() > 0 ) {
			for ( size_t i = 0; i < parentsNext.size(); i ++ ) {
				parentsCurr.push_back(parentsNext[i]);
			}
		} else break;
		level++;
		printf( "DBSCAN Level: %d\n", level );
		fflush( stdout );
	}
	
	return level;
}

// DBSCAN (Comparer between epsilon box and quadrant)
// Check the epsilon box is in a quadrant first
size_t compareEBinQ(std::vector<PointDBSCAN> &dataset, size_t index, Quadrant *root) {
	int numPoints = 0;
	if ( dataset[index].northEastern.lat >= root->southWestern.lat && 
	     dataset[index].northEastern.lat <= root->northEastern.lat &&
	     dataset[index].northEastern.lon >= root->southWestern.lon &&
	     dataset[index].northEastern.lon <= root->northEastern.lon ) numPoints++;
	if ( dataset[index].northWestern.lat >= root->southWestern.lat && 
	     dataset[index].northWestern.lat <= root->northEastern.lat &&
	     dataset[index].northWestern.lon >= root->southWestern.lon &&
	     dataset[index].northWestern.lon <= root->northEastern.lon ) numPoints++;
	if ( dataset[index].southEastern.lat >= root->southWestern.lat && 
	     dataset[index].southEastern.lat <= root->northEastern.lat &&
	     dataset[index].southEastern.lon >= root->southWestern.lon &&
	     dataset[index].southEastern.lon <= root->northEastern.lon ) numPoints++;
	if ( dataset[index].southWestern.lat >= root->southWestern.lat && 
	     dataset[index].southWestern.lat <= root->northEastern.lat &&
	     dataset[index].southWestern.lon >= root->southWestern.lon &&
	     dataset[index].southWestern.lon <= root->northEastern.lon ) numPoints++;
	return numPoints;
}
// Check the quadrant is in epsilon box
size_t compareQinEB(std::vector<PointDBSCAN> &dataset, size_t index, Quadrant *root) {
	int numPoints = 0;
	if ( root->northEastern.lat >= dataset[index].southWestern.lat && 
	     root->northEastern.lat <= dataset[index].northEastern.lat &&
	     root->northEastern.lon >= dataset[index].southWestern.lon &&
	     root->northEastern.lon <= dataset[index].northEastern.lon ) numPoints++;
	if ( root->northWestern.lat >= dataset[index].southWestern.lat && 
	     root->northWestern.lat <= dataset[index].northEastern.lat &&
	     root->northWestern.lon >= dataset[index].southWestern.lon &&
	     root->northWestern.lon <= dataset[index].northEastern.lon ) numPoints++;
	if ( root->southEastern.lat >= dataset[index].southWestern.lat && 
	     root->southEastern.lat <= dataset[index].northEastern.lat &&
	     root->southEastern.lon >= dataset[index].southWestern.lon &&
	     root->southEastern.lon <= dataset[index].northEastern.lon ) numPoints++;
	if ( root->southWestern.lat >= dataset[index].southWestern.lat && 
	     root->southWestern.lat <= dataset[index].northEastern.lat &&
	     root->southWestern.lon >= dataset[index].southWestern.lon &&
	     root->southWestern.lon <= dataset[index].northEastern.lon ) numPoints++;
	return numPoints;
}
// Check epsilon box and quadrant are overlapped each other
size_t compareOverlap(std::vector<PointDBSCAN> &dataset, size_t index, Quadrant *root) {
	int numPoints = 0;
	if ( dataset[index].northEastern.lat >= root->southWestern.lat &&
	     dataset[index].northEastern.lat <= root->northEastern.lat &&
	     root->northEastern.lon >= dataset[index].southWestern.lon &&
	     root->northEastern.lon <= dataset[index].northEastern.lon) numPoints++;
	if ( root->northEastern.lat >= dataset[index].southWestern.lat &&
	     root->northEastern.lat <= dataset[index].northEastern.lat &&
	     dataset[index].northEastern.lon >= root->southWestern.lon &&
	     dataset[index].northEastern.lon <= root->northEastern.lon) numPoints++;
	return numPoints;
}

// DBSCAN (Do haversine calculation based on candidate list)
void candidateListCalculator(std::vector<PointDBSCAN> &dataset, size_t index, std::vector<size_t> &borders, Quadrant *root, size_t epsilon, double startTime) {
	if ( root->done == 0 ) {
		for ( size_t i = 0; i < root->child.size(); i ++ ) {
			candidateListCalculator(dataset, index, borders, root->child[i], epsilon, startTime);
		}
	} else {
		if ( root->diagonal <= epsilon ) {
			for ( size_t i = 0; i < root->cities.size(); i ++ ) {
				if ( haversine(dataset[index].point, root->cities[i].point, startTime) <= epsilon ) {
					for ( size_t j = 0; j < root->cities.size(); j ++ ) {
						borders.push_back(root->cities[j].datasetID);
					}
					break;
				}
			}
		} else {
			if ( root->same == 1 ) {
				if ( haversine(dataset[index].point, root->cities[0].point, startTime) <= epsilon ) {
					for ( size_t j = 0; j < root->cities.size(); j ++ ) {
						borders.push_back(root->cities[j].datasetID);
					}
				}
			} else {
				for ( size_t i = 0; i < root->cities.size(); i ++ ) {
					if ( haversine(dataset[index].point, root->cities[i].point, startTime) <= epsilon ) {
						borders.push_back(root->cities[i].datasetID);
					}
				}
			}
		}
	}
}

// DBSCAN (Do make candidate list in case of quadrant is in epsilon box)
void findQuadrantsQinEB(std::vector<PointDBSCAN> &dataset, size_t index, std::vector<size_t> &borders, Quadrant *root, size_t epsilon, double startTime) {
	size_t resultQinEB = compareQinEB(dataset, index, root);
	if ( resultQinEB == 1 || resultQinEB == 2 || resultQinEB == 3 ) {
		if ( root->done == 0 ) {
			for ( size_t i = 0; i < root->child.size(); i ++ ) {
				findQuadrantsQinEB(dataset, index, borders, root->child[i], epsilon, startTime);
			}
		} else {
			candidateListCalculator(dataset, index, borders, root, epsilon, startTime);
		}
	} else if ( resultQinEB == 4 ) {
		candidateListCalculator(dataset, index, borders, root, epsilon, startTime);
	}
}

// DBSCAN (Do make candidate list in case of epsilon box is in quadrant)
void findQuadrantsEBinQ(std::vector<PointDBSCAN> &dataset, size_t index, std::vector<size_t> &borders, Quadrant *root, size_t epsilon, double startTime) {
	size_t resultEBinQ = compareEBinQ(dataset, index, root);
	if ( resultEBinQ != 0 ) {
		if ( root->done == 1 ) {
			candidateListCalculator(dataset, index, borders, root, epsilon, startTime);
		} else {
			size_t resultQinEB = compareQinEB(dataset, index, root);
			if ( resultQinEB == 0 ) {
				for ( size_t i = 0; i < root->child.size(); i ++ ) {
					findQuadrantsEBinQ(dataset, index, borders, root->child[i], epsilon, startTime);
				}
			} else findQuadrantsQinEB(dataset, index, borders, root, epsilon, startTime);
		}
	} else {
		size_t resultQinEB = compareQinEB(dataset, index, root);
		if ( resultQinEB != 0 ) {
			findQuadrantsQinEB(dataset, index, borders, root, epsilon, startTime);
		} else {
			size_t resultPart3 = compareOverlap(dataset, index, root);
			if ( resultPart3 != 0 ) {
				if ( root->done == 1 ) candidateListCalculator(dataset, index, borders, root, epsilon, startTime);
				else {
					for ( size_t i = 0; i < root->child.size(); i ++ ) {
						findQuadrantsEBinQ(dataset, index, borders, root->child[i], epsilon, startTime);
					}
				}
			}
		}
	}
}

// DBSCAN (Border Point Finder of Core Point)
void borderFinderCore(std::vector<PointDBSCAN> &dataset, size_t corePoint, std::vector<size_t> &bordersCore, Quadrant *root, size_t epsilon, double startTime) {
	bordersCore.clear();
	bordersCore.shrink_to_fit();
	findQuadrantsEBinQ(dataset, corePoint, bordersCore, root, epsilon, startTime);
}

// DBSCAN (Border Point Finder of Border Point)
void borderFinderBorder(std::vector<PointDBSCAN> &dataset, size_t borderPoint, std::vector<size_t> &bordersBorder, Quadrant *root, size_t epsilon, double startTime) {
	bordersBorder.clear();
	bordersBorder.shrink_to_fit();
	findQuadrantsEBinQ(dataset, borderPoint, bordersBorder, root, epsilon, startTime);
}

// DBSCAN (Cluster Expander)
int clusterExpander(std::vector<PointDBSCAN> &dataset, size_t index, size_t clusterID, Quadrant *root, size_t epsilon, double startTime) {
	std::vector<size_t> bordersCore;
	std::vector<size_t> bordersBorder;
	borderFinderCore(dataset, index, bordersCore, root, epsilon, startTime);

	if ( bordersCore.size() < MINIMUM_POINTS ) {
		dataset[index].clusterID = NOISE;
		std::vector<size_t>().swap(bordersCore);
		return FAILURE;
	} else {
		for ( size_t i = 0; i < bordersCore.size(); i ++ ) {
			size_t borderPoint = bordersCore[i];
			dataset[borderPoint].clusterID = clusterID;
			numComplete++;
		}

		for ( size_t i = 0; i < bordersCore.size(); i ++ ) {
			size_t borderPoint = bordersCore[i];
			if ( (dataset[borderPoint].point.lat == dataset[index].point.lat) && 
			     (dataset[borderPoint].point.lon == dataset[index].point.lon) ) {
				continue;
			} else {
				borderFinderBorder(dataset, borderPoint, bordersBorder, root, epsilon, startTime);

				if ( bordersBorder.size() >= MINIMUM_POINTS ) {
					for ( size_t j = 0; j < bordersBorder.size(); j ++ ) {
						size_t neighbourPoint = bordersBorder[j];
						if ( dataset[neighbourPoint].clusterID == (size_t)UNCLASSIFIED ||
						     dataset[neighbourPoint].clusterID == (size_t)NOISE ) {
							if ( dataset[neighbourPoint].clusterID == (size_t)UNCLASSIFIED ) {
								bordersCore.push_back(neighbourPoint);
							}
							dataset[neighbourPoint].clusterID = clusterID;
							numComplete++;
						}
					}
				}
			}
		}
		std::vector<size_t>().swap(bordersCore);
		std::vector<size_t>().swap(bordersBorder);
		return SUCCESS;
	}
}

// DBSCAN (Main)
int dbscan(std::vector<PointDBSCAN> &dataset, Quadrant *root, size_t epsilon, double startTime) {
	size_t clusterID = 1;
	for ( size_t i = 0; i < dataset.size(); i ++ ) {
		if ( dataset[i].clusterID == (size_t)UNCLASSIFIED ) {
			if ( clusterExpander(dataset, i, clusterID, root, epsilon, startTime) != FAILURE ) {
				clusterID += 1;
				printf( "Generating cluster %lu done!\n", clusterID-1 );
				fflush( stdout );
			}
		}
	}

	return clusterID-1;
}

// Function for printing results
void printResults(std::vector<PointDBSCAN> &dataset) {
	printf(" x       y       cluster_id\n"
	       "---------------------------\n");

	for ( int i = 0; i < (int)dataset.size() - 1; i ++ ) {
		for ( int j = 0; j < (int)dataset.size() - 1 - i; j ++ ) {
			if ( dataset[j].clusterID > dataset[j+1].clusterID ) {
				PointDBSCAN temp = dataset[j];
				dataset[j] = dataset[j+1];
				dataset[j+1] = temp;
			}
		}
	}

	for ( int i = 0; i < (int)dataset.size(); i ++ ) {
		printf("%f, %f: %lu\n", dataset[i].point.lat, dataset[i].point.lon, dataset[i].clusterID);
	}

	printf( "--------------------------\n" );
}

// Main
int main(int argc, char **argv) {
	size_t numCities = 44691;//700968*160;
	//size_t numCities = 2899550649;
	size_t epsilon = atoi(argv[1]);
	size_t pointsCondition = atoi(argv[2]);
	printf( "Epsilon Value                              : %lu\n", epsilon );
	fflush( stdout );
	printf( "Elements in a Leaf Node                    : %lu\n", pointsCondition );
	fflush( stdout );

	Quadrant *root = new Quadrant;

	// Read point data
	char benchmark_filename[] = "../../../worldcities.bin";
	readBenchmarkDataDBSCAN(dataset, benchmark_filename, numCities);
	printf( "Read File Done![First]\n" );
	fflush( stdout );
	readBenchmarkDataQuadTree(root, benchmark_filename, numCities);
	printf( "Read File Done![Second]\n" );
	fflush( stdout );

	// Initialize
	initialize(root, epsilon);

	// Get four edge points of epsilon box of each data point
	printf( "Finding Four Edge Points of Epsilon Box of The World Cities Start!\n" );
	fflush( stdout );
	double processStartStep1 = timeCheckerCPU();
	findEdgePointsEpsilonBox(dataset, epsilon);
	double processFinishStep1 = timeCheckerCPU();
	double processTimeStep1 = processFinishStep1 - processStartStep1;
	printf( "Finding Four Edge Points of Epsilon Box of The World Cities Done!\n" );
	printf( "\n" );
	fflush( stdout );
	printf( "Elapsed Time [Step1] [Epsilon Box] (CPU)   : %.8f\n", processTimeStep1 );
	fflush( stdout );

	// Quadtree
	printf( "Quadtree for The World Cities Start!\n" );
	fflush( stdout );
	double processStartStep2 = timeCheckerCPU();
	int level = quadtree(root, epsilon, pointsCondition);
	double processFinishStep2 = timeCheckerCPU();
	double processTimeStep2 = processFinishStep2 - processStartStep2;
	printf( "Quadtree for The World Cities Done!\n" );
	printf( "\n" );
	fflush( stdout );
	float notCompressSize = (float)numNormal * 64.00;
	float compressSize16b = (float)numPointsCondition * 32.00;
	float compressSize12b = (float)numPointsCondition * 24.00;
	float compressSize08b = (float)numPointsCondition * 16.00;
	float dataSize16b = (((((notCompressSize + compressSize16b) / 8) / 1024) / 1024) / 1024);
	float dataSize12b = (((((notCompressSize + compressSize12b) / 8) / 1024) / 1024) / 1024);
	float dataSize08b = (((((notCompressSize + compressSize08b) / 8) / 1024) / 1024) / 1024);
	printf( "Elapsed Time [Step2] [Quadtree] (CPU)      : %.8f\n", processTimeStep2 );
	printf( "The Number of Data Points [Total]          : %ld\n", numDataPoints );
	printf( "The Number of Data Points [Compress]       : %ld\n", numPointsCondition );
	printf( "The Number of Data Points [Not Compress]   : %ld\n", numNormal );
	printf( "Data Structure Size [16 bits]              : %.8f\n", dataSize16b );
	printf( "Data Structure Size [12 bits]              : %.8f\n", dataSize12b );
	printf( "Data Structure Size [08 bits]              : %.8f\n", dataSize08b );
	printf( "The Maximum of Tree Level                  : %d\n", level );
	printf( "The Number of Quadrants                    : %ld\n", numQuadrants );


	// DBSCAN
	printf( "Quadtree-based DBSCAN Clustering for The World Cities Start!\n" );
	double processStartStep3 = timeCheckerCPU();
	int maxClusterID = dbscan(dataset, root, epsilon, processStartStep3);
	double processFinishStep3 = timeCheckerCPU();
	double processTimeStep3 = processFinishStep3 - processStartStep3;
	printf( "Quadtree-based DBSCAN Clustering for The World Cities Done!\n" );
	printf( "\n" );
	fflush( stdout );
/*
	// Save the result dataset as bin file
	char resultDataset_filename[] = "../worldcities_result.bin";
	printf( "Saving The Result Dataset Start!\n" );
	writeBenchmarkDataDBSCAN(dataset, resultDataset_filename, numCities);
	printf( "Saving The Result Dataset Done!\n" );
	printf( "\n" );
	fflush( stdout );
*/	
	// Result of Quadtree-based DBSCAN algorithm
	//printResults(dataset);
	double haversinePercent_1 = (double)numHaversine / (double)7677439641;
	double haversinePercent_2 = haversinePercent_1 * (double)216.46002400;
	double haversinePercent_3 = haversinePercent_2 / (double)processTimeStep3;
	double haversinePercent_4 = haversinePercent_3 * (double)100.00;
	double haversinePercent_5 = haversinePercent_4 + (double)31.00;

	printf( "Elapsed Time [Step3] [DBSCAN] (CPU)        : %.8f\n", processTimeStep3 );
	printf( "The Number of Haversine                    : %ld\n", numHaversine );
	printf( "Haversine Percentage of Whole Elased Time  : %.8f\n", haversinePercent_5 );
	printf( "Max Cluster ID                             : %d\n", maxClusterID );

	delete root;	
	return 0;
}
