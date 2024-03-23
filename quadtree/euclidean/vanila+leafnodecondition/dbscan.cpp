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


uint64_t numDataPoints = 0;
uint64_t numQuadrants = 0;
uint64_t numEuclidean = 0;
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

// Euclidean
float euclidean(const Point pointCore, const Point pointTarget) {
	numEuclidean++;
	float sub_lat = pointCore.lat - pointTarget.lat;
	float sub_lon = pointCore.lon - pointTarget.lon;

	float before_sqrt = pow(sub_lat, 2) + pow(sub_lon, 2);

	return sqrt(before_sqrt);
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

// Function for finding a length of diagonal euclidean distance
void findDiagonal(Quadrant *root) {
	root->diagonal = euclidean(root->northEastern, root->southWestern);
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

	// Diagonal euclidean distance
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
size_t compareEBinQ(std::vector<PointDBSCAN> &dataset, size_t index, Quadrant *root, double epsilon) {
	int numPoints = 0;
	// northEastern
	if ( (dataset[index].point.lat + epsilon) >= root->southWestern.lat && 
	     (dataset[index].point.lat + epsilon) <= root->northEastern.lat &&
	     (dataset[index].point.lon + epsilon) >= root->southWestern.lon &&
	     (dataset[index].point.lon + epsilon) <= root->northEastern.lon ) numPoints++;
	// northWestern
	if ( (dataset[index].point.lat - epsilon) >= root->southWestern.lat && 
	     (dataset[index].point.lat - epsilon) <= root->northEastern.lat &&
	     (dataset[index].point.lon + epsilon) >= root->southWestern.lon &&
	     (dataset[index].point.lon + epsilon) <= root->northEastern.lon ) numPoints++;
	// southEastern
	if ( (dataset[index].point.lat + epsilon) >= root->southWestern.lat && 
	     (dataset[index].point.lat + epsilon) <= root->northEastern.lat &&
	     (dataset[index].point.lon - epsilon) >= root->southWestern.lon &&
	     (dataset[index].point.lon - epsilon) <= root->northEastern.lon ) numPoints++;
	// southWestern
	if ( (dataset[index].point.lat - epsilon) >= root->southWestern.lat && 
	     (dataset[index].point.lat - epsilon) <= root->northEastern.lat &&
	     (dataset[index].point.lon - epsilon) >= root->southWestern.lon &&
	     (dataset[index].point.lon - epsilon) <= root->northEastern.lon ) numPoints++;
	return numPoints;
}
// Check the quadrant is in epsilon box
size_t compareQinEB(std::vector<PointDBSCAN> &dataset, size_t index, Quadrant *root, double epsilon) {
	int numPoints = 0;
	if ( root->northEastern.lat >= (dataset[index].point.lat - epsilon) && 
	     root->northEastern.lat <= (dataset[index].point.lat + epsilon) &&
	     root->northEastern.lon >= (dataset[index].point.lon - epsilon) &&
	     root->northEastern.lon <= (dataset[index].point.lon + epsilon) ) numPoints++;
	if ( root->northWestern.lat >= (dataset[index].point.lat - epsilon) && 
	     root->northWestern.lat <= (dataset[index].point.lat + epsilon) &&
	     root->northWestern.lon >= (dataset[index].point.lon - epsilon) &&
	     root->northWestern.lon <= (dataset[index].point.lon + epsilon) ) numPoints++;
	if ( root->southEastern.lat >= (dataset[index].point.lat - epsilon) && 
	     root->southEastern.lat <= (dataset[index].point.lat + epsilon) &&
	     root->southEastern.lon >= (dataset[index].point.lon - epsilon) &&
	     root->southEastern.lon <= (dataset[index].point.lon + epsilon) ) numPoints++;
	if ( root->southWestern.lat >= (dataset[index].point.lat - epsilon) && 
	     root->southWestern.lat <= (dataset[index].point.lat + epsilon) &&
	     root->southWestern.lon >= (dataset[index].point.lon - epsilon) &&
	     root->southWestern.lon <= (dataset[index].point.lon + epsilon) ) numPoints++;
	return numPoints;
}
// Check epsilon box and quadrant are overlapped each other
size_t compareOverlap(std::vector<PointDBSCAN> &dataset, size_t index, Quadrant *root, double epsilon) {
	int numPoints = 0;
	if ( (dataset[index].point.lat + epsilon) >= root->southWestern.lat &&
	     (dataset[index].point.lat + epsilon) <= root->northEastern.lat &&
	     root->northEastern.lon >= (dataset[index].point.lon - epsilon) &&
	     root->northEastern.lon <= (dataset[index].point.lon + epsilon) ) numPoints++;
	if ( root->northEastern.lat >= (dataset[index].point.lat - epsilon) &&
	     root->northEastern.lat <= (dataset[index].point.lat + epsilon) &&
	     (dataset[index].point.lon + epsilon) >= root->southWestern.lon &&
	     (dataset[index].point.lon + epsilon) <= root->northEastern.lon) numPoints++;
	return numPoints;
}

// DBSCAN (Do euclidean calculation based on candidate list)
void candidateListCalculator(std::vector<PointDBSCAN> &dataset, size_t index, std::vector<size_t> &borders, Quadrant *root, size_t epsilon, double startTime) {
	if ( root->done == 0 ) {
		for ( size_t i = 0; i < root->child.size(); i ++ ) {
			candidateListCalculator(dataset, index, borders, root->child[i], epsilon, startTime);
		}
	} else {
		if ( root->diagonal <= epsilon ) {
			for ( size_t i = 0; i < root->cities.size(); i ++ ) {
				if ( euclidean(dataset[index].point, root->cities[i].point) <= epsilon ) {
					for ( size_t j = 0; j < root->cities.size(); j ++ ) {
						borders.push_back(root->cities[j].datasetID);
					}
					break;
				}
			}
		} else {
			if ( root->same == 1 ) {
				if ( euclidean(dataset[index].point, root->cities[0].point) <= epsilon ) {
					for ( size_t j = 0; j < root->cities.size(); j ++ ) {
						borders.push_back(root->cities[j].datasetID);
					}
				}
			} else {
				for ( size_t i = 0; i < root->cities.size(); i ++ ) {
					if ( euclidean(dataset[index].point, root->cities[i].point) <= epsilon ) {
						borders.push_back(root->cities[i].datasetID);
					}
				}
			}
		}
	}
}

// DBSCAN (Do make candidate list in case of quadrant is in epsilon box)
void findQuadrantsQinEB(std::vector<PointDBSCAN> &dataset, size_t index, std::vector<size_t> &borders, Quadrant *root, size_t epsilon, double startTime) {
	size_t resultQinEB = compareQinEB(dataset, index, root, epsilon);
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
	size_t resultEBinQ = compareEBinQ(dataset, index, root, epsilon);
	if ( resultEBinQ != 0 ) {
		if ( root->done == 1 ) {
			candidateListCalculator(dataset, index, borders, root, epsilon, startTime);
		} else {
			size_t resultQinEB = compareQinEB(dataset, index, root, epsilon);
			if ( resultQinEB == 0 ) {
				for ( size_t i = 0; i < root->child.size(); i ++ ) {
					findQuadrantsEBinQ(dataset, index, borders, root->child[i], epsilon, startTime);
				}
			} else findQuadrantsQinEB(dataset, index, borders, root, epsilon, startTime);
		}
	} else {
		size_t resultQinEB = compareQinEB(dataset, index, root, epsilon);
		if ( resultQinEB != 0 ) {
			findQuadrantsQinEB(dataset, index, borders, root, epsilon, startTime);
		} else {
			size_t resultPart3 = compareOverlap(dataset, index, root, epsilon);
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
			/*
			numComplete++;
			if ( numComplete % 10000 ) {
				double middleTime = timeCheckerCPU();
				double time = middleTime - startTime;
				printf( "Euclidean    : %ld\n", numEuclidean );
				printf( "Complete     : %ld\n", numComplete );
				printf( "Elapsed Time : %.8f\n", time );
				printf( "\n" );
				fflush( stdout );
			}
			*/
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
							/*
							numComplete++;
							if ( numComplete % 10000 ) {
								double middleTime = timeCheckerCPU();
								double time = middleTime - startTime;
								printf( "Euclidean    : %ld\n", numEuclidean );
								printf( "Complete     : %ld\n", numComplete );
								printf( "Elapsed Time : %.8f\n", time );
								printf( "\n" );
								fflush( stdout );
							}
							*/
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
				//printf( "Generating cluster %lu done!\n", clusterID-1 );
				//fflush( stdout );
			}
		}
	}

	return clusterID-1;
}

// Main
int main(int argc, char **argv) {
	size_t numCities = 44691; //(2899550649);
	double epsilon = 0.7;
	size_t pointsCondition = atoi(argv[1]);
	printf( "Epsilon Value                              : %f\n", epsilon );
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
	fflush( stdout );


	// DBSCAN
	printf( "Quadtree-based DBSCAN Clustering for The World Cities Start!\n" );
	fflush( stdout );
	double processStartStep3 = timeCheckerCPU();
	int maxClusterID = dbscan(dataset, root, epsilon, processStartStep3);
	double processFinishStep3 = timeCheckerCPU();
	double processTimeStep3 = processFinishStep3 - processStartStep3;
	printf( "Quadtree-based DBSCAN Clustering for The World Cities Done!\n" );
	printf( "\n" );
	fflush( stdout );
	
	// Result of Quadtree-based DBSCAN algorithm
	//printResults(dataset);
	printf( "Elapsed Time [Step3] [DBSCAN] (CPU)        : %.8f\n", processTimeStep3 );
	printf( "The Number of Euclidean                    : %ld\n", numEuclidean );
	printf( "Max Cluster ID                             : %d\n", maxClusterID );

	delete root;	
	return 0;
}
