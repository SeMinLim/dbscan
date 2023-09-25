#include <sys/resource.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>


// DBSCAN
#define CORE_POINT 1
#define BORDER_POINT 2

#define SUCCESS 0
#define UNCLASSIFIED -1
#define NOISE -2
#define FAILURE -3

#define MINIMUM_POINTS 2
#define EPSILON 0.1

#define NUMPOINTS (87621*2) // (700968 * 160), 112154880


typedef struct Point {
	float lat, lon; // latitude & longitude
	int clusterID;  // clustered ID
}Point;


Point m_points[NUMPOINTS+1];
int clusterSeeds[NUMPOINTS+1];
int clusterNeighbours[NUMPOINTS+1];

int clusterSeedsSize = 0;
int clusterNeighboursSize = 0;
uint64_t numManhattan = 0;


// Elapsed time checker
static inline double timeCheckerCPU(void) {
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

// Function for reading benchmark file
void readBenchmarkData(char* filename) {
	FILE *data = fopen(filename, "rb");
	if( data == NULL ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	for ( int i = 0; i < NUMPOINTS; i ++ ) {
		fread(&m_points[i].lat, sizeof(float), 1, data);
		fread(&m_points[i].lon, sizeof(float), 1, data);
		m_points[i].clusterID = UNCLASSIFIED;
	}

	fclose(data);
}

// Manhattan
float manhattan(const Point pointCore, const Point pointTarget) {
	numManhattan++;
	float sub_lat = pointCore.lat - pointTarget.lat;
	float sub_lon = pointCore.lon - pointTarget.lon;
	
	return abs(sub_lat) + abs(sub_lon);
}

// Border point finder for core point
void calculateClusterSeeds(Point point) {
	clusterSeedsSize = 0;
	
	for( int i = 0; i < NUMPOINTS; i ++ ) {
		if ( manhattan(point, m_points[i]) <= EPSILON ) {
			clusterSeeds[clusterSeedsSize++] = i;
		}
	}
}

// Border point finder for border point
void calculateClusterNeighbours(Point point) {
	clusterNeighboursSize = 0;

	for( int i = 0; i < NUMPOINTS; i ++ ) {
		if ( manhattan(point, m_points[i]) <= EPSILON ) {
			clusterNeighbours[clusterNeighboursSize++] = i;
		}
	}
}

// Cluster expander
int expandCluster(Point point, int clusterID) {    
	calculateClusterSeeds(point);
	
	if ( clusterSeedsSize < MINIMUM_POINTS ) {
		point.clusterID = NOISE;
		return FAILURE;
	} else {
		for( int i = 0; i < clusterSeedsSize; i ++ ) {
			int pointBorder = clusterSeeds[i];
			m_points[pointBorder].clusterID = clusterID;
		}
			
		for( int i = 0; i < clusterSeedsSize; i ++ ) {
			int pointBorder = clusterSeeds[i];
			if ( m_points[pointBorder].lat == point.lat &&
			     m_points[pointBorder].lon == point.lon ) {
				continue;
			} else {
				calculateClusterNeighbours(m_points[pointBorder]);
			
				if ( clusterNeighboursSize >= MINIMUM_POINTS ) {
					for ( int j = 0; j < clusterNeighboursSize; j ++ ) {
						int pointNeighbour = clusterNeighbours[j];
						if ( m_points[pointNeighbour].clusterID == UNCLASSIFIED || 
						     m_points[pointNeighbour].clusterID == NOISE ) {
							if ( m_points[pointNeighbour].clusterID == UNCLASSIFIED ) {
								clusterSeeds[clusterSeedsSize++] = pointNeighbour;
							}
							m_points[pointNeighbour].clusterID = clusterID;
						}
					}
				}
			}
		}
		return SUCCESS;
	}
}

// DBSCAN main
int run() {
	int clusterID = 1;
	for( int i = 0; i < NUMPOINTS; i ++ ) {
		if ( m_points[i].clusterID == UNCLASSIFIED ) {
			if ( expandCluster(m_points[i], clusterID) != FAILURE ) {
				clusterID += 1;
				printf( "Generating cluster %d done!\n", clusterID-1 );
			}
		}
	}

	return clusterID-1;
}

// Function for printing results
void printResults() {
	int i = 0;

	printf(" x       y       cluster_id\n"
	       "---------------------------\n");

	while (i < NUMPOINTS) {
		printf("%f, %f: %d\n", m_points[i].lat, m_points[i].lon, m_points[i].clusterID);
		++i;
	}

	printf( "--------------------------\n" );
}

int main() {
	char benchmark_filename[] = "../../worldcities_augmented.bin";

	// read point data
	readBenchmarkData(benchmark_filename);

	// main loop
	printf( "Pure DBSCAN Clustering for The World Cities Start!\n" );
	double processStart = timeCheckerCPU();
	int maxClusterID = run();
	double processFinish = timeCheckerCPU();
	double processTime = processFinish - processStart;
	printf( "Pure DBSCAN Clustering for The World Cities Done!\n" );
	printf( "\n" );

	// result of DBSCAN algorithm
	//printResults();    
	printf( "Elapsed Time [DBSCAN] (CPU): %.8f\n", processTime );
	printf( "The Number of Data Points  : %d\n", NUMPOINTS );
	printf( "The Number of Manhattan    : %ld\n", numManhattan );
	printf( "Max Cluster ID             : %d\n", maxClusterID );

	return 0;
}
