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
#define EPSILON 0.7

// ZFP
#define BIT_BUDGET 13
// Exponent of single is 8 bit signed integer (-126 to +127)
#define EXP_MAX ((1<<(8-1))-1)
// Exponent of the minimum granularity value expressible via single
#define ZFP_MIN_EXP -149
#define ZFP_MAX_PREC 32
#define EBITS (8+1)


uint64_t numDataPoints = 0;
uint64_t numQuadrants = 0;
uint64_t numManhattan = 0;


typedef struct Point {
	float lat, lon;
}Point;

typedef struct PointQuadTree {
	Point point;
	int datasetID;
}PointQuadTree;

typedef struct PointDBSCAN {
	Point point;
	Point northEastern;
	Point northWestern;
	Point southEastern;
	Point southWestern;
	int clusterID_org;
	int clusterID_new;
}PointDBSCAN;

typedef struct Quadrant {
	std::vector<Quadrant*> child;
	std::vector<PointQuadTree> cities;
	Point northEastern;
	Point northWestern;
	Point southEastern;
	Point southWestern;
	Point center;
	float diagonal;
	bool done;
}Quadrant;

class BitBuffer {
public:
	BitBuffer(int bytes);
	~BitBuffer();
	void EncodeBit(uint32_t src);
	void EncodeBits(uint32_t src, int bits);
	int BitCount() { return curbits; };

	void DecodeBit(uint32_t* dst);
	void DecodeBits(uint32_t* dst, int bits);

	uint8_t* buffer;
	int bufferbytes;
	int curbits;
	int decoff;
};
BitBuffer::BitBuffer(int bytes) {
	this->buffer = (uint8_t*)malloc(bytes);
	for ( int i = 0; i < bytes; i++ ) this->buffer[i] = 0;

	this->bufferbytes = bytes;
	this->curbits = 0;
	this->decoff = 0;
}
BitBuffer::~BitBuffer() {
	free(this->buffer);
}
void BitBuffer::EncodeBit(uint32_t src) {
	int byteoff = (curbits/8);
	int bitoff = curbits%8;
	buffer[byteoff] |= ((src&1)<<bitoff);

	curbits++;
}
void BitBuffer::EncodeBits(uint32_t src, int bits) {
	for ( int i = 0; i < bits; i++ ) {
		EncodeBit(src);
		src >>= 1;
	}
}	
void BitBuffer::DecodeBit(uint32_t* dst) {
	int byteoff = (decoff/8);
	int bitoff = decoff%8;
	uint8_t buf = buffer[byteoff];
	*dst = (buf>>bitoff) & 1;
	decoff++;
}
void BitBuffer::DecodeBits(uint32_t* dst, int bits) {
	*dst = 0;
	for ( int i = 0; i < bits; i++ ) {
		uint32_t o = 0;
		DecodeBit(&o);
		*dst |= (o<<i);
	}
}


// Convert negabinary uint to int
static int32_t uint2int_uint32(uint32_t x) {
	return (int32_t)((x ^ 0xaaaaaaaa) - 0xaaaaaaaa);
}

// Convert int to negabinary uint
static uint32_t int2uint_int32(int32_t x) {
	return ((uint32_t)x + 0xaaaaaaaa) ^ 0xaaaaaaaa;
}

// Forward lifting
static void fwd_lift_int32(int32_t* p, uint s) {
	int32_t x, y, z, w;
	x = *p; p += s;
	y = *p; p += s;
	z = *p; p += s;
	w = *p; p += s;
	
	x += w; x >>= 1; w -= x;
	z += y; z >>= 1; y -= z;
	x += z; x >>= 1; z -= x;
	w += y; w >>= 1; y -= w;
	w += y >> 1; y -= w >> 1;
	
	p -= s; *p = w;
	p -= s; *p = z;
	p -= s; *p = y;
	p -= s; *p = x;
}

// Inverse lifting
static void inv_lift_int32(int32_t* p, uint32_t s) {
	int32_t x, y, z, w;
	x = *p; p += s;
	y = *p; p += s;
	z = *p; p += s;
	w = *p; p += s;
	
	y += w >> 1; w -= y >> 1;
	y += w; w <<= 1; w -= y;
	z += x; x <<= 1; x -= z;
	y += z; z <<= 1; z -= y;
	w += x; x <<= 1; x -= w;
	
	p -= s; *p = w;
	p -= s; *p = z;
	p -= s; *p = y;
	p -= s; *p = x;
}

// Compressor
void compress_1d(float origin[4], BitBuffer* decompressed, int bit_budget) {
	int exp_max = -EXP_MAX;
	for ( int i = 0; i < 4; i++ ) {
		if ( origin[i] != 0 ) {
			int exp = 0;
			std::frexp(origin[i], &exp);
			if ( exp > exp_max ) exp_max = exp;
		}
	}
	int dimension = 1;
	int precision_max = std::min(ZFP_MAX_PREC, std::max(0, exp_max - ZFP_MIN_EXP + (2*(dimension+1))));
	if ( precision_max != 0 ) {
		int e = exp_max + EXP_MAX;
		int32_t idata[4];
		for ( int i = 0; i < 4; i++ ) {
			idata[i] = (int32_t)(origin[i]*(pow(2, 32-2 - exp_max)));
		}
		
		// perform lifting
		fwd_lift_int32(idata, 1);
		
		// convert to negabinary
		uint32_t udata[4];
		for ( int i = 0; i < 4; i++ ) { 
			udata[i] = int2uint_int32(idata[i]);
		}

		int total_bits = EBITS;
		decompressed->EncodeBits(e, EBITS);

		for ( int i = 0; i < 4; i++ ) {
			uint32_t u = udata[i];
			if ( (u>>28) == 0 ) {
				decompressed->EncodeBit(0);
				decompressed->EncodeBits(u>>(32-bit_budget-4), bit_budget);
				total_bits += bit_budget + 1;
			} else {
				decompressed->EncodeBit(1);
				decompressed->EncodeBits(u>>(32-bit_budget), bit_budget);
				total_bits += bit_budget + 1;
			}
		}
	}
}

// Decompressor
void decompress_1d(uint8_t comp[9], BitBuffer* compressed, float* output, int bit_budget) {
	uint32_t e;
	for ( int i = 0; i < 9; i ++ ) {
		compressed->EncodeBits(comp[i], 8);
	}
	compressed->EncodeBits(0, 1);
	compressed->DecodeBits(&e, EBITS);
	int exp_max = ((int)e) - EXP_MAX;
	int dimension = 1;
	int precision_max = std::min(ZFP_MAX_PREC, std::max(0, exp_max - ZFP_MIN_EXP + (2*(dimension+1))));

	uint32_t udata[4] = {0,};
	for ( int i = 0; i < 4; i++ ) {
		uint32_t flag;
		uint32_t bits;
		compressed->DecodeBit(&flag);
		compressed->DecodeBits(&bits, bit_budget);

		if ( flag == 0 ) {
			udata[i] = (bits<<(32-bit_budget-4));
		} else {
			udata[i] = (bits<<(32-bit_budget));
		}
	}

	int32_t idata[4];

	for ( int i = 0; i < 4; i++ ) {
		idata[i] = uint2int_uint32(udata[i]);
	}
	
	inv_lift_int32(idata, 1);

	for ( int i = 0; i < 4; i++ ) {
		float q = pow(2,exp_max-(32-2));
		output[i] = ((float)idata[i]*q);
	}
}

// Compressor
void compressor(Point pointCore, Point pointTarget, uint8_t *compPoint) {
	BitBuffer* output = new BitBuffer(4*sizeof(float));
	float originPoint[4];

	originPoint[0] = pointCore.lat;
	originPoint[1] = pointCore.lon;
	originPoint[2] = pointTarget.lat;
	originPoint[3] = pointTarget.lon;
	
	compress_1d(originPoint, output, BIT_BUDGET);
	
	for ( int c = 0; c < 9; c ++ ) {
		compPoint[c] = output->buffer[c];
	}

	delete output;
}

// Decompressor
void decompressor(uint8_t compPoint[9], Point &pointCore, Point &pointTarget) {
	float decompPoint[4];
	BitBuffer* compressed = new BitBuffer(4*sizeof(float));
	
	decompress_1d(compPoint, compressed, decompPoint, BIT_BUDGET);
	delete compressed;

	pointCore.lat = decompPoint[0];
	pointCore.lon = decompPoint[1];
	pointTarget.lat = decompPoint[2];
	pointTarget.lon = decompPoint[3];
}

// Elapsed time checker
static inline double timeCheckerCPU(void) {
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

// Function for reading benchmark file
// Quadrant
void readBenchmarkDataQuadTree(Quadrant *root, char* filename, int length) {
	FILE *f_data = fopen(filename, "rb");
	if( f_data == NULL ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	for ( int i = 0; i < length; i ++ ) {
		PointQuadTree temp;
		int clusterID_org = 0;
		int datasetID = i;
		float lat = 0.00;
		float lon = 0.00;
		fread(&lat, sizeof(float), 1, f_data);
		fread(&lon, sizeof(float), 1, f_data);
		fread(&clusterID_org, sizeof(int), 1, f_data);
		temp.point.lat = lat;
		temp.point.lon = lon;
		temp.datasetID = datasetID;
		root->cities.push_back(temp);
	}

	root->done = 0;
	
	fclose(f_data);
}
// Data points
void readBenchmarkDataDBSCAN(std::vector<PointDBSCAN> &dataset, char* filename, int length) {
	FILE *f_data = fopen(filename, "rb");
	if( f_data == NULL ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	for ( int i = 0; i < length; i ++ ) {
		PointDBSCAN temp;
		int clusterID_org = 0;
		int clusterID_new = UNCLASSIFIED;
		float lat = 0.00;
		float lon = 0.00;
		fread(&lat, sizeof(float), 1, f_data);
		fread(&lon, sizeof(float), 1, f_data);
		fread(&clusterID_org, sizeof(int), 1, f_data);
		temp.point.lat = lat;
		temp.point.lon = lon;
		temp.clusterID_org = clusterID_org;
		temp.clusterID_new = clusterID_new;
		dataset.push_back(temp);
	}
	
	fclose(f_data);
}

// Manhattan
float manhattan(const Point pointCore, const Point pointTarget) {
	numManhattan++;
	float sub_lat = pointCore.lat - pointTarget.lat;
	float sub_lon = pointCore.lon - pointTarget.lon;

	return abs(sub_lat) + abs(sub_lon);
}

// Function for four edge points of square
void findEdgePointsEpsilonBox(std::vector<PointDBSCAN> &dataset) {
	for ( int i = 0; i < (int)dataset.size(); i ++ ) {
		dataset[i].northEastern.lat = dataset[i].point.lat + EPSILON;
		dataset[i].northEastern.lon = dataset[i].point.lon + EPSILON;
		dataset[i].northWestern.lat = dataset[i].point.lat - EPSILON;
		dataset[i].northWestern.lon = dataset[i].point.lon + EPSILON;
		dataset[i].southEastern.lat = dataset[i].point.lat + EPSILON;
		dataset[i].southEastern.lon = dataset[i].point.lon - EPSILON;
		dataset[i].southWestern.lat = dataset[i].point.lat - EPSILON;
		dataset[i].southWestern.lon = dataset[i].point.lon - EPSILON;
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
	float totalX = 0.00;
	float totalY = 0.00;
	for ( int i = 0; i < (int)root->cities.size(); i ++ ) {
		totalX = totalX + root->cities[i].point.lat;
		totalY = totalY + root->cities[i].point.lon;
	}
	root->center.lat = totalX / (float)root->cities.size();
	root->center.lon = totalY / (float)root->cities.size();
}

// Function for finding a length of diagonal manhattan distance
void findDiagonal(Quadrant *root) {
	root->diagonal = manhattan(root->northEastern, root->southWestern);
}

// Function for initialization
void initialize(Quadrant *root) {
	// Highest and lowest
	for ( int i = 0; i < (int)root->cities.size(); i ++ ) {
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

	// Diagonal manhattan distance
	findDiagonal(root);

	if ( root->diagonal <= EPSILON ) root->done = 1;
}

// Quadtree (Divide parent quadrant to 4 childrent quadrant)
void divideQuad(Quadrant *root) {
	for ( int i = 0; i < (int)root->cities.size(); i ++ ) {
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
void getInfoQuad(Quadrant *root) {
	for ( int i = 0; i < (int)root->child.size(); ) {
		if ( root->child[i]->cities.size() > 1 ) {
			findCenterMass(root->child[i]);
			findDiagonal(root->child[i]);
			if ( root->child[i]->diagonal <= EPSILON ) {
				root->child[i]->done = 1;
				numDataPoints = numDataPoints + (int)root->child[i]->cities.size();
			} else root->child[i]->done = 0;
			i++;
		} else if ( root->child[i]->cities.size() == 1 ) {
			findCenterMass(root->child[i]);
			findDiagonal(root->child[i]);
			root->child[i]->done = 1;
			numDataPoints = numDataPoints + (int)root->child[i]->cities.size();
			i++;
		} else {
			delete root->child[i];
			root->child.erase(root->child.begin() + i);
		}
	}
	
	// Count the total number of quadrants & delete parent's data
	if ( (int)root->child.size() > 0 ) {
		numQuadrants = numQuadrants + (int)root->child.size();
		std::vector<PointQuadTree>().swap(root->cities);
	}
}

// Quadtree (Insert new child quadrant to parent quadrant)
void insertQuad(Quadrant *root) {
	if ( root->done == 0 ) {
		// Generate child quadrants first
		root->child.resize(4);
		for ( int i = 0; i < 4; i ++ ) {
			root->child[i] = new Quadrant;
			root->child[i]->diagonal = 0.00; 
			root->child[i]->done = 0;
		}

		// Divide
		divideQuad(root);

		// Highest and lowest values of each quadrant
		findEdgePointsQuadrant(root);

		// Center mass value and diagonal distance of each quadrant
		getInfoQuad(root);
	} else return;
}

// Quadtree (Main)
int quadtree(Quadrant *root) {
	std::vector<Quadrant*> parentsCurr;
	std::vector<Quadrant*> parentsNext;
	parentsCurr.clear();
	parentsCurr.shrink_to_fit();
	parentsNext.clear();
	parentsNext.shrink_to_fit();
	int level = 0;

	// Root quadrant
	insertQuad(root);
	if ( (int)root->child.size() > 0 ) {
		for ( int i = 0; i < (int)root->child.size(); i ++ ) {
			if ( root->child[i]->done == 0 ) parentsCurr.push_back(root->child[i]);
		}
	}
	level++;

	// Iteration until meeting the terminate conditions
	while (1) {
		parentsNext.clear();
		parentsNext.shrink_to_fit();
		for ( int i = 0; i < (int)parentsCurr.size(); i ++ ) {
			insertQuad(parentsCurr[i]);
			if ( (int)parentsCurr[i]->child.size() > 0 ) {
				for ( int j = 0; j < (int)parentsCurr[i]->child.size(); j ++ ) {
					if ( parentsCurr[i]->child[j]->done == 0 ) parentsNext.push_back(parentsCurr[i]->child[j]);
				}
			}
		}
		parentsCurr.clear();
		parentsCurr.shrink_to_fit();
		if ( (int)parentsNext.size() > 0 ) {
			for ( int i = 0; i < (int)parentsNext.size(); i ++ ) {
				parentsCurr.push_back(parentsNext[i]);
			}
		} else break;
		level++;
	}
	
	return level;
}

// DBSCAN (Comparer between epsilon box and quadrant)
// Check the epsilon box is in a quadrant first
int compareEBinQ(std::vector<PointDBSCAN> &dataset, int index, Quadrant *root) {
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
int compareQinEB(std::vector<PointDBSCAN> &dataset, int index, Quadrant *root) {
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
int compareOverlap(std::vector<PointDBSCAN> &dataset, int index, Quadrant *root) {
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

// DBSCAN (Do manhattan calculation based on candidate list)
void candidateListCalculator(std::vector<PointDBSCAN> &dataset, int index, std::vector<int> &borders, Quadrant *root) {
	if ( root->done == 0 ) {
		for ( int i = 0; i < (int)root->child.size(); i ++ ) {
			candidateListCalculator(dataset, index, borders, root->child[i]);
		}
	} else {
		if ( root->diagonal <= EPSILON ) {
			for ( int i = 0; i < (int)root->cities.size(); i ++ ) {
				uint8_t compPoint[9] = {0,};
				Point pointCore = dataset[index].point;
				Point pointTarget = root->cities[i].point;
				compressor(pointCore, pointTarget, &compPoint[0]);
				decompressor(&compPoint[0], pointCore, pointTarget);
				if ( manhattan(pointCore, pointTarget) <= EPSILON ) {
					for ( int j = 0; j < (int)root->cities.size(); j ++ ) {
						borders.push_back(root->cities[j].datasetID);
					}
					break;
				}
			}
		} else {
			uint8_t compPoint[9] = {0,};
			Point pointCore = dataset[index].point;
			Point pointTarget = root->cities[0].point;
			compressor(pointCore, pointTarget, &compPoint[0]);
			decompressor(&compPoint[0], pointCore, pointTarget);
			//printf( "%f %f\n", dataset[index].point.lat, dataset[index].point.lon );
			//printf( "%f %f\n", pointCore.lat, pointCore.lon );	
			if ( manhattan(pointCore, pointTarget) <= EPSILON ) {
				borders.push_back(root->cities[0].datasetID);
			}
		}
	}
}

// DBSCAN (Do make candidate list in case of quadrant is in epsilon box)
void findQuadrantsQinEB(std::vector<PointDBSCAN> &dataset, int index, std::vector<int> &borders, Quadrant *root) {
	int resultQinEB = compareQinEB(dataset, index, root);
	if ( resultQinEB == 1 || resultQinEB == 2 || resultQinEB == 3 ) {
		if ( root->done == 0 ) {
			for ( int i = 0; i < (int)root->child.size(); i ++ ) {
				findQuadrantsQinEB(dataset, index, borders, root->child[i]);
			}
		} else {
			candidateListCalculator(dataset, index, borders, root);
		}
	} else if ( resultQinEB == 4 ) {
		candidateListCalculator(dataset, index, borders, root);
	}
}

// DBSCAN (Do make candidate list in case of epsilon box is in quadrant)
void findQuadrantsEBinQ(std::vector<PointDBSCAN> &dataset, int index, std::vector<int> &borders, Quadrant *root) {
	int resultEBinQ = compareEBinQ(dataset, index, root);
	if ( resultEBinQ != 0 ) {
		if ( root->done == 1 ) {
			candidateListCalculator(dataset, index, borders, root);
		} else {
			int resultQinEB = compareQinEB(dataset, index, root);
			if ( resultQinEB == 0 ) {
				for ( int i = 0; i < (int)root->child.size(); i ++ ) {
					findQuadrantsEBinQ(dataset, index, borders, root->child[i]);
				}
			} else findQuadrantsQinEB(dataset, index, borders, root);
		}
	} else {
		int resultQinEB = compareQinEB(dataset, index, root);
		if ( resultQinEB != 0 ) {
			findQuadrantsQinEB(dataset, index, borders, root);
		} else {
			int resultPart3 = compareOverlap(dataset, index, root);
			if ( resultPart3 != 0 ) {
				if ( root->done == 1 ) candidateListCalculator(dataset, index, borders, root);
				else {
					for ( int i = 0; i < (int)root->child.size(); i ++ ) {
						findQuadrantsEBinQ(dataset, index, borders, root->child[i]);
					}
				}
			}
		}
	}
}

// DBSCAN (Border Point Finder of Core Point)
void borderFinderCore(std::vector<PointDBSCAN> &dataset, int corePoint, std::vector<int> &bordersCore, Quadrant *root) {
	bordersCore.clear();
	bordersCore.shrink_to_fit();
	findQuadrantsEBinQ(dataset, corePoint, bordersCore, root);
}

// DBSCAN (Border Point Finder of Border Point)
void borderFinderBorder(std::vector<PointDBSCAN> &dataset, int borderPoint, std::vector<int> &bordersBorder, Quadrant *root) {
	bordersBorder.clear();
	bordersBorder.shrink_to_fit();
	findQuadrantsEBinQ(dataset, borderPoint, bordersBorder, root);
}

// DBSCAN (Cluster Expander)
int clusterExpander(std::vector<PointDBSCAN> &dataset, int index, int clusterID, Quadrant *root) {
	std::vector<int> bordersCore;
	std::vector<int> bordersBorder;
	borderFinderCore(dataset, index, bordersCore, root);

	if ( bordersCore.size() < MINIMUM_POINTS ) {
		dataset[index].clusterID_new = NOISE;
		std::vector<int>().swap(bordersCore);
		return FAILURE;
	} else {
		for ( int i = 0; i < (int)bordersCore.size(); i ++ ) {
			int borderPoint = bordersCore[i];
			dataset[borderPoint].clusterID_new = clusterID;
		}

		for ( int i = 0; i < (int)bordersCore.size(); i ++ ) {
			int borderPoint = bordersCore[i];
			if ( (dataset[borderPoint].point.lat == dataset[index].point.lat) && 
			     (dataset[borderPoint].point.lon == dataset[index].point.lon) ) {
				continue;
			} else {
				borderFinderBorder(dataset, borderPoint, bordersBorder, root);

				if ( bordersBorder.size() >= MINIMUM_POINTS ) {
					for ( int j = 0; j < (int)bordersBorder.size(); j ++ ) {
						int neighbourPoint = bordersBorder[j];
						if ( dataset[neighbourPoint].clusterID_new == UNCLASSIFIED ||
						     dataset[neighbourPoint].clusterID_new == NOISE ) {
							if ( dataset[neighbourPoint].clusterID_new == UNCLASSIFIED ) {
								bordersCore.push_back(neighbourPoint);
							}
							dataset[neighbourPoint].clusterID_new = clusterID;
						}
					}
				}
			}
		}
		std::vector<int>().swap(bordersCore);
		std::vector<int>().swap(bordersBorder);
		return SUCCESS;
	}
}

// DBSCAN (Main)
int dbscan(std::vector<PointDBSCAN> &dataset, Quadrant *root) {
	int clusterID = 1;
	for ( int i = 0; i < (int)dataset.size(); i ++ ) {
		if ( dataset[i].clusterID_new == UNCLASSIFIED ) {
			if ( clusterExpander(dataset, i, clusterID, root) != FAILURE ) {
				clusterID += 1;
				printf( "Generating cluster %d done!\n", clusterID-1 );
			}
		}
	}

	return clusterID-1;
}

// Main
int main() {
	int numCities = 700968*1;

	std::vector<PointDBSCAN> dataset;
	Quadrant *root = new Quadrant;

	// Read point data
	char benchmark_filename[] = "../worldcities_result.bin";
	readBenchmarkDataDBSCAN(dataset, benchmark_filename, numCities);
	readBenchmarkDataQuadTree(root, benchmark_filename, numCities);

	// Initialize
	initialize(root);

	// Get four edge points of epsilon box of each data point
	printf( "Finding Four Edge Points of Epsilon Box of The World Cities Start!\n" );
	double processStartStep1 = timeCheckerCPU();
	findEdgePointsEpsilonBox(dataset);
	double processFinishStep1 = timeCheckerCPU();
	double processTimeStep1 = processFinishStep1 - processStartStep1;
	printf( "Finding Four Edge Points of Epsilon Box of The World Cities Done!\n" );
	printf( "\n" );
	fflush( stdout );

	// Quadtree
	printf( "Quadtree for The World Cities Start!\n" );
	double processStartStep2 = timeCheckerCPU();
	int level = quadtree(root);
	double processFinishStep2 = timeCheckerCPU();
	double processTimeStep2 = processFinishStep2 - processStartStep2;
	printf( "Quadtree for The World Cities Done!\n" );
	printf( "\n" );
	fflush( stdout );

	// DBSCAN
	printf( "Quadtree-based DBSCAN Clustering for The World Cities Start!\n" );
	double processStartStep3 = timeCheckerCPU();
	int maxClusterID = dbscan(dataset, root);
	double processFinishStep3 = timeCheckerCPU();
	double processTimeStep3 = processFinishStep3 - processStartStep3;
	printf( "Quadtree-based DBSCAN Clustering for The World Cities Done!\n" );
	printf( "\n" );
	fflush( stdout );

	// ZFP Accuracy
	double accuracyCnt = 0;
	for ( int i = 0; i < (int)dataset.size(); i ++ ) {
		if ( dataset[i].clusterID_org == dataset[i].clusterID_new ) accuracyCnt = accuracyCnt + 1;;
	}

	// Result of Quadtree-based DBSCAN algorithm
	//printResults(dataset);
	printf( "Elapsed Time [Step1] [Epsilon Box] (CPU): %.8f\n", processTimeStep1 );
	printf( "Elapsed Time [Step2] [Quadtree] (CPU)   : %.8f\n", processTimeStep2 );
	printf( "Elapsed Time [Step3] [DBSCAN] (CPU)     : %.8f\n", processTimeStep3 );
	printf( "The Number of Data Points               : %ld\n", numDataPoints );
	printf( "The Maximum of Tree Level               : %d\n", level );
	printf( "The Number of Quadrants                 : %ld\n", numQuadrants );
	printf( "The Number of Manhattan                 : %ld\n", numManhattan );
	printf( "Max Cluster ID                          : %d\n", maxClusterID );
	printf( "The Number of Accuracy                  : %f\n", accuracyCnt );
	printf( "ZFP Accuracy [Bit Budget: %d]           : %.8f\n", BIT_BUDGET, (accuracyCnt/(double)numCities)*(double)100 );

	delete root;
	return 0;
}
