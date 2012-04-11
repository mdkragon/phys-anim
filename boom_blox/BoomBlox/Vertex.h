#ifndef VERTEX_H
#define VERTEX_H

#include <OgreMath.h>
#include <vector>
using namespace std;

class Vertex;

class Vertex{
private: 
	Vector3 location; 
	vector<Vertex *> neighbors;

public:
	Vertex(Vector3 location);
	Vector3 getLocation();

	void addNeighbor(Vertex * neighbor);
};

#endif