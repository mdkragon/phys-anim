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
	int id;

public:
	Vertex(Vector3 location, int id);
	Vector3 getLocation();

	int getId();
	void addNeighbor(Vertex * neighbor);
	vector<Vertex *> getNeighbor();
};

#endif