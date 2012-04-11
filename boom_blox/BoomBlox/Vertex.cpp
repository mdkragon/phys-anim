#include "Vertex.h"

Vertex::Vertex(Vector3 location) {
	this->location = location;
}

Vector3 Vertex::getLocation() {
	return this->location;
}

void Vertex::addNeighbor(Vertex * neighbor){
	neighbors.push_back(neighbor);
}