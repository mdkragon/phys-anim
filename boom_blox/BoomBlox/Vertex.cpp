#include "Vertex.h"

Vertex::Vertex(Vector3 location, int id) {
	this->location = location;
	this->id = id;
	std::cout << "creating vertex at " << this->location << endl;
}

Vector3 Vertex::getLocation() {
	return this->location;
}

void Vertex::addNeighbor(Vertex * neighbor){
	neighbors.push_back(neighbor);
	std::cout << "creating edge between verticies: " << this->id << " and " << neighbor->id << endl;
}

int Vertex::getId() {
	return this->id;
}