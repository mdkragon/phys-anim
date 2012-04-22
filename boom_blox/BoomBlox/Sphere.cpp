#include "Sphere.h"
#include "Material.h"
#include <gl/glut.h>


Sphere::Sphere(float radius) : m_radius(radius)
{
	meshify(5); // can't take arguement less than 2
//	cout << getK() << endl;
}

Sphere* Sphere::Clone() const
{
	return new Sphere(*this);
}

Vector3 Sphere::GetMinExtents() const
{
	return Vector3(-GetRadius(), -GetRadius(), -GetRadius()) + GetPosition();
}

Vector3 Sphere::GetMaxExtents() const
{
	return Vector3(GetRadius(), GetRadius(), GetRadius()) + GetPosition();
}

Vector3 Sphere::GetNormalLocalInertialTensor() const
{
	return Vector3::UNIT_SCALE * (0.4f * GetNormalMass() * GetRadius() * GetRadius());
}

float Sphere::GetNormalMass() const
{
	return GetMaterial()->density * (4.0f/3.0f) * float(M_PI) * GetRadius() * GetRadius();
}

//float Sphere::SignedDistance(Vector3 const& point) const
//{
//	float distFromCenter = (point - GetPosition()).length();
//	return distFromCenter - GetRadius();
//}

void Sphere::DoRender() const
{
	glutSolidSphere(GetRadius(), 20, 20);
}

float Sphere::GetRadius() const
{
	return m_radius;
}


// gets k matrix
Eigen::MatrixXd Sphere::getK(){
	int dimension = this->verticies.size(); // there are this many verticies
	Eigen::MatrixXd K = Eigen::MatrixXd::Zero(dimension, dimension); // creates matrix
	for (int i = 0; i < dimension; i ++) {
		vector<Vertex *> neighbors = this->verticies.at(i)->getNeighbor();
		int a = this->verticies.at(i)->getId(); // a is first index
		for (int j = 0; j < neighbors.size(); j++) {
			int b = neighbors.at(j)->getId(); // b is the neighbor's id

			K(a, b) = K(a, b) + 1; 
			K(a, a) = K(a, a) - 1;
			K(b, b) = K(b, b) - 1;
		}
	}
	return K;
}

// makes a mesh
void Sphere::meshify(int divide){
	int id = 0;

	int circle_resolution = divide;
	int radius = this->m_radius;

	// push back top cap
	Vector3 base = Vector3(0, radius, 0);

	for (int i = 0; i <= circle_resolution; i++) {
		// rotates x from 0 to end
		float rot_x = (M_PI / circle_resolution) * i;
		Matrix3 rotation_x = Matrix3(1, 0, 0,
									0, cos(rot_x), -1* sin(rot_x),
									0, sin(rot_x), cos(rot_x));
		Vector3 base_rotation = rotation_x * base;

		if (i == 0) {
			// push back the base as the end cap
			verticies.push_back(new Vertex(base, id)); id = id+1;
		} else if (i == circle_resolution) {
			// push back the end cap
			Vertex * last_vertex = new Vertex(base_rotation, id); id = id+1;
			verticies.push_back(last_vertex);
			// add neighbors
			for (int j = 1 ; j <= circle_resolution; j ++) {
				last_vertex->addNeighbor(verticies.at(verticies.size() - (j+1)));
				verticies.at(verticies.size() - (j+1))->addNeighbor(last_vertex);
			}
		} else {
			// rotation around the y axis
			for (int j = 0; j < circle_resolution; j ++){
				// for each y
				float rot_y = (2 * M_PI / circle_resolution) * j;
				Matrix3 rotation_y = Matrix3 ( cos(rot_y), 0, sin(rot_y), 
												0, 1, 0,
												-1 * sin(rot_y), 0, cos(rot_y));
				Vertex * current_vertex = new Vertex(rotation_y * base_rotation, id); id = id+1;
				verticies.push_back(current_vertex);

				// adding above neighbors
				if (i==1) {
					// then we only need to add the first one
					current_vertex->addNeighbor(verticies.at(0));
					verticies.at(0)->addNeighbor(current_vertex);
				} else {
					// else we need to add the one above current
					current_vertex->addNeighbor(verticies.at( verticies.size() - (circle_resolution+1)));
					verticies.at(verticies.size() - (circle_resolution+1)) ->addNeighbor(current_vertex);
				}
				
				// one to the left
				if (j == 0) { // if first node
					// skip
				} else if (j ==circle_resolution - 1) {
					// one to the left
					current_vertex->addNeighbor(verticies.at(verticies.size() - 2));
					verticies.at(verticies.size() -2)->addNeighbor(current_vertex);

					// one to the right
					current_vertex->addNeighbor(verticies.at(verticies.size() - circle_resolution));
					verticies.at(verticies.size() - circle_resolution)->addNeighbor(current_vertex);
				} else {
					// just add the one to the left
					current_vertex->addNeighbor(verticies.at(verticies.size() - 2));
					verticies.at(verticies.size() -2)->addNeighbor(current_vertex);
				}
					
				
			}
		}
	}
}
