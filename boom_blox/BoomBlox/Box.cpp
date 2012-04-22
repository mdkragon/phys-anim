#include "Box.h"
#include "Intersection.h"
#include <gl/glut.h>
#include <TestHarness.h>
#include "Material.h"

Box::Box(Vector3 const& halfSize) : m_halfSize(halfSize)
{
	m_radius = halfSize.length();
	meshify(20);
	//cout << getK() << endl;
}

Box* Box::Clone() const
{
	return new Box(*this);
}

Vector3 Box::GetMinExtents() const
{
	// not implemented // LOOK if you need this, you'll need to write it yourself
	//assert(false);
	//return Vector3::ZERO;
	
	return GetPosition() + Vector3(-m_radius, -m_radius, -m_radius);
}

Vector3 Box::GetMaxExtents() const
{
	// not implemented // LOOK if you need this, you'll need to write it yourself
	//assert(false);
	//return Vector3::ZERO;
	return GetPosition() + Vector3(m_radius, m_radius, m_radius);
}

TEST(GetNormalLocalInertialTensor, Box)
{
	Material m(3,0,0,Vector3::ZERO);
	Box b(Vector3(1,3,2));
	b.SetMaterial(&m);
	Vector3 tensor = b.GetNormalLocalInertialTensor();
	CHECK(tensor.positionEquals(Vector3(624, 240, 480)));
}

Vector3 Box::GetNormalLocalInertialTensor() const
{
	// TODO: is this correct? I am not sure what the 'Normal' Inertial Tensor is
	// interial matrix for box is
	//				 | (y^2 + z^2)      0           0      |
	//	I = (M/12) * |      0      (x^2 + z^2)      0      |
	//				 |      0           0      (x^2 + y^2) |
	Vector3 size = 2*GetHalfSize();
	float x2 = size.x * size.x;
	float y2 = size.y * size.y;
	float z2 = size.z * size.z;
	return (GetNormalMass()/12.0) * Vector3((y2+z2), (x2+z2), (x2+y2));
}

float Box::GetNormalMass() const
{
	return GetMaterial()->density * 8 * m_halfSize.x * m_halfSize.y * m_halfSize.z;
}

void Box::DoRender() const
{
	glScalef(GetHalfSize().x, GetHalfSize().y, GetHalfSize().z);
	glutSolidCube(2);
	glPushAttrib(GL_LIGHTING_BIT);
	glDisable(GL_LIGHTING);
	glColor3f(1,1,1);
	glutWireCube(2);
	glPopAttrib();
}

Vector3 const& Box::GetHalfSize() const
{
	return m_halfSize;
}

TEST(GetExtentsAlongAxis, Box)
{
	Box box(Vector3(1,2,3));
	float a, b;

	box.GetExtentsAlongAxis(Vector4(1,0,0,0), a, b);
	FLOATS_EQUAL(-1, a);
	FLOATS_EQUAL(1, b);

	box.SetOrientation(Quaternion(M_PI_2, Vector3(0,0,1)));
	box.GetExtentsAlongAxis(Vector4(1,0,0,0), a, b);
	FLOATS_EQUAL(-2, a);
	FLOATS_EQUAL(2, b);

	box.SetOrientation(Quaternion::IDENTITY);
	box.GetExtentsAlongAxis(Vector4(1,0,0,3), a, b);
	FLOATS_EQUAL(2, a);
	FLOATS_EQUAL(4, b);

	box.GetExtentsAlongAxis(Vector4(sqrtf(0.5f),sqrtf(0.5f),0,0), a, b);
	FLOATS_EQUAL(-3*sqrtf(0.5f), a);
	FLOATS_EQUAL(3*sqrtf(0.5f), b);
}

void Box::GetExtentsAlongAxis(Vector4 const& v, float & a, float & b) const
{
	// The implementation of this function is optimized to the point of opacity.
	// If you're interested in its inner workings, contact me directly.

	// transform v into object space.
	Vector4 ov = TransformPlane(v, GetTransformation().inverseAffine());
	Vector3 ov3(ov.x, ov.y, ov.z);
	// note: ov3 is unit length

	// distance from center of box to plane
	float dc = ov.w;

	// distances along ov from center to each of the positive three faces, abs'ed and summed
	float extDist = GetHalfSize().absDotProduct(ov3);

	a = dc - extDist;
	b = dc + extDist;
}

TEST(GetClosestPoint, Box)
{
	Box box(Vector3(1,2,3));
	CHECK(box.GetClosestPoint(Vector3(0,0,0)).positionEquals(Vector3(0,0,0)));
	CHECK(box.GetClosestPoint(Vector3(1,0,0)).positionEquals(Vector3(1,0,0)));
	CHECK(box.GetClosestPoint(Vector3(2,0,0)).positionEquals(Vector3(1,0,0)));
	CHECK(box.GetClosestPoint(Vector3(5,5,0)).positionEquals(Vector3(1,2,0)));
	CHECK(box.GetClosestPoint(Vector3(5,5,5)).positionEquals(Vector3(1,2,3)));
	CHECK(box.GetClosestPoint(Vector3(-5,-5,-5)).positionEquals(Vector3(-1,-2,-3)));

	box.SetOrientation(Quaternion(M_PI_2, Vector3(0,0,1)));
	CHECK(box.GetClosestPoint(Vector3(0,0,0)).positionEquals(Vector3(0,0,0)));
	CHECK(box.GetClosestPoint(Vector3(1,0,0)).positionEquals(Vector3(1,0,0)));
	CHECK(box.GetClosestPoint(Vector3(3,0,0)).positionEquals(Vector3(2,0,0)));
	CHECK(box.GetClosestPoint(Vector3(5,5,0)).positionEquals(Vector3(2,1,0)));
	CHECK(box.GetClosestPoint(Vector3(5,5,5)).positionEquals(Vector3(2,1,3)));
	CHECK(box.GetClosestPoint(Vector3(-5,-5,-5)).positionEquals(Vector3(-2,-1,-3)));
}

Ogre::Vector3 Box::GetClosestPoint(Vector3 const& p) const
{
	// transform into local space
	Vector3 op = GetTransformation().inverseAffine()*p;
	Vector3 hs = GetHalfSize();
	if(fabs(op.x) < hs.x && fabs(op.y) < hs.y && fabs(op.z) < hs.z)
	{
		// point inside box
		return p;
	}

	Vector3 oc = Vector3(
		std::min(std::max(op.x, -hs.x), hs.x),
		std::min(std::max(op.y, -hs.y), hs.y),
		std::min(std::max(op.z, -hs.z), hs.z));

	return GetTransformation()*oc;
}


void Box::meshify(int divide) {
	// positive z is into the page
	int id = 0;

	Vector3 bound = GetHalfSize();
	float x_base = -1 * bound[0];
	float y_base = -1 * bound[1];
	float z_base = -1 * bound[2];

	Vertex * P111 = new Vertex(Vector3(x_base, bound[1], z_base), id); id = id + 1;
	Vertex * P211 = new Vertex(Vector3(x_base, y_base, z_base), id); id = id + 1;
	Vertex * P121 = new Vertex(Vector3(bound[0], bound[1], z_base), id); id = id + 1;
	Vertex * P221 = new Vertex(Vector3(bound[0], y_base, z_base), id); id = id+1;
	Vertex * P112 = new Vertex(Vector3(x_base, bound[1], bound[2]), id); id = id+1;
	Vertex * P122 = new Vertex(Vector3(bound[0], bound[1], bound[2]), id); id = id + 1;
	Vertex * P212 = new Vertex(Vector3(x_base, y_base, bound[2]), id); id = id + 1;
	Vertex * P222 = new Vertex(Vector3(bound[0], y_base, bound[2]), id); id = id+1;

	P111->addNeighbor(P112);
	P111->addNeighbor(P211);
	P111->addNeighbor(P121);

	P211->addNeighbor(P111);
	P211->addNeighbor(P212);
	P211->addNeighbor(P221);

	P121->addNeighbor(P111);
	P121->addNeighbor(P221);
	P121->addNeighbor(P122);

	P221->addNeighbor(P211);
	P221->addNeighbor(P222);
	P221->addNeighbor(P121);

	P112->addNeighbor(P122);
	P112->addNeighbor(P212);
	P112->addNeighbor(P111);

	P122->addNeighbor(P112);
	P122->addNeighbor(P222);
	P122->addNeighbor(P121);

	P222->addNeighbor(P122);
	P222->addNeighbor(P221);
	P222->addNeighbor(P212);

	P212->addNeighbor(P222);
	P212->addNeighbor(P112);
	P212->addNeighbor(P211);

	verticies.push_back(P111);
	verticies.push_back(P121);
	verticies.push_back(P221);
	verticies.push_back(P211);
	verticies.push_back(P112);
	verticies.push_back(P122);
	verticies.push_back(P222);
	verticies.push_back(P212);
}


// gets k matrix
Eigen::MatrixXd Box::getK(){
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