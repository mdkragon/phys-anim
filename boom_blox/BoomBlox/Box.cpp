#include "Box.h"
#include "Intersection.h"
#include <gl/glut.h>
#include <TestHarness.h>
#include "Material.h"

Box::Box(Vector3 const& halfSize) : m_halfSize(halfSize)
{
	m_radius = halfSize.length();
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

