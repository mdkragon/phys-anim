#include "Extent.h"
#include <OgreMath.h>


Extent::Extent(RigidBody* body, int bodyInd, int axisInd, bool isBeginning) {
	m_body = body;
	m_bodyInd = bodyInd;
	m_axisInd = axisInd;
	m_isBeginning = isBeginning;
}

Extent::~Extent() {
}

const float Extent::get_extent() {
	Vector3 p;
	if (m_isBeginning) {
		p = m_body->GetMinExtents();
	} else {
		p = m_body->GetMaxExtents();
	}
	return p[m_axisInd];
}

int Extent::get_bodyInd() {
	return m_bodyInd;
}

bool Extent::is_beginning() {
	return m_isBeginning;
}
