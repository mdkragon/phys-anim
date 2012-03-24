#ifndef EXTENT_H
#define EXTENT_H

#include <vector>
#include <string>
#include "RigidBody.h"


class RigidBody;

class Extent
{
public:
	Extent(RigidBody* body, int bodyInd, int axisInd, bool isBeginning);
	~Extent();

	const float get_extent();
	int Extent::get_bodyInd();
	bool Extent::is_beginning();

	inline bool operator < ( Extent& rhs ) {
		if (get_extent() < rhs.get_extent()) {
			return true;
		} else {
			return false;
		}		
	}

	inline bool operator > ( Extent& rhs ) {
		if (get_extent() > rhs.get_extent()) {
			return true;
		} else {
			return false;
		}		
	}

private:
	RigidBody *m_body;
	int m_bodyInd;
	int m_axisInd;
	bool m_isBeginning;
};

#endif