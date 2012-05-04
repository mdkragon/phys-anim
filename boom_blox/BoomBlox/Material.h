#ifndef MATERIAL_H
#define MATERIAL_H

#include <OgreMath.h>

struct Material
{
	Material(float density, float friction, float restitution, Vector3 const& color) :
		density(density), friction(friction), restitution(restitution), color(color),
			thickness(0.1f), youngsModulus(900.0f), fluidDamping(0.00001f), viscoelasticDamping(0.1f)
	{
		// compute particle mass = density * thickness * area
		pmass = density * thickness;
	}

	Material(float density, float friction, float restitution, Vector3 const& color,
				float thickness, float youngsModulus, float fluidDamping, float viscoelasticDamping) :
		density(density), friction(friction), restitution(restitution), color(color),
			thickness(thickness), youngsModulus(youngsModulus), fluidDamping(fluidDamping), viscoelasticDamping(viscoelasticDamping)
	{
		// compute particle mass = density * thickness * area
		pmass = density * thickness;
	}

	float density;
	float friction;
	float restitution;
	Vector3 color;


	// sound parameters
	float thickness;
	float youngsModulus;
	float fluidDamping;
	float viscoelasticDamping;
	// particle mass
	float pmass; 
	


	// Built-in materials, useful for testing stuff
	static Material DIRT;
	static Material BLUE_PLASTIC;
	static Material GREEN_PLASTIC;
	static Material MISSILE_STUFF;
};

#endif