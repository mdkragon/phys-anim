#include <OgreMath.h>
#include <algorithm>
#include <queue>
#include <map>

#include "World.h"

#undef min
#undef max

SoundManager Sound_Manager;

World::World()
{
	m_useSweepAndPrune = true;
}

World::~World()
{
	Clear();
}

World::World(World const& other)
{
	FillFrom(other);
}

World & World::operator=(World const& other)
{
	if(this == &other) return *this;
	Clear();
	FillFrom(other);
	return *this;
}

void World::FillFrom(World const& other)
{
	std::map<Material const*, Material*> materialMapping;

	for(int i=0; i<other.GetNumMaterials(); i++)
	{
		Material const* oldMaterial = other.GetMaterial(i);
		Material* newMaterial = new Material(*oldMaterial);
		materialMapping[oldMaterial] = newMaterial;
		AddMaterial(newMaterial);
	}

	for(int i=0; i<other.GetNumBodies(); i++)
	{
		RigidBody const* oldBody = other.GetBody(i);
		RigidBody* newBody = oldBody->Clone();
		newBody->SetMaterial(materialMapping[oldBody->GetMaterial()]);
		AddBody(newBody);
		newBody->Sound_Manager = &Sound_Manager;
	}
}

void World::AddMaterial(Material* material)
{
	m_materials.push_back(material);
}

int World::GetNumMaterials() const
{
	return int(m_materials.size());
}

Material const* World::GetMaterial(int material) const
{
	return m_materials[material];
}

void World::Clear()
{
	for(int i=0; i<GetNumBodies(); i++)
	{
		delete GetBody(i);
	}

	for(int i=0; i<GetNumMaterials(); i++)
	{
		delete GetMaterial(i);
	}

	m_bodies.clear();
	m_materials.clear();

	m_xExtents.clear();
	m_yExtents.clear();
	m_zExtents.clear();
}

void World::AddBody(RigidBody* body)
{
	m_bodies.push_back(body);

	// add extents
	int bodyInd = m_bodies.size()-1;
	m_xExtents.push_back(Extent(body, bodyInd, 0, true));
	m_xExtents.push_back(Extent(body, bodyInd, 0, false));

	m_yExtents.push_back(Extent(body, bodyInd, 1, true));
	m_yExtents.push_back(Extent(body, bodyInd, 1, false));

	m_zExtents.push_back(Extent(body, bodyInd, 2, true));
	m_zExtents.push_back(Extent(body, bodyInd, 2, false));
}

int World::GetNumBodies() const
{
	return int(m_bodies.size());
}

RigidBody* World::GetBody(int body)
{
	return m_bodies[body];
}

RigidBody const* World::GetBody(int body) const
{
	return m_bodies[body];
}

int World::FindBody(RigidBody* body) const
{
	for(int i=0; i<GetNumBodies(); i++)
	{
		if(GetBody(i) == body)
		{
			return i;
		}
	}

	return -1;
}

void World::RemoveBody(int body)
{
	m_bodies.erase(m_bodies.begin()+body);

	// remove extents
	for (int i = m_xExtents.size()-1; i >= 0; i--) {
		if (m_xExtents[i].get_bodyInd() == body) {
			m_xExtents.erase(m_xExtents.begin() + i);
		}
	}
	for (int i = m_yExtents.size()-1; i >= 0; i--) {
		if (m_yExtents[i].get_bodyInd() == body) {
			m_yExtents.erase(m_yExtents.begin() + i);
		}
	}
	for (int i = m_zExtents.size()-1; i >= 0; i--) {
		if (m_zExtents[i].get_bodyInd() == body) {
			m_zExtents.erase(m_zExtents.begin() + i);
		}
	}
}

void World::Render() const
{
	for(int i=0; i<GetNumBodies(); i++)
	{
		GetBody(i)->Render();
	}
}

void World::Simulate(float dT)
{
	// LOOK this is the basis of simulation
	ResolveCollisions(dT);
	AdvanceVelocities(dT);
	ResolveContacts(dT);
	AdvancePositions(dT);
}

void World::ResolveCollisions(float dT)
{
	const int NUM_ITERATIONS = 5; // LOOK many iterations. Some things are easier to troubleshoot with only one.
	for(int ii=0; ii<NUM_ITERATIONS; ii++)
	{
		// test iterating forward
		SaveState();
		AdvanceVelocities(dT);
		AdvancePositions(dT);
		std::vector<Intersection> intersections;
		FindIntersections(intersections);
		RestoreState();
		if(intersections.empty())
		{
			break;
		}

		for(int ic=0; ic<int(intersections.size()); ic++)
		{
			Intersection &i = intersections[ic];
			float epsilon = std::min(i.bodyA->GetMaterial()->restitution, i.bodyB->GetMaterial()->restitution);
			ResolveIntersection(i, epsilon);
		}

		bool hadImpulses = false;
		for(int i=0; i<GetNumBodies(); i++)
		{
			hadImpulses = GetBody(i)->ApplyQueuedImpulses() || hadImpulses;
		}
		if(!hadImpulses)
		{
			break; // We stop early if nothing changed
		}
	}
}

TEST(AdvanceVelocities, World)
{
	Box* b = new Box(Vector3(1,1,1));
	b->SetVelocity(Vector3(0,10,0));
	World w;
	w.AddBody(b);

	w.AdvanceVelocities(1);
	CHECK(b->GetVelocity().positionEquals(Vector3(0, 0.2f, 0)));
}

void World::AdvanceVelocities(float dT)
{
	Vector3 g(0.0, -9.8, 0.0);

	// TODO compute new linear velocities and angular momentum
	for(int i = 0; i < GetNumBodies(); i++)
	{
		RigidBody *body = GetBody(i);
		if (!body->HasInfiniteMass()) {
			// P: LINEAR MOMENTUM
			// change in linear momentum is the force acting on the object
			// The only forces we have are from gravity and collisions
			// TODO: update velocity based on collisions
			body->SetVelocity(body->GetVelocity() + dT*g);

			// L: ANGULAR MOMEMTUM
			// change in angular momentum is the torque acting on the object
			//	Torques on the object will be a result of collisions
			// TODO: update angular velocity based on collision
			//body->ApplyQueuedImpulses();
		}
	}
}

void World::ResolveContacts(float dT)
{
	const int NUM_RELAXATION_STEPS = 50;
	const int NUM_PARTIAL_RELAXATION_STEPS = 10;

	for(int ir=0; ir<NUM_RELAXATION_STEPS; ir++)
	{
		SaveState();
		AdvancePositions(dT);
		std::vector<Intersection> intersections;
		FindIntersections(intersections);
		RestoreState();
		if(intersections.empty())
		{
			break;
		}

		for(int ic=0; ic<int(intersections.size()); ic++)
		{
			Intersection &i = intersections[ic];
			ResolveIntersection(i, -1 + std::min(1.0f, float(ir)/NUM_PARTIAL_RELAXATION_STEPS), true);
		}
	}
}

TEST(AdvancePositions, World)
{
	Box* b = new Box(Vector3(1,1,1));
	b->SetPosition(Vector3(0,50,0));
	b->SetVelocity(Vector3(0,-2,0));
	b->SetAngularVelocity(Vector3(0,0,0.01f));
	World w;
	w.AddBody(b);

	w.AdvancePositions(2);
	CHECK(b->GetPosition().positionEquals(Vector3(0, 46, 0)));
	float angle;
	Vector3 axis;
	b->GetOrientation().ToAngleAxis(angle, axis);
	CHECK(axis.positionEquals(Vector3(0,0,1)));
	FLOATS_EQUAL(0.02f, angle);
}

void World::AdvancePositions(float dT)
{
	for(int i=0; i<GetNumBodies(); i++)
	{
		GetBody(i)->AdvanceTransformation(dT);
	}
}

namespace
{
	Matrix3 CalcK(RigidBody const& body, Vector3 const& r)
	{
		if(body.HasInfiniteMass())
		{
			return Matrix3::ZERO;
		}
		else
		{
			Matrix3 rstar(
				0, -r.z, r.y,
				r.z, 0, -r.x,
				-r.y, r.x, 0);

			Matrix3 Iinv = body.GetInverseInertialTensor();

			return Matrix3::IDENTITY * body.GetInverseMass() + rstar.Transpose() * Iinv * rstar;
		}
	}
}

TEST(ResolveIntersection, World)
{
	Material m(1, 50, 1, Vector3(1,1,1)); // high friction
	Sphere s1(1); 
	s1.SetMaterial(&m);
	s1.SetPosition(Vector3(-3,0,0));
	s1.SetVelocity(Vector3(8,0,0));
	Sphere s2(2); 
	s2.SetMaterial(&m);
	Intersection i(&s1, &s2, Vector3(-2,0,0), Vector3(-1,0,0));
	World::ResolveIntersection(i, 1, true);

	CHECK(s1.GetVelocity().positionEquals(Vector3(-4.8f,0,0)));
	CHECK(s2.GetVelocity().positionEquals(Vector3(3.2f,0,0)));

	s1.SetPosition(Vector3(0,1,0));
	s1.SetVelocity(Vector3(1,-1,0));
	Ground g;
	g.SetMaterial(&m);
	Intersection i2(&s1, &g, Vector3(0,0,0), Vector3(0,1,0));
	World::ResolveIntersection(i2, 1, true);
	CHECK(s1.GetVelocity().positionEquals(Vector3(0.71428f,1,0)));
	CHECK(s1.GetAngularVelocity().positionEquals(Vector3(0,0,-0.71428f)));

	s1.SetVelocity(Vector3(1,-1,0));
	s1.SetAngularVelocity(Vector3::ZERO);
	m.friction = 0; // zero friction
	World::ResolveIntersection(i2, 1, true);
	CHECK(s1.GetVelocity().positionEquals(Vector3(1,1,0)));
	CHECK(s1.GetAngularVelocity().positionEquals(Vector3(0,0,0)));

	s1.SetVelocity(Vector3(1,-1,0));
	s1.SetAngularVelocity(Vector3::ZERO);
	m.friction = 0.1f; // low friction
	World::ResolveIntersection(i2, 1, true);
	CHECK(s1.GetVelocity().positionEquals(Vector3(0.8f,1,0)));
	CHECK(s1.GetAngularVelocity().positionEquals(Vector3(0,0,-0.5f)));

}

void World::ResolveIntersection(Intersection &i, float epsilon, bool immediate)
{
	// TODO I've shown you the types of the variables and a good order to calculate them in.
	RigidBody &a = *i.bodyA;
	RigidBody &b = *i.bodyB;
	Vector3 N = i.outVector;
	Vector3 loc = i.contactPoint;
	Vector3 locA = i.bodyA->GetTransformation()*i.contactPointA;
	Vector3 locB = i.bodyB->GetTransformation()*i.contactPointB;


	// relative velocity of the collision point
	// TODO: is the relative velocity correct?
	//Vector3 urel = a.GetVelocity() - b.GetVelocity();
	Vector3 urel = a.GetVelocityAtPoint(locA) - b.GetVelocityAtPoint(locB);

	float ureln = urel.dotProduct(N);
	// if ureln is positive the objects are moving away from each other?
	if(ureln > 0) {
		//std::cout << "Skipping intersection ureln: " << ureln << std::endl;
		return;
	}

	// tangential component of the relative velocity
	Vector3 urelT = urel - ureln * N;
	// T is the normalized vector of the tangential velocity
	Vector3 T = urelT.normalisedCopy();

	Matrix3 Ka = CalcK(a, locA-a.GetPosition());
	Matrix3 Kb = CalcK(b, locB-b.GetPosition());
	// TODO: I think KT is just the sum of the K's
	Matrix3 KT = Ka + Kb;
	// mu is the coefficient of friction?
	float mu = std::max(a.GetMaterial()->friction, b.GetMaterial()->friction);

	// try with sticking collision (zero tangential motion after collision)
	//  i.e. uprime_relt = 0
	//		then uprime_rel = -epsilon * ureln * N
	Vector3 uprimerel1 = -epsilon * ureln * N;

	// impulse j is calculated as u'rel = urel + KT*j
	//	(Paragraph 7, section 1)
	Vector3 j = KT.Inverse() * (uprimerel1 - urel);
	float jdotN = j.dotProduct(N);
	if((j-jdotN*N).squaredLength() <= mu*mu*jdotN*jdotN)	{
		// sticking collision; j is acceptable
	
		// play sound
		calcSound(a, b);
	} else {
		// jn = -(epsilon + 1) * ureln / (transpose(N) * K_T * (N - mu * T));		
		float jn = -(epsilon + 1) * ureln / (N.dotProduct(KT * (N - mu*T)));
		j = jn * N - mu * jn * T;	

		// play sound
		calcSound(a, b);
	}

	/*
	if ((a.id == 3 && b.id == 2) || (a.id == 3 && b.id == 2)) {
		printf("resolving intersection between %d and %d\n", a.id, b.id);
		std::cout << "j: " << j << "\tlocA: " << locA << std::endl;
		fflush(stdout);
	}
	*/

	if(immediate) {
		a.ApplyImpulse(j, locA);
		b.ApplyImpulse(-j, locB);
	} else {
		a.QueueImpulse(j, locA);
		b.QueueImpulse(-j, locB);
	}
}

void World::SaveState()
{
	for(int i=0; i<GetNumBodies(); i++)
	{
		GetBody(i)->SaveTransformation();
		GetBody(i)->SaveVelocity();
	}
}

void World::RestoreState()
{
	for(int i=0; i<GetNumBodies(); i++)
	{
		GetBody(i)->RestoreTransformation();
		GetBody(i)->RestoreVelocity();
	}
}

TEST(FindIntersections, World)
{
	Sphere s1(1); s1.SetPosition(Vector3(-2.5f,0,0));
	Sphere s2(1); s2.SetPosition(Vector3(-1.25f,-1,0));
	Sphere s3(1); s3.SetPosition(Vector3(0,0,0));
	Sphere s4(1); s4.SetPosition(Vector3(1.25f,-1,0));
	Sphere s5(1); s5.SetPosition(Vector3(2.5f,0,0));
	Sphere s6(5); s6.SetPosition(Vector3(0,0,0));
	World w;
	w.AddBody(s1.Clone());
	w.AddBody(s2.Clone());
	w.AddBody(s3.Clone());
	w.AddBody(s4.Clone());
	w.AddBody(s5.Clone());
	w.AddBody(s6.Clone());

	std::vector<Intersection> intersections;
	w.FindIntersections(intersections);
	LONGS_EQUAL(9, intersections.size());
}

void World::FindIntersections(std::vector<Intersection> & intersections)
{
	// TODO

	// LOOK this method is slow and not acceptable for a final project, but it works
	// (slowly) and may help you test other parts

	if (m_useSweepAndPrune) {
		SweepAndPrune(intersections);
	} else {
		for(int i=0; i<GetNumBodies(); i++)
		{
			for(int j=i+1; j<GetNumBodies(); j++)
			{
				FindIntersection(GetBody(i), GetBody(j), intersections);
			}
		}
	}	
}

void World::SweepAndPrune(std::vector<Intersection> & intersections) {
	// sweep and prune

	// sort the extent arrays
	// TODO: use bubble/insertion? not sure what this uses
	//printf("extent array: %d\n", m_xExtents.size());
	std::sort(m_xExtents.begin(), m_xExtents.end());
	std::sort(m_yExtents.begin(), m_yExtents.end());
	std::sort(m_zExtents.begin(), m_zExtents.end());
	//PrintExtentVector("x", m_xExtents);
	//PrintExtentVector("y", m_yExtents);
	//PrintExtentVector("z", m_zExtents);


	// create intersection array
	static Eigen::MatrixXi possibleIntersect = Eigen::MatrixXi::Zero(GetNumBodies(), GetNumBodies());
	if (possibleIntersect.cols() != GetNumBodies()) {
		possibleIntersect = Eigen::MatrixXi::Zero(GetNumBodies(), GetNumBodies());
	} else {
		possibleIntersect.setZero();
	}
	
	// update possible collision matrix from extent overlaps
	TestExtentIntersection(m_xExtents, possibleIntersect);
	TestExtentIntersection(m_yExtents, possibleIntersect);
	TestExtentIntersection(m_zExtents, possibleIntersect);

	// any pair with all axis as possible intersections need to be collision checked
	//	i.e. their matrix value is 3
	for (int i = 0; i < GetNumBodies(); i++) {
		for (int j = i+1; j < GetNumBodies(); j++) {
			if (possibleIntersect(i,j) == 3) {
				FindIntersection(GetBody(i), GetBody(j), intersections);
			}
		}
	}
}

void World::TestExtentIntersection(std::vector<Extent> extents, Eigen::MatrixXi &pintersect) {
	// increment the value for each possible intersection along each axis
	std::vector<int> exStack;
	for (int i = 0; i < extents.size(); i++) {
		Extent extent = extents[i];
		int bodyInd = extent.get_bodyInd();

		if (extent.is_beginning()) {
			// push onto stack
			exStack.push_back(bodyInd);
		} else {
			// there is a possible collision with any objects that are on the stack
			int begInd = 0;
			for (int j = 0; j < exStack.size(); j++) {
				if (exStack[j] == bodyInd) {
					begInd = j;
				} else {
					pintersect(bodyInd,exStack[j]) += 1;
					pintersect(exStack[j],bodyInd) += 1;
				}
			}
			// pop the beginning from the stack
			exStack.erase(exStack.begin() + begInd);
		}
	}
}

void World::PrintExtentVector(const char *label, std::vector<Extent> extents) {
	printf("%s:: ", label);
	for (int i = 0; i < extents.size(); i++) {
		Extent extent = extents[i];
		printf("%s:%d:%5.2f  ", extent.is_beginning() ? "b" : "e", 
								extent.get_bodyInd(),
								extent.get_extent());		
	}
	printf("\n");
}

void World::SetUseSweepAndPrune(bool useSweepAndPrune) {
	m_useSweepAndPrune = useSweepAndPrune;
}

void World::calcSound(RigidBody &a, RigidBody &b) {
	Ground * g; 

	g = dynamic_cast<Ground *>(&a);
	if (g == 0) {  // g == 0 then cast not successful
		a.calculateSound(&Sound_Manager);
	} 
	g = dynamic_cast<Ground * > (&b);
	if (g== 0) {
		b.calculateSound(&Sound_Manager);
	}
	
	//b.calculateSound(Sound_Manager);
	/* sound_manager now saves the FMOD::Sound objects/structs in a vector
	(it is slow right now)
	 and when you do play user create sounds, it plays everything in the vector

	 I did this to test that the sound object is what stores the information. 
	 And that it can handle multiple sound objects. 

	 I haven't had a chance to hack pcmreadcallback yet but from what I understand:

	 FMOD_CREATESOUNDEXINFO is the system that they save the all the information related
	 to the user generated sound.

	 result = system->createSound(0, mode, &createsoundexinfo, &tmpsound);

	 here, createSound will take the info from createsoundexinfo, and put
	 the sounddata into an object of type FMOD::SOUND. The location of that
	 object is the last parameter (in this case tmpsound)
	
	 to create the sound, createSound knows to use pcmreadcallback. 
	 pcmreadcallback(FMOD_SOUND *sound, void *data, unsigned int datalen)

	 you get all these params, but you don't necessarily need to use it.
	 data is where the sound buffer lives. datalen is how long the data stream is

	 so if somehow we can put into void *data, the information generated by 
	 the math, then we are good to go. (in other words, we don't need to create new 
	 instances of SoundManager. We only need the program to generate new FMOD::Sound 
	 objects

	 I'm thinking either some kind of global or maybe soundmanager can have a 
	 current_sound array pointer or something. 
	 (this'll probably work since we are only processing one colission at a time).

	 I know you mentioned something about not haveing to process the whole sound at one time.
	 I'm not sure how to do that. But I'm not sure if we want to deal with it, since 
	 sound manager seems to take care of that for us. 

	 I think the next step could be to take that file we generated, and put that 
	 into data* and see if it sounds similar. 

	 let me know what you think
	 */


	// creates 3 instances of user-gen sound.
	//Sound_Manager.InitUserCreatedSound();
	//Sound_Manager.InitUserCreatedSound();
	//Sound_Manager.InitUserCreatedSound();
	// plays all of them
	//Sound_Manager.PlayUserCreatedSound();

	//TODO , figure out when to release sound streams
	// ie, maybe each object can only have max n sound at a time or something like that
	// or some kind of global max sounds thing
}
