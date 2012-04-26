#include "RigidBody.h"
#include <limits>
#include <gl/glut.h>
#include "Material.h"
#include <TestHarness.h>

namespace
{
	float inf = std::numeric_limits<float>::infinity();
}

RigidBody::RigidBody() :
	m_transformation(Matrix4::IDENTITY),
	m_velocity(Vector3::ZERO),
	m_angularVelocity(Vector3::ZERO),
	m_material(NULL),
	m_queuedDeltaVelocity(Vector3::ZERO),
	m_queuedDeltaAngularVelocity(Vector3::ZERO),
	m_hasInfiniteMass(false),
	id(0)
{
}

void RigidBody::SetID(int id_)
{
	id = id_;
}

Material const* RigidBody::GetMaterial() const
{
	return m_material;
}

void RigidBody::SetMaterial(Material const* material)
{
	m_material = material;
}

Vector3 RigidBody::GetNormalLocalInverseInertialTensor() const
{
	Vector3 tensor = GetNormalLocalInertialTensor();
	return Vector3(1/tensor.x, 1/tensor.y, 1/tensor.z);
}

Vector3 RigidBody::GetPosition() const
{
	return m_transformation * Vector3::ZERO;
}

void RigidBody::SetPosition(Vector3 const& position)
{
	m_transformation.setTrans(position);
}

Quaternion RigidBody::GetOrientation() const
{
	return m_transformation.extractQuaternion();
}

void RigidBody::SetOrientation(Quaternion const& orientation)
{
	m_transformation.makeTransform(GetPosition(), Vector3::UNIT_SCALE, orientation);
}


Matrix4 const& RigidBody::GetTransformation() const
{
	return m_transformation;
}

void RigidBody::SetTransformation(Matrix4 const& transformation)
{
	m_transformation = transformation;
}

Vector3 const& RigidBody::GetVelocity() const
{
	return m_velocity;
}

void RigidBody::SetVelocity(Vector3 const& velocity)
{
	m_velocity = velocity;
}

Ogre::Vector3 RigidBody::GetVelocityAtPoint(Vector3 const& p) const
{
	Vector3 v = p-GetPosition();
	return GetAngularVelocity().crossProduct(v) + GetVelocity();
}

Vector3 const& RigidBody::GetAngularVelocity() const
{
	return m_angularVelocity;
}

void RigidBody::SetAngularVelocity(Vector3 const& angularVelocity)
{
	m_angularVelocity = angularVelocity;
}

void RigidBody::QueueImpulse(Vector3 const& impulse, Vector3 const& p)
{
	if(HasInfiniteMass()) return;

	Vector3 deltaVelocity = impulse * GetInverseMass();
	m_queuedDeltaVelocity += deltaVelocity;

	Vector3 deltaOmega = GetInverseInertialTensor() * (p-GetPosition()).crossProduct(impulse);
	m_queuedDeltaAngularVelocity += deltaOmega;
}

bool RigidBody::ApplyQueuedImpulses()
{
	if(HasInfiniteMass()) return false;

	bool significant = (m_queuedDeltaVelocity.squaredLength() > 1e-8 || m_queuedDeltaAngularVelocity.squaredLength() > 1e-8);

	if(significant)
	{
		SetVelocity(GetVelocity() + m_queuedDeltaVelocity);
		SetAngularVelocity(GetAngularVelocity() + m_queuedDeltaAngularVelocity);
	}


	m_queuedDeltaVelocity = Vector3::ZERO;
	m_queuedDeltaAngularVelocity = Vector3::ZERO;

	return significant;
}

void RigidBody::ApplyImpulse(Vector3 const& impulse, Vector3 const& point)
{
	if(HasInfiniteMass()) return;

	Vector3 deltaVelocity = impulse * GetInverseMass();
	SetVelocity(GetVelocity() + deltaVelocity);

	Vector3 deltaOmega = GetInverseInertialTensor() * (point-GetPosition()).crossProduct(impulse);
	SetAngularVelocity(GetAngularVelocity() + deltaOmega);
}

void RigidBody::SaveTransformation()
{
	m_savedTransformation = m_transformation;
}

void RigidBody::RestoreTransformation()
{
	m_transformation = m_savedTransformation;
}

void RigidBody::SaveVelocity()
{
	m_savedVelocity = m_velocity;
	m_savedAngularVelocity = m_angularVelocity;
}

void RigidBody::RestoreVelocity()
{
	m_velocity = m_savedVelocity;
	m_angularVelocity = m_savedAngularVelocity;
}

void RigidBody::AdvanceTransformation(float dT)
{
	Vector3 newPosition = GetPosition() + GetVelocity()*dT;
	float angularSpeed = GetAngularVelocity().length();
	Vector3 axis = GetAngularVelocity().normalisedCopy();
	Quaternion newOrientation = Quaternion(angularSpeed * dT, axis) * GetOrientation();
	newOrientation.normalise();

	m_transformation.makeTransform(newPosition, Vector3::UNIT_SCALE, newOrientation);
}

void RigidBody::AdvanceVelocity(Vector3 const& F, float dT)
{
	if(!HasInfiniteMass())
	{
		SetVelocity(GetVelocity() + F * dT);
	}
}

void RigidBody::SetInfiniteMass()
{
	m_hasInfiniteMass = true;
}

void RigidBody::UnsetInfiniteMass()
{
	m_hasInfiniteMass = false;
}

bool RigidBody::HasInfiniteMass() const
{
	return m_hasInfiniteMass;
}

Vector3 RigidBody::GetLocalInertialTensor() const
{
	return HasInfiniteMass() ? Vector3(inf, inf, inf) : GetNormalLocalInertialTensor();
}

Vector3 RigidBody::GetLocalInverseInertialTensor() const
{
	return HasInfiniteMass() ? Vector3(0, 0, 0) : GetNormalLocalInverseInertialTensor();
}

float RigidBody::GetMass() const
{
	return HasInfiniteMass() ? inf : GetNormalMass();
}

float RigidBody::GetInverseMass() const
{
	return HasInfiniteMass() ? inf : GetNormalInverseMass();
}

Matrix3 RigidBody::GetInertialTensor() const
{
	Vector3 diagonal = GetLocalInertialTensor();
	Matrix3 localTensor = Matrix3::ZERO;
	localTensor[0][0] = diagonal.x;
	localTensor[1][1] = diagonal.y;
	localTensor[2][2] = diagonal.z;

	Matrix3 rotation;
	GetTransformation().extract3x3Matrix(rotation);

	return rotation * localTensor * rotation.Transpose();
}

Matrix3 RigidBody::GetInverseInertialTensor() const
{
	Vector3 diagonal = GetLocalInverseInertialTensor();
	Matrix3 localInverseTensor = Matrix3::ZERO;
	localInverseTensor[0][0] = diagonal.x;
	localInverseTensor[1][1] = diagonal.y;
	localInverseTensor[2][2] = diagonal.z;

	Matrix3 rotation;
	GetTransformation().extract3x3Matrix(rotation);

	return rotation.Transpose() * localInverseTensor * rotation;
}

void RigidBody::Render() const
{
	glColor3f(GetMaterial()->color.x, GetMaterial()->color.y, GetMaterial()->color.z);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	
	Vector3 pos = GetPosition();
	glTranslatef(pos.x, pos.y, pos.z);
	float angle;
	Vector3 axis;
	GetOrientation().ToAngleAxis(angle, axis);
	glRotatef(angle*180.0f/M_PI, axis.x, axis.y, axis.z);

	DoRender();
	glPopMatrix();
	int err = glGetError();
	err = err;
}


void RigidBody::getK(Eigen::MatrixXd &K){
	int dimension = this->verticies.size(); // there are this many verticies
	Eigen::MatrixXd K_Matrix = Eigen::MatrixXd::Zero(dimension, dimension); // creates matrix
	for (int i = 0; i < dimension; i ++) {
		vector<Vertex *> neighbors = this->verticies.at(i)->getNeighbor();
		int a = this->verticies.at(i)->getId(); // a is first index
		for (int j = 0; j < neighbors.size(); j++) {
			int b = neighbors.at(j)->getId(); // b is the neighbor's id

			// twice diagaonal 
			// positive K- constant
			K_Matrix(a, b) = K_Matrix(a, b) - 1; 
			K_Matrix(a, a) = K_Matrix(a, a) + 1;
			K_Matrix(b, b) = K_Matrix(b, b) + 1;
		}
	}

	
	// now the K is created for one dimention
	// so we need to make the K matrix 3N x 3N

	Eigen::MatrixXd K_expand = Eigen::MatrixXd::Zero(3*dimension, 3*dimension); 
	K_expand << K_Matrix , Eigen::MatrixXd::Zero(dimension, dimension), Eigen::MatrixXd::Zero(dimension, dimension),
		Eigen::MatrixXd::Zero(dimension, dimension), K_Matrix, Eigen::MatrixXd::Zero(dimension, dimension),
		Eigen::MatrixXd::Zero(dimension, dimension), Eigen::MatrixXd::Zero(dimension, dimension), K_Matrix;

	// now we need to multiply it by the k-constant

	float thickness = 5;
	float Y = 200; // young's modulus of steel
	float k_constant = thickness * Y;

	K_expand = k_constant * K_expand;

	K = K_expand;

	return;
}

// returns the mass matrix
void RigidBody::getM(Eigen::MatrixXd &M) {
	int dimension = this->verticies.size(); // there are this many verticies
	dimension = 3 * dimension;
	// m = density * thickness * area; 
	// for now assume homogenous object
	// density for steel is 7.85 g/cm^3

	double thickness = 5;
	double mass = 0.785 * 3 * thickness;
	mass = 0.1; // over write for now

	Eigen::MatrixXd mass_matrix = mass * Eigen::MatrixXd::Identity(dimension, dimension);

	M = mass_matrix;

	return;
}

void RigidBody::diagonalizeK(const Eigen::MatrixXd &K, Eigen::MatrixXd &G,
										Eigen::MatrixXd &D, Eigen::MatrixXd &Ginv) {
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(K);
	if (eigensolver.info() != Eigen::Success) abort();

	D = eigensolver.eigenvalues().asDiagonal(); // eigenvalues.  G
	G = eigensolver.eigenvectors(); // eigen vectors. D
	Ginv = eigensolver.eigenvectors().inverse(); // Ginv

	return;
}

void RigidBody::constructW(Eigen::VectorXcd &W_plus, Eigen::VectorXcd &W_minus){
	double gama = 0.00001; // fluid damping
	double n_eta = 0.1; // viscolastic damping

	int dimension = verticies.size();
	dimension = dimension * 3;

	W_plus = Eigen::VectorXcd::Zero(dimension);
	W_minus = Eigen::VectorXcd::Zero(dimension);

	// imaginary numbers lol ... 
	// puts complex cast every where and hope it works

	// w_i = ( -(gama*lambda_i + n) +/- sqrt ((gamma * lambda_i + n)^2 - 4*lambda_i))/2
	for (int i = 0; i < dimension; i ++) {
		double lambda_i = D(i,i);
		complex<double> p1 = -1 * (gama * lambda_i + n_eta);
		complex<double> p2 = sqrt( complex<double> (pow(gama * lambda_i + dimension, 2) - 4 * lambda_i )); 
		
		W_plus(i) = complex<double> ((p1 + p2)/2.0);
		W_minus(i) = complex<double> ((p1 + p2)/2.0);
	}
}


void RigidBody::initSoundScene() {
	// call meshify
	this->meshify(5);
	// first calcualte K
	this->getK(K);
	//std::cout << K << endl;
	// then get M
	this->getM(M);
	this->diagonalizeK(K, G, D, Ginv);
	// construct the W
	constructW(W_plus, W_minus);
	std::cout << "this is plus "<< W_plus << endl;
	std::cout << "this is minus " << W_minus << endl;
}
