#include "jelloMesh.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <algorithm>

#define PRINT_NORMAL_DIST(n, d) {printf("normal: (%0.3f, %0.3f, %0.3f) :: dist = %0.5f\n", n[0], n[1], n[2], d); fflush(stdout); }

#define USE_DIAG_SHEAR false 

double JelloMesh::g_structuralKs = 1000.0; 
double JelloMesh::g_structuralKd = 20.0; 
double JelloMesh::g_attachmentKs = 0.0;
double JelloMesh::g_attachmentKd = 0.0;
double JelloMesh::g_shearKs = 800.0;
double JelloMesh::g_shearKd = 10.0;
double JelloMesh::g_bendKs = 400.0;
double JelloMesh::g_bendKd = 20.0;
double JelloMesh::g_penaltyKs = 200000.0;
double JelloMesh::g_penaltyKd = 100.0;
double JelloMesh::g_shearDiagKs = 500.0;
double JelloMesh::g_shearDiagKd = 1.0;


double JelloMesh::contactThres = 0.05;

JelloMesh::JelloMesh() :     
    m_integrationType(JelloMesh::RK4), m_drawflags(MESH | STRUCTURAL),
    m_cols(0), m_rows(0), m_stacks(0), m_width(0.0), m_height(0.0), m_depth(0.0) {
  SetSize(1.0, 1.0, 1.0);
  SetGridSize(6, 6, 6);
  //SetGridSize(20, 20, 20);
}

JelloMesh::~JelloMesh() { }

void JelloMesh::Reset() {
  InitJelloMesh();
}

JelloMesh::Particle& JelloMesh::GetParticle(JelloMesh::ParticleGrid& grid, int i, int j, int k) {
  return grid[i][j][k];
}

JelloMesh::Particle& JelloMesh::GetParticle(JelloMesh::ParticleGrid& grid, int idx) {
  int i,j,k;
  GetCell(idx, i, j, k);
  return GetParticle(grid, i,j,k);
}

const JelloMesh::Particle& JelloMesh::GetParticle(const JelloMesh::ParticleGrid& grid, int i, int j, int k) const {
  return grid[i][j][k];
}

const JelloMesh::Particle& JelloMesh::GetParticle(const JelloMesh::ParticleGrid& grid, int idx) const {
  int i,j,k;
  GetCell(idx, i, j, k);
  return GetParticle(grid, i,j,k);
}

bool JelloMesh::isInterior(const JelloMesh::Spring& s) const {
  int i1,j1,k1,i2,j2,k2;
  GetCell(s.m_p1, i1, j1, k1);
  GetCell(s.m_p2, i2, j2, k2);
  return isInterior(i1,j1,k1) || isInterior(i2,j2,k2);
}


bool JelloMesh::isInterior(int idx) const {
  int i,j,k;
  GetCell(idx, i, j, k);
  return isInterior(i,j,k);
}

bool JelloMesh::isInterior(int i, int j, int k) const {
  return (i*j*k*(m_rows-i)*(m_cols-j)*(m_stacks-k) != 0);
}

void JelloMesh::SetGridSize(int cols, int rows, int stacks) {
  m_cols = cols;
  m_rows = rows;
  m_stacks = stacks;

  if (m_cols > 0 && m_rows > 0 && m_stacks > 0) {
    m_vparticles.resize(m_rows+1);
    for (int i = 0; i < m_rows+1; i++) {
      m_vparticles[i].resize(m_cols+1);
      for (int j = 0; j < m_cols+1; j++) {
        m_vparticles[i][j].resize(m_stacks+1);
      }
    }
  }
  InitJelloMesh();
}

int JelloMesh::GetGridCols() const {
  return m_cols;
}

int JelloMesh::GetGridRows() const {
  return m_rows;
}

int JelloMesh::GetGridStacks() const {
  return m_stacks;
}

void JelloMesh::SetSize(float width, float height, float depth) {
  m_width = width;
  m_height = height;
  m_depth = depth;
  InitJelloMesh();
}

float JelloMesh::GetWidth() const {
  return m_width;
}

float JelloMesh::GetHeight() const {
  return m_height;
}

float JelloMesh::GetDepth() const {
  return m_depth;
}

int JelloMesh::GetIndex(int i, int j, int k) const {
  int cols = j;
  int rows = i*(m_cols+1);
  int stacks = k*(m_cols+1)*(m_rows+1);
  return cols + rows + stacks;
}

#define ROUND(x) (floor(x + 0.5))
#define FLOOR(x) (floor(x))
#define FRACT(x) (x - FLOOR(x))
void JelloMesh::GetCell(int idx, int& i, int &j, int& k) const {
  float rows = m_rows+1;
  float cols = m_cols+1;
  float stacks = m_stacks+1;

  // derived from idx = cols*(rows*k + i) + j
  float tmp = FLOOR(idx/cols);
  j = (int) ROUND(cols*(FRACT(idx/cols)));
  i = (int) ROUND(rows*(FRACT(tmp/rows)));
  k = (int) FLOOR(tmp/rows);
}

void JelloMesh::InitJelloMesh() {
  m_vsprings.clear();

  if (m_width < 0.01 || m_height < 0.01 || m_depth < 0.01) return;
  if (m_cols < 1 || m_rows < 1 || m_stacks < 1) return;

  // Init particles
  float wcellsize = m_width / m_cols;
  float hcellsize = m_height / m_rows;
  float dcellsize = m_depth / m_stacks;
  float mass = (1.0 * 7*7*7) / ((m_cols+1)* (m_rows+1) * (m_stacks+1));
  vec3 zero(0, 0, 0);
    
  for (int i = 0; i < m_rows+1; i++) {
    for (int j = 0; j < m_cols+1; j++) {
      for (int k = 0; k < m_stacks+1; k++) {
        float x = -m_width*0.5f + wcellsize*i;
        float y = hcellsize*j; 
        float z = -m_depth*0.5f + dcellsize*k;
 
        // rotate about x
        /*
        float ry = y * cos(M_PI/4) - z * sin(M_PI/4);
        float rz = y * sin(M_PI/4) + z * cos(M_PI/4);
        y = ry;
        z = rz;
        */

        // translate off ground
        //y += 1.0;
        x += 0.0;
        y += 0.5;
        m_vparticles[i][j][k] = Particle(GetIndex(i,j,k), vec3(x, y, z), zero, mass);
      }
    }
  }

  // store copy of particles as previous
  m_vparticlesPrev = m_vparticles;


  ParticleGrid& g = m_vparticles;
  
  // Setup structural/bend springs
  for (int d = 1; d < max(max(m_rows, m_cols), m_stacks); d++) {
    for (int i = 0; i < m_rows+1; i++) {
      for (int j = 0; j < m_cols+1; j++) {
        for (int k = 0; k < m_stacks+1; k++) {
          if (j < m_cols+1-d) {
            AddStructuralSpring(GetParticle(g,i,j,k), GetParticle(g,i,j+d,k));
          }
          if (i < m_rows+1-d) { 
            AddStructuralSpring(GetParticle(g,i,j,k), GetParticle(g,i+d,j,k));
          }
          if (k < m_stacks+1-d) {
            AddStructuralSpring(GetParticle(g,i,j,k), GetParticle(g,i,j,k+d));
          }
        }
      }
    }
  }

  // Setup shear springs
  //for (int d = 1; d < max(max(m_rows, m_cols), m_stacks); d++) {
  for (int d = 1; d < max(max(m_rows, m_cols), m_stacks); d++) {
    for (int i = 0; i < m_rows+1; i++) {
      for (int j = 0; j < m_cols+1; j++) {
        for (int k = 0; k < m_stacks+1; k++) {
          // xy plane springs
          if (j < m_cols+1-d && i < m_rows+1-d) {
            // add right side spring
            AddShearSpring(GetParticle(g,i,j,k), GetParticle(g,i+d,j+d,k));
          }
          if (j > 0-1+d && i < m_rows+1-d) {
            // add left side spring
            AddShearSpring(GetParticle(g,i,j,k), GetParticle(g,i+d,j-j,k));
          }

          // yz plane
          if (k < m_stacks+1-d && j < m_cols+1-d) {
            // add right side spring
            AddShearSpring(GetParticle(g,i,j,k), GetParticle(g,i,j+d,k+d));
          }
          if (k > 0-1+d && j < m_cols+1-d) {
            // add right side spring
            AddShearSpring(GetParticle(g,i,j,k), GetParticle(g,i,j+d,k-d));
          }

          // xz plane
          if (i < m_rows+1-d && k < m_stacks+1-d) {
            // add right side spring
            AddShearSpring(GetParticle(g,i,j,k), GetParticle(g,i+d,j,k+d));
          }
          if (i > 0-1+d && k < m_stacks+1-d) {
            // add right side spring
            AddShearSpring(GetParticle(g,i,j,k), GetParticle(g,i-d,j,k+d));
          }
        }
      }
    }
  }

  if (USE_DIAG_SHEAR) {
    // Setup diag shear springs
    //for (int d = 1; d < max(max(m_rows, m_cols), m_stacks); d++) {
    for (int d = 1; d < 2; d++) {
      for (int i = 0; i < m_rows+1; i++) {
        for (int j = 0; j < m_cols+1; j++) {
          for (int k = 0; k < m_stacks+1; k++) {
            // right forward
            if (j < m_cols+1-d && i < m_rows+1-d && k < m_stacks+1-d) {
              // add right side spring
              AddShearDiagSpring(GetParticle(g,i,j,k), GetParticle(g,i+d,j+d,k+d));
            }
            // left forward
            if (j > 0-1+d && i < m_rows+1-d && k < m_stacks+1-d) {
              // add right side spring
              AddShearDiagSpring(GetParticle(g,i,j,k), GetParticle(g,i+d,j-d,k+d));
            }

            // right back 
            if (j < m_cols+1-d && i > 0-1+d && k < m_stacks+1-d) {
              // add right side spring
              AddShearDiagSpring(GetParticle(g,i,j,k), GetParticle(g,i-d,j+d,k+d));
            }
            // left back 
            if (j > 0-1+d && i > 0-1+d && k < m_stacks+d-1) {
              // add right side spring
              AddShearDiagSpring(GetParticle(g,i,j,k), GetParticle(g,i-d,j-d,k+d));
            }

          }
        }
      }
    }
  }


  // Init mesh geometry
  m_mesh.clear();
  m_mesh.push_back(FaceMesh(*this,XLEFT));
  m_mesh.push_back(FaceMesh(*this,XRIGHT));
  m_mesh.push_back(FaceMesh(*this,YTOP));
  m_mesh.push_back(FaceMesh(*this,YBOTTOM));
  m_mesh.push_back(FaceMesh(*this,ZFRONT));
  m_mesh.push_back(FaceMesh(*this,ZBACK));
}

void JelloMesh::AddStructuralSpring(Particle& p1, Particle& p2) {
  double restLen = (p1.position - p2.position).Length();
  m_vsprings.push_back(Spring(STRUCTURAL, p1.index, p2.index, g_structuralKs, g_structuralKd, restLen));
}

void JelloMesh::AddBendSpring(JelloMesh::Particle& p1, JelloMesh::Particle& p2) {
  double restLen = (p1.position - p2.position).Length();
  m_vsprings.push_back(Spring(BEND, p1.index, p2.index, g_bendKs, g_bendKd, restLen));
}

void JelloMesh::AddShearSpring(JelloMesh::Particle& p1, JelloMesh::Particle& p2) {
  double restLen = (p1.position - p2.position).Length();
  m_vsprings.push_back(Spring(SHEAR, p1.index, p2.index, g_shearKs, g_shearKd, restLen));
}

void JelloMesh::AddShearDiagSpring(JelloMesh::Particle& p1, JelloMesh::Particle& p2) {
  double restLen = (p1.position - p2.position).Length();
  m_vsprings.push_back(Spring(SHEAR, p1.index, p2.index, g_shearDiagKs, g_shearDiagKd, restLen));
}

void JelloMesh::SetIntegrationType(JelloMesh::IntegrationType type) {
  m_integrationType = type;
}

JelloMesh::IntegrationType JelloMesh::GetIntegrationType() const {
  return m_integrationType;
}

void JelloMesh::SetDrawFlags(unsigned int flags) {
  m_drawflags = flags;
}

unsigned int JelloMesh::GetDrawFlags() const {
  return m_drawflags;
}

void JelloMesh::DrawMesh(const vec3& eyePos) {
  const ParticleGrid& g = m_vparticles;
  float red[4] = {1.0,0.4,0.4,0.8};
  float white[4] = {1.0,1.0,1.0,1.0};
  float pink[4] = {0.5,0.0,0.0,1.0};
  float black[4] = {0.0,0.0,0.0,1.0};
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, red);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, black);
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, black);
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, pink);

  for (unsigned int i = 0; i < m_mesh.size(); i++) {
    m_mesh[i].CalcDistToEye(*this, eyePos);
  }
  std::sort(m_mesh.begin(), m_mesh.end(), FaceMesh::compare);
  for (unsigned int i = 0; i < m_mesh.size(); i++) {        
    m_mesh[i].Draw(*this);
  }

  //glDisable(GL_LIGHTING);
  //for (unsigned int i = 0; i < m_mesh.size(); i++)
  //{        
  //   m_mesh[i].DrawNormals(*this);
  //}

}

void JelloMesh::DrawSprings(double a) {
  const ParticleGrid& g = m_vparticles;
  glBegin(GL_LINES);
  for (unsigned int i = 0; i < m_vsprings.size(); i++) {
    if (!(m_vsprings[i].m_type & m_drawflags)) continue;
    if (isInterior(m_vsprings[i])) continue;

    switch (m_vsprings[i].m_type) {
      case BEND:       
        glColor4f(1.0, 1.0, 0.0, a); 
        break;
      case STRUCTURAL: 
        glColor4f(1.0, 1.0, 0.0, a); 
        break;
      case SHEAR:      
        glColor4f(0.0, 1.0, 1.0, a); 
        break;
    };

    vec3 p1 = GetParticle(g, m_vsprings[i].m_p1).position;
    vec3 p2 = GetParticle(g, m_vsprings[i].m_p2).position;
    glVertex3f(p1[0], p1[1], p1[2]);
    glVertex3f(p2[0], p2[1], p2[2]);
  }

  glEnd();
}

void JelloMesh::DrawCollisionNormals() {
  const ParticleGrid& g = m_vparticles;
  glBegin(GL_LINES);
  glColor3f(0.0, 1.0, 0.0);
  for(unsigned int i = 0; i < m_vcollisions.size(); i++) {
    Intersection intersection = m_vcollisions[i];
    if (isInterior(intersection.m_p)) {
      continue;
    }

    const Particle& pt = GetParticle(g, intersection.m_p);
    vec3 normal = intersection.m_normal;
    vec3 end = pt.position + 0.3 * normal;
    glVertex3f(pt.position[0], pt.position[1], pt.position[2]);
    glVertex3f(end[0], end[1], end[2]);
  }     
  /*
  for(unsigned int i = 0; i < m_vcontacts.size(); i++) {
    Intersection intersection = m_vcontacts[i];
    if (isInterior(intersection.m_p)) {
      continue;
    }

    const Particle& pt = GetParticle(g, intersection.m_p);
    vec3 normal = intersection.m_normal;
    vec3 end = pt.position + 1.0 * normal;
    glVertex3f(pt.position[0], pt.position[1], pt.position[2]);
    glVertex3f(end[0], end[1], end[2]);
  }     
  */
  glEnd();
}

void JelloMesh::DrawForces() {
  glBegin(GL_LINES);
  glColor3f(1.0, 0.0, 0.0);
  for (int i = 0; i < m_rows+1; i++) {
    for (int j = 0; j < m_cols+1; j++) {
      for (int k = 0; k < m_stacks+1; k++) {
        Particle p = m_vparticles[i][j][k];
        if (isInterior(i,j,k)) {
          continue;
        }

        vec3 normal = p.force.Normalize();
        vec3 end = p.position + 0.1 * normal;
        glVertex3f(p.position[0], p.position[1], p.position[2]);
        glVertex3f(end[0], end[1], end[2]);
      }
    }
  }     
  glEnd();
}

void JelloMesh::Draw(const vec3& eyePos) {
  if (m_drawflags & MESH) DrawMesh(eyePos);

  if (m_drawflags & (STRUCTURAL|BEND|SHEAR)) {
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glLineWidth(1.0);
    DrawSprings(0.2);
    glLineWidth(1.5);
    glEnable(GL_DEPTH_TEST);
    DrawSprings(0.4);
  }

  if (m_drawflags & NORMALS) {
    DrawCollisionNormals();
  }
  if (m_drawflags & FORCES) {
    DrawForces();
  }

  glEnable(GL_LIGHTING);
}

void JelloMesh::Update(double dt, const World& world, const vec3& externalForces) {
  m_externalForces = externalForces;

	CheckForCollisions(m_vparticles, world);
	ComputeForces(m_vparticles);
  /*
	ResolveCollisions(m_vparticles);
	ResolveContacts(m_vparticles);
  */
	ResolveContacts2(m_vparticles);
	ResolveCollisions2(m_vparticles);

  // copy the current points
  ParticleGrid currParticlesCopy = m_vparticles;

  switch (m_integrationType) {
    case EULER: 
      EulerIntegrate(dt); 
      break;
    case MIDPOINT: 
      MidPointIntegrate(dt); 
      break;
    case RK4: 
      RK4Integrate(dt); 
      break;
    case VERLET: 
      VerletIntegrate(dt); 
      break;
  }

  // set copy of the original particle positions as the previous
  m_vparticlesPrev = currParticlesCopy;
}

void JelloMesh::CheckForCollisions(ParticleGrid& grid, const World& world) {
  m_vcontacts.clear();
  m_vcollisions.clear();

  for (int i = 0; i < m_rows+1; i++) {
    for (int j = 0; j < m_cols+1; j++) {
      for (int k = 0; k < m_stacks+1; k++) {
        Particle& p = GetParticle(grid, i,j,k);

        // 1. Check collisions with world objects 
        for (unsigned int i = 0; i < world.m_shapes.size(); i++) {
          Intersection intersection;

          if (world.m_shapes[i]->GetType() == World::CYLINDER 
              && CylinderIntersection2(p, (World::Cylinder*) world.m_shapes[i], intersection)) {
            if (intersection.m_type == CONTACT) {
              m_vcontacts.push_back(intersection);
            } else if (intersection.m_type == COLLISION) {
              m_vcollisions.push_back(intersection);
            }
          } else if (world.m_shapes[i]->GetType() == World::SPHERE
              && SphereIntersection(p, (World::Sphere*) world.m_shapes[i], intersection)) {
            if (intersection.m_type == CONTACT) {
              m_vcontacts.push_back(intersection);
            } else if (intersection.m_type == COLLISION) {
              m_vcollisions.push_back(intersection);
            }
          } else if (world.m_shapes[i]->GetType() == World::CUBE
              && CubeIntersection(p, (World::Cube*) world.m_shapes[i], intersection)) {
            if (intersection.m_type == CONTACT) {
              m_vcontacts.push_back(intersection);
            } else if (intersection.m_type == COLLISION) {
              m_vcollisions.push_back(intersection);
            }
          } else if (world.m_shapes[i]->GetType() == World::GROUND 
                      && FloorIntersection(p, intersection)) {
            if (intersection.m_type == CONTACT) {
              m_vcontacts.push_back(intersection);
            } else if (intersection.m_type == COLLISION) {
              m_vcollisions.push_back(intersection);
            }
          }
        }
      }
    }
  }
  //sleep(0.2);
}

void JelloMesh::ComputeForces(ParticleGrid& grid) {
  // Add external froces to all points
  for (int i = 0; i < m_rows+1; i++) {
    for (int j = 0; j < m_cols+1; j++) {
      for (int k = 0; k < m_stacks+1; k++) {
        Particle& p = GetParticle(grid, i,j,k);
        p.force = m_externalForces * p.mass;
      }
    }
  }

  // Update springs
  for(unsigned int i = 0; i < m_vsprings.size(); i++) {
    Spring& spring = m_vsprings[i];
    Particle& a = GetParticle(grid, spring.m_p1);
    Particle& b = GetParticle(grid, spring.m_p2);

    // calculate spring force from hooks law
    // F = -(ks * (|l| - r) + kd * (ldot*l)/|l|) * (l/|l|)
    // l = a-b
    // ldot = va - vb
    vec3 l = a.position - b.position;
    float lnorm = l.Length();
    vec3 ldot = a.velocity - b.velocity;

    // proportional term
    double prop = spring.m_Ks * (lnorm - spring.m_restLen);
    // damping force
    double damp = spring.m_Kd * (Dot(ldot, l)/lnorm);
    // combined force
    vec3 force = -(prop + damp) * (l/lnorm);

    // Fa = f; Fb = -Fa;
    a.force += force;
    b.force += -force;
  }
}

void JelloMesh::ResolveContacts(ParticleGrid& grid) {
  // coefficient of restitution [0.0 - 1.0]
  // 1 - perfect bounciness (no loss of energy)
  // 0 - no bounciness (hit and stick)
  double coefOfRestitution = 0.6;

  for (unsigned int i = 0; i < m_vcontacts.size(); i++) {
    const Intersection& contact = m_vcontacts[i];
    Particle& p = GetParticle(grid, contact.m_p);
    vec3 normal = contact.m_normal; 

    double vdotn = Dot(p.velocity, normal);

    // reflect particle based on collision normal
    // project velocity onto normal
    vec3 Nproj = vdotn * normal;
    // reflect velocity (factor in bounciness of object, coefficient of restitution)
    vec3 vprime = p.velocity - 2 * Nproj * coefOfRestitution;

    // set particles velocity
    p.velocity = vprime;
  }
}

void JelloMesh::ResolveCollisions(ParticleGrid& grid) {
  for(unsigned int i = 0; i < m_vcollisions.size(); i++) {
    Intersection result = m_vcollisions[i];
    Particle& pt = GetParticle(grid, result.m_p);
    vec3 normal = result.m_normal;
    float dist = result.m_distance;

    // for now just move the particle to surface
    pt.position += dist * normal;
	}
}

void JelloMesh::ResolveContacts2(ParticleGrid& grid) {
  // coefficient of restitution [0.0 - 1.0]
  // 1 - perfect bounciness (no loss of energy)
  // 0 - no bounciness (hit and stick)
  double coefOfRestitution = 0.6;

  for (unsigned int i = 0; i < m_vcontacts.size(); i++) {
    const Intersection& contact = m_vcontacts[i];
    Particle& p = GetParticle(grid, contact.m_p);
    vec3 normal = contact.m_normal; 

    double vdotn = Dot(p.velocity, normal);

    if (vdotn < 0) {
      // reflect particle based on collision normal
      // project velocity onto normal
      vec3 Nproj = vdotn * normal;
      // reflect velocity (factor in bounciness of object, coefficient of restitution)
      vec3 vprime = p.velocity - 2 * Nproj * coefOfRestitution;

      // set particles velocity
      p.velocity = vprime;
    }

    // subtract out force into obstacle
    /*
    double ndotf = Dot(normal, p.force);
    if (ndotf < 0) {
      vec3 forceProj = ndotf * normal;
      p.force -= forceProj;
    }
    */
  }
}

void JelloMesh::ResolveCollisions2(ParticleGrid& grid) {
  for(unsigned int i = 0; i < m_vcollisions.size(); i++) {
    double coefOfRestitution = 0.6;
    Intersection result = m_vcollisions[i];
    Particle& pt = GetParticle(grid, result.m_p);
    vec3 normal = result.m_normal;
    float dist = result.m_distance;

    double vdotn = Dot(pt.velocity, normal);

    // subtract out force into obstacle
    double ndotf = Dot(normal, pt.force);
    if (ndotf < 0) {
      vec3 forceProj = ndotf * normal;
      //pt.force -= forceProj;
    }


    // add repulsion force

    // calculate spring force from hooks law
    // F = -(ks * (|l| - r) + kd * (ldot*l)/|l|) * (l/|l|)
    // l = a-b
    // ldot = va - vb
    vec3 l = dist*normal;
    double llen = l.Length();

    // proportional term
    double prop = g_penaltyKs * (llen - 0);
    // damping force
    double damp = g_penaltyKd * (Dot(pt.velocity, l)/llen);
    // combined force
    vec3 force = -(prop + damp) * (l/llen);

    // add spring collision force
    /*
    printf("------------------------------\n");
    printf("vel: (%0.3f, %0.3f, %0.3f)\n", pt.velocity[0], pt.velocity[1], pt.velocity[2]);
    printf("prop: %0.3f\n", prop);
    printf("damp: %0.3f\n", damp);
    printf("force: (%0.3f, %0.3f, %0.3f)\n", force[0], force[1], force[2]);
    printf("normal: (%0.3f, %0.3f, %0.3f)\n", normal[0], normal[1], normal[2]);
    printf("dist: %0.3f\n", dist);
    fflush(stdout);
    */
    pt.force += force;
	}
}

bool JelloMesh::FloorIntersection(Particle& p, Intersection& intersection) {
  // test if point is intersecting with the floor
  // up direction is Y and floor is at Y = 0;
  //float contactThres = 0.01;

  if (p.position[1] > contactThres) {
    return false;
  } else if (p.position[1] >= 0.0) {
    // CONTACT
    intersection.m_type = CONTACT;
  } else {
    // COLLISION
    intersection.m_type = COLLISION;
  }

  // store particle index
  intersection.m_p = p.index;
  // set collision normal
  // for floor it is just up vector
  intersection.m_normal = vec3(0, 1, 0);
  // compute distance from edge
  intersection.m_distance = p.position[1];


  return true;
}

bool JelloMesh::CubeIntersection(Particle& p, World::Cube* cube, JelloMesh::Intersection& intersection) {
  // test if point is intersecting with a box
  //float contactThres = 0.01;

  // only supports boxes aligned with the world axes
  vec3 sideHalfLen = vec3(cube->hx, cube->hy, cube->hz);
  // vector from particle to box center
  vec3 d = p.position - cube->pos;

  // outside box?
  if (abs(d[0]) > sideHalfLen[0] + contactThres
      || abs(d[1]) > sideHalfLen[1] + contactThres
      || abs(d[2]) > sideHalfLen[2] + contactThres) {
    return false;
  }

  // determine closest wall
  double ax = abs(abs(d[0]) - sideHalfLen[0]);
  double ay = abs(abs(d[1]) - sideHalfLen[1]);
  double az = abs(abs(d[2]) - sideHalfLen[2]);
  int aind = 0;
  if (ay < ax && ay < az) {
    aind = 1;
  } else if (az < ax && az < ay) {
    aind = 2;
  }

  // collsion or contact
  if (abs(d[aind]) - sideHalfLen[aind] >= 0.0) {
    // CONTACT
    intersection.m_type = CONTACT;
  } else {
    // COLLISION
    intersection.m_type = COLLISION;
  }

  // store particle index
  intersection.m_p = p.index;
  // determine wall normal
  intersection.m_normal = vec3(0,0,0);
  intersection.m_normal[aind] = d[aind]/abs(d[aind]);
  // compute distance from edge
  intersection.m_distance = abs(d[aind]) - sideHalfLen[aind];

  return true;
}



bool JelloMesh::SphereIntersection(Particle& p, World::Sphere* sphere, JelloMesh::Intersection& intersection) {
  // test if point is intersecting with the floor
  // up direction is Y and floor is at Y = 0;
  //float contactThres = 0.01;

  double r = sphere->r;
  double d = (p.position - sphere->pos).Length();

  if (d > r + contactThres) {
    return false;
  } 
  if (d >= r) {
    // CONTACT
    intersection.m_type = CONTACT;
  } else {
    // COLLISION
    intersection.m_type = COLLISION;
  }

  // store particle index
  intersection.m_p = p.index;
  // set collision normal
  // for floor it is just up vector
  intersection.m_normal = (p.position - sphere->pos).Normalize();
  // compute distance from edge
  intersection.m_distance = d - r;


  return true;
}

bool JelloMesh::CylinderIntersection(Particle& p, World::Cylinder* cylinder, JelloMesh::Intersection& intersection) {
  vec3 start = cylinder->start;
  vec3 end = cylinder->end;
  vec3 axis = end - start;
  double r = cylinder->r; 
  double cylinderLen = axis.Length();
  double sqcylinderLen = axis.SqrLength();

  // normalized axis
  vec3 naxis = axis / cylinderLen;

  //double contactThres = 0.01;
  double sqLenPlusContact = pow((cylinderLen + contactThres),2);

  // projection of p onto cylinder axis
  double pdota = Dot((p.position - start), naxis);
  vec3 pproj = pdota * naxis;
  double projLen = pproj.Length();
  // rejection of p from cylinder axis
  vec3 prej = (p.position - start) - pproj;
  // d is the distance from the particle to the axis
  double d = prej.Length();
  vec3 nprej = prej / d;

  int i, j, k;
  GetCell(p.index, i, j, k); 

  // if p is outside the radius then no contact
  if (d > r + contactThres) {
    return false;
  }
  // is p outside end caps if pdota is greater than the square of the length or negative
  if (pdota < 0.0 - pow(contactThres,2) || pdota > sqLenPlusContact) {
    return false;
  }

  // There is a contact/collision
  // store particle index
  intersection.m_p = p.index;

  // determine if it is an end cap collision
  vec3 pstart = p.position - start;
  vec3 pend = p.position - end;
  double pstartLen = pstart.Length();
  double pendLen = pend.Length();

  // simplify it as if it is in a sphere of radius r around start/end
  // then it is and endcap collision
  if (projLen < 0.0 + 10*contactThres || pdota > cylinderLen - 10*contactThres) {
    if (pstartLen < pendLen && pstartLen < r) {
      // start cap collision
      int distsign = +1;
      if (pdota < 0.0) {
        intersection.m_type = CONTACT;
        distsign = -1;
      } else {
        intersection.m_type = COLLISION;
        distsign = +1;
      }

      // set normal (negative of axis)
      intersection.m_normal = -axis / cylinderLen;
      // compute distance from edge
      intersection.m_distance = distsign * pproj.Length();
      //PRINT_NORMAL_DIST(intersection.m_normal, intersection.m_distance);
      return true;
    } else if (pendLen < pstartLen && pendLen < r) {
      // end cap collision
      if (pdota > sqcylinderLen) {
        intersection.m_type = CONTACT;
      } else {
        intersection.m_type = COLLISION;
      }

      // set normal (positive of axis)
      intersection.m_normal = axis / cylinderLen;
      // compute distance from edge
      intersection.m_distance = -(pproj.Length() - cylinderLen);
      //PRINT_NORMAL_DIST(intersection.m_normal, intersection.m_distance);
      //PRINT_NORMAL_DIST(pproj, intersection.m_distance);
      return true;
    }
  }

  // collision with wall
  if (d > r) {
    intersection.m_type = CONTACT;
  } else {
    intersection.m_type = COLLISION;
  }

  // contact with edge
  // set collision normal (rejection vector)
  intersection.m_normal = nprej;
  // compute distance from edge
  intersection.m_distance = -(d - r);

  return true;
}

bool JelloMesh::CylinderIntersection2(Particle& p, World::Cylinder* cylinder, JelloMesh::Intersection& intersection) {
  // attempt at making cylinder intersection better
  vec3 start = cylinder->start;
  vec3 end = cylinder->end;
  vec3 axis = end - start;
  double r = cylinder->r; 
  double cylinderLen = axis.Length();
  //double sqcylinderLen = axis.SqrLength();

  // normalized axis
  vec3 naxis = axis / cylinderLen;

  //double contactThres = 0.01;
  //double sqLenPlusContact = pow((cylinderLen + contactThres),2);

  // projection of p onto cylinder axis
  double pdota = Dot((p.position - start), naxis);
  vec3 pproj = pdota * naxis;
  double projLen = pproj.Length();
  // rejection of p from cylinder axis
  vec3 prej = (p.position - start) - pproj;
  // d is the distance from the particle to the axis
  double d = prej.Length();
  vec3 nprej = prej / d;

  int i, j, k;
  GetCell(p.index, i, j, k); 

  // if p is outside the radius then no contact
  if (d > r + contactThres) {
    return false;
  }
  // is p outside end caps if pdota is greater than the square of the length or negative
  if (pdota < 0.0 - contactThres || pdota > cylinderLen) {
    return false;
  }
  /*
  printf("pos(%0.2f, %0.2f, %0.2f): pdota: %0.3f :: d: %0.3f\n", p.position[0], p.position[1], p.position[2],
                                                                                pdota, d);
  fflush(stdout);
  */

  // There is a contact/collision
  // store particle index
  intersection.m_p = p.index;


  vec3 pstart = p.position - start;

  // determine which wall the particle is closest to (end cap or cylinder wall)
  // distance to cylinder edge (circular part)
  double distToEdge = abs(d - r);
  // start is at (0,0,0) in cylinder frame 
  double distToStart = pproj.Length();
  // distance from end to particle in cylinder frame
  double distToEnd = (pproj - (start - end)).Length();


  // edge contact/collision
  if (distToEdge < distToStart && distToEdge < distToEnd) {
    // collision with wall
    if (d >= r) {
      intersection.m_type = CONTACT;
    } else {
      intersection.m_type = COLLISION;
      /*
      printf("-----------------------\n");
      printf("pos: (%0.2f, %0.2f, %0.2f): pdota: %0.3f :: d: %0.3f\n", p.position[0], p.position[1], p.position[2],
                                                                                pdota, d);
      printf("start: (%0.2f, %0.2f, %0.2f): end: (%0.2f, %0.2f, %0.2f)\n", start[0], start[1], start[2], end[0], end[1], end[2]);
      printf("axis: (%0.2f, %0.2f, %0.2f): paxis: (%0.2f, %0.2f, %0.2f)\n", axis[0], axis[1], axis[2], pstart[0], pstart[1], pstart[2]);
      printf("proj: (%0.2f, %0.2f, %0.2f): rej: (%0.2f, %0.2f, %0.2f)\n", pproj[0],  pproj[1], pproj[2], prej[0], prej[1], prej[2]);
      printf("distToEdge: %0.3f, distToStart: %0.3f, distToEnd: %0.3f\n", distToEdge, distToStart, distToEnd);
      fflush(stdout);
      */
    }

    // contact with edge
    // set collision normal (rejection vector)
    intersection.m_normal = nprej;
    // compute distance from edge
    intersection.m_distance = d - r;
    return true;

  } else {
    // end cap contact/collision 
    if (distToStart < distToEnd) {
      // start cap collision
      int distsign = +1;
      if (pdota <= 0.0) {
        intersection.m_type = CONTACT;
        distsign = -1;
      } else {
        intersection.m_type = COLLISION;
        distsign = +1;
      }

      // set normal (negative of axis)
      intersection.m_normal = -axis / cylinderLen;
      // compute distance from edge
      intersection.m_distance = -pproj.Length();
      return true;
    } else {
      // end cap collision
      if (pdota >= cylinderLen) {
        intersection.m_type = CONTACT;
      } else {
        intersection.m_type = COLLISION;
      }

      // set normal (positive of axis)
      intersection.m_normal = axis / cylinderLen;
      // compute distance from edge
      intersection.m_distance = pproj.Length() - cylinderLen;
      return true;
    }
  }

  // should never reach here
  assert(false);

  return true;
}

void JelloMesh::EulerIntegrate(double dt) {
  ParticleGrid& particles = m_vparticles;  // source is a ptr!

  // Step 1
  for (int i = 0; i < m_rows+1; i++) {
    for (int j = 0; j < m_cols+1; j++) {
      for (int k = 0; k < m_stacks+1; k++) {
        // get current particle
        Particle& p = GetParticle(particles, i,j,k);

        vec3 dv = dt * p.force * (1/p.mass); 
        vec3 vel = p.velocity + dv;
        vec3 pos = p.position + dt * p.velocity;

        p.velocity = vel;
        p.position = pos;
      }
    }
  }
}

void JelloMesh::MidPointIntegrate(double dt) {
  // pointer to current particles
  ParticleGrid& particles = m_vparticles;

  // First step -- h/2 integration
  // get copy of particle grid for storage
  ParticleGrid midpoint = m_vparticles;
  for (int i = 0; i < m_rows+1; i++) {
    for (int j = 0; j < m_cols+1; j++) {
      for (int k = 0; k < m_stacks+1; k++) {
        // get current particle
        Particle& p = GetParticle(particles, i,j,k);

        vec3 dv = (dt/2) * p.force * (1/p.mass); 
        vec3 vel = p.velocity + dv;
        vec3 pos = p.position + (dt/2) * p.velocity;
        
        // store new computed values in midpoint grid
        Particle& m = GetParticle(midpoint, i,j,k);
        m.velocity = vel;
        m.position = pos;
      }
    }
  }

  // compute new forces on midpoint matrix
  ComputeForces(midpoint);

  for (int i = 0; i < m_rows+1; i++) {
    for (int j = 0; j < m_cols+1; j++) {
      for (int k = 0; k < m_stacks+1; k++) {
        // get current particle
        Particle& m = GetParticle(midpoint, i,j,k);
        Particle& p = GetParticle(particles, i,j,k);

        vec3 dv = dt * m.force * (1/m.mass); 
        vec3 vel = p.velocity + dv;
        vec3 pos = p.position + dt * m.velocity;
        
        // store final paritcle velocity and position
        p.velocity = vel;
        p.position = pos;
      }
    }
  }


}

void JelloMesh::RK4Integrate(double dt) {
    ParticleGrid target = m_vparticles;  // target is a copy!
    ParticleGrid& source = m_vparticles;  // source is a ptr!

    // Step 1
    ParticleGrid accum1 = m_vparticles;
    for (int i = 0; i < m_rows+1; i++)
    {
        for (int j = 0; j < m_cols+1; j++)
        {
            for (int k = 0; k < m_stacks+1; k++)
            {
                Particle& s = GetParticle(source, i,j,k);

                Particle& k1 = GetParticle(accum1, i,j,k);
                k1.force = dt * s.force * 1/s.mass;
                k1.velocity = dt * s.velocity;

                Particle& t = GetParticle(target, i,j,k);
                t.velocity = s.velocity + k1.force * 0.5;
                t.position = s.position + k1.velocity * 0.5;
            }
        }
    }

    ComputeForces(target);

    // Step 2
    ParticleGrid accum2 = m_vparticles;
    for (int i = 0; i < m_rows+1; i++)
    {
        for (int j = 0; j < m_cols+1; j++)
        {
            for (int k = 0; k < m_stacks+1; k++)
            {
                Particle& t = GetParticle(target, i,j,k);
                Particle& k2 = GetParticle(accum2, i,j,k);

                k2.force = dt * t.force * 1/t.mass;
                k2.velocity = dt * t.velocity;

                Particle& s = GetParticle(source, i,j,k);
                t.velocity = s.velocity + k2.force * 0.5;
                t.position = s.position + k2.velocity * 0.5;
            }
        }
    }

    ComputeForces(target);

    // Step 3
    ParticleGrid accum3 = m_vparticles;
    for (int i = 0; i < m_rows+1; i++)
    {
        for (int j = 0; j < m_cols+1; j++)
        {
            for (int k = 0; k < m_stacks+1; k++)
            {
            Particle& t = GetParticle(target, i,j,k);
            Particle& k3 = GetParticle(accum3, i,j,k);

            k3.force = dt * t.force * 1/t.mass;
            k3.velocity = dt * t.velocity;

            Particle& s = GetParticle(source, i,j,k);
            t.velocity = s.velocity + k3.force;
            t.position = s.position + k3.velocity;
            }
        }
    }
    ComputeForces(target);

    // Step 4
    ParticleGrid accum4 = m_vparticles;
    for (int i = 0; i < m_rows+1; i++)
    {
        for (int j = 0; j < m_cols+1; j++)
        {
            for (int k = 0; k < m_stacks+1; k++)
            {
                Particle& t = GetParticle(target, i,j,k);
                Particle& k4 = GetParticle(accum4, i,j,k);

                k4.force = dt * t.force * 1/t.mass;
                k4.velocity = dt * t.velocity;
            }
        }
    }

    // Put it all together
    double asixth = 1/6.0;
    double athird = 1/3.0;
    for (int i = 0; i < m_rows+1; i++)
    {
        for (int j = 0; j < m_cols+1; j++)
        {
            for (int k = 0; k < m_stacks+1; k++)
            {
                Particle& p = GetParticle(m_vparticles, i,j,k);
                Particle& k1 = GetParticle(accum1, i,j,k);
                Particle& k2 = GetParticle(accum2, i,j,k);
                Particle& k3 = GetParticle(accum3, i,j,k);
                Particle& k4 = GetParticle(accum4, i,j,k);
                
                p.velocity = p.velocity + asixth * k1.force + 
                    athird * k2.force + athird * k3.force + asixth * k4.force;

                p.position = p.position + asixth * k1.velocity + 
                    athird * k2.velocity + athird * k3.velocity + asixth * k4.velocity;
            }
        }
    }
}

void JelloMesh::VerletIntegrate(double dt) {
  ParticleGrid& particles = m_vparticles;  // source is a ptr!
  ParticleGrid& prev = m_vparticlesPrev;  // source is a ptr!

  for (int i = 0; i < m_rows+1; i++) {
    for (int j = 0; j < m_cols+1; j++) {
      for (int k = 0; k < m_stacks+1; k++) {
        // get previous particle
        Particle& pp = GetParticle(prev, i,j,k);
        // get current particle
        Particle& cp = GetParticle(particles, i,j,k);

        // compute central difference
        vec3 a = cp.force * (1/cp.mass);
        vec3 npos = 2 * cp.position - pp.position + a * pow(dt,2);

        cp.position = npos;


        // compute velocity
        vec3 vel = (cp.position - pp.position) / (2 * dt);
        cp.velocity = vel;
      }
    }
  }
}

void JelloMesh::VelocityVerletIntegrate(double dt) {
  ParticleGrid& particles = m_vparticles;  // source is a ptr!
  ParticleGrid& prev = m_vparticlesPrev;  // source is a ptr!

  ParticleGrid target = m_vparticles;

  // Velocity Verlet
  for (int i = 0; i < m_rows+1; i++) {
    for (int j = 0; j < m_cols+1; j++) {
      for (int k = 0; k < m_stacks+1; k++) {
        // get current particle
        Particle& cp = GetParticle(particles, i,j,k);
        Particle& t = GetParticle(target, i,j,k);

        // compute central difference
        vec3 a = cp.force * (1/cp.mass);
        vec3 pos = cp.position + cp.velocity * dt + (0.5) * a * pow(dt,2);

        t.position = pos;
      }
    }
  }

  // compute new forces
  ComputeForces(target);

  // compute new velocities
  for (int i = 0; i < m_rows+1; i++) {
    for (int j = 0; j < m_cols+1; j++) {
      for (int k = 0; k < m_stacks+1; k++) {
        // get current particle
        Particle& cp = GetParticle(particles, i,j,k);
        Particle& t = GetParticle(target, i,j,k);

        // compute positiion
        vec3 a = cp.force * (1/cp.mass);
        vec3 pos = cp.position + cp.velocity * dt + (0.5) * a * pow(dt,2);

        // compute accel at t+1
        vec3 a1 = t.force * (1/t.mass);
        vec3 vel = cp.velocity + (0.5) * dt * (a + a1);

        cp.position = pos;
        cp.velocity = vel;
      }
    }
  }
}



//---------------------------------------------------------------------
// Spring
//---------------------------------------------------------------------
JelloMesh::Spring::Spring() :
    m_type(JelloMesh::STRUCTURAL), 
    m_p1(-1), 
    m_p2(-1), 
    m_Ks(1.0), m_Kd(1.0), m_restLen(1.0)
{
}

JelloMesh::Spring::Spring(const JelloMesh::Spring& p) :
    m_type(p.m_type), m_p1(p.m_p1), m_p2(p.m_p2),
    m_Ks(p.m_Ks), m_Kd(p.m_Kd), m_restLen(p.m_restLen)
{
}

JelloMesh::Spring& JelloMesh::Spring::operator=(const JelloMesh::Spring& p)
{
    if (&p == this) return *this;

    m_type = p.m_type;
    m_p1 = p.m_p1;
    m_p2 = p.m_p2;
    m_Ks = p.m_Ks;
    m_Kd = p.m_Kd;
    m_restLen = p.m_restLen;
    return *this;
}

JelloMesh::Spring::Spring(JelloMesh::SpringType t, 
    int p1, int p2, double Ks, double Kd, double restLen) :
    m_type(t), m_Ks(Ks), m_Kd(Kd), m_p1(p1), m_p2(p2), m_restLen(restLen)
{
}

//---------------------------------------------------------------------
// Particle
//---------------------------------------------------------------------

JelloMesh::Particle JelloMesh::Particle::EMPTY;

JelloMesh::Particle::Particle(int idx, const vec3& p, const vec3& v, double m)
{
    index = idx;
    position = p;
    velocity = v;
    force = vec3(0,0,0);
    mass = m;
}

JelloMesh::Particle::Particle() : index(-1), position(0,0,0), velocity(0,0,0), force(0,0,0), mass(1.0)
{
}

JelloMesh::Particle::Particle(const JelloMesh::Particle& p) : 
    index(p.index), position(p.position), velocity(p.velocity), force(p.force), mass(p.mass)
{
}

JelloMesh::Particle& JelloMesh::Particle::operator=(const JelloMesh::Particle& p)
{
    if (&p == this) return *this;

    index = p.index;
    position = p.position;
    velocity = p.velocity;
    force = p.force;
    mass = p.mass;
    return *this;
}

//---------------------------------------------------------------------
// Intersection
//---------------------------------------------------------------------

JelloMesh::Intersection::Intersection() : 
    m_p(-1), m_normal(0,0,0), m_distance(0) , m_type(CONTACT)
{
}

JelloMesh::Intersection::Intersection(const JelloMesh::Intersection& p) :
    m_p(p.m_p), m_normal(p.m_normal), m_distance(p.m_distance), m_type(p.m_type)
{
}

JelloMesh::Intersection& JelloMesh::Intersection::operator=(const JelloMesh::Intersection& p)
{
    if (&p == this) return *this;
    m_p = p.m_p;
    m_normal = p.m_normal;
    m_distance = p.m_distance;
    m_type = p.m_type;
    return *this;
}

JelloMesh::Intersection::Intersection(IntersectionType type, int p, const vec3& normal, double d) :
    m_p(p), m_normal(normal), m_distance(d), m_type(type)
{
}


//---------------------------------------------------------------------
// Drawing
//---------------------------------------------------------------------

void JelloMesh::FaceMesh::Draw(const JelloMesh& m)
{
    const ParticleGrid& g = m.m_vparticles;
    for (unsigned int strip = 0; strip < m_strips.size(); strip++)
    {
        const std::vector<int>& points = m_strips[strip];

        glBegin(GL_TRIANGLE_STRIP);
        for (unsigned int pi = 0; pi < points.size(); pi++)
        {
            int idx = points[pi];
            vec3 p = m.GetParticle(g, idx).position;

            vec3 n(0,0,0);
            const std::vector<int>& neighbors = m_neighbors[idx];
            if (neighbors.size() > 0)
            {
                vec3 pup = m.GetParticle(g, neighbors[0]).position;
                vec3 pdown = m.GetParticle(g, neighbors[1]).position;
                vec3 pleft = m.GetParticle(g, neighbors[2]).position;
                vec3 pright = m.GetParticle(g, neighbors[3]).position;

                vec3 n1 = -((pright - p) ^ (pup - p));
                vec3 n2 = -((pdown - p) ^ (pright - p));
                vec3 n3 = -((pleft - p) ^ (pdown - p));
                vec3 n4 = -((pup - p) ^ (pleft - p));

                n = n1 + n2 + n3 + n4;
                n = n.Normalize();
            }

            glNormal3f(n[0], n[1], n[2]);
            glVertex3f(p[0], p[1], p[2]);
        }
        glEnd();
    }
}

void JelloMesh::FaceMesh::DrawNormals(const JelloMesh& m)
{
    glDisable(GL_LIGHTING);

    glBegin(GL_LINES);
    glColor3f(0.0, 1.0, 0.0);

    const ParticleGrid& g = m.m_vparticles;
    for (unsigned int strip = 0; strip < m_strips.size(); strip++)
    {
        const std::vector<int>& points = m_strips[strip];
        for (unsigned int pi = 0; pi < points.size(); pi++)
        {
            int idx = points[pi];
            vec3 p = m.GetParticle(g, idx).position;

            const std::vector<int>& neighbors = m_neighbors[idx];
            if (neighbors.size() == 0) continue;

            vec3 pup = m.GetParticle(g, neighbors[0]).position;
            vec3 pdown = m.GetParticle(g, neighbors[1]).position;
            vec3 pleft = m.GetParticle(g, neighbors[2]).position;
            vec3 pright = m.GetParticle(g, neighbors[3]).position;

            vec3 n1 = -((pright - p) ^ (pup - p));
            vec3 n2 = -((pdown - p) ^ (pright - p));
            vec3 n3 = -((pleft - p) ^ (pdown - p));
            vec3 n4 = -((pup - p) ^ (pleft - p));

            vec3 n = n1 + n2 + n3 + n4;
            n = n.Normalize();

            vec3 end = p + 0.2 * n;
            glVertex3f(p[0], p[1], p[2]);
            glVertex3f(end[0], end[1], end[2]);
        }
    }

    glEnd();
    glEnable(GL_LIGHTING);
}

#define R(i) max(0, min(i, m.m_rows)) // CLAMP row index
#define C(j) max(0, min(j, m.m_cols)) // CLAMP col index
#define D(j) max(0, min(j, m.m_stacks)) // CLAMP stack index
JelloMesh::FaceMesh::FaceMesh(const JelloMesh& m, JelloMesh::Face f)
{
    const ParticleGrid& g = m.m_vparticles;
    switch(f)
    {
    case ZFRONT:
        m_strips.resize(m.m_rows);
        for (int i = 0; i < m.m_rows+1; i++)
            for (int j = 0; j < m.m_cols+1; j++)
            {
                if (i < m.m_rows)
                {
                    m_strips[i].push_back(m.GetIndex(i+1,j,0));
                    m_strips[i].push_back(m.GetIndex(i,j,0));
                }

                std::vector<int> neighbors;
                neighbors.push_back(m.GetIndex(R(i), C(j+1), D(0)));
                neighbors.push_back(m.GetIndex(R(i), C(j-1), D(0)));
                neighbors.push_back(m.GetIndex(R(i-1), C(j), D(0)));
                neighbors.push_back(m.GetIndex(R(i+1), C(j), D(0)));
                m_neighbors[m.GetIndex(i,j,0)] = neighbors;
            }
        break;
    case ZBACK:
        m_strips.resize(m.m_rows);
        for (int i = 0; i < m.m_rows+1; i++)
            for (int j = 0; j < m.m_cols+1; j++)
            {
                if (i < m.m_rows)
                {
                    m_strips[i].push_back(m.GetIndex(i+1,j,m.m_stacks));
                    m_strips[i].push_back(m.GetIndex(i,j,m.m_stacks));
                }

                std::vector<int> neighbors;
                neighbors.push_back(m.GetIndex(R(i+1), C(j), D(m.m_stacks)));
                neighbors.push_back(m.GetIndex(R(i-1), C(j), D(m.m_stacks)));
                neighbors.push_back(m.GetIndex(R(i), C(j-1), D(m.m_stacks)));
                neighbors.push_back(m.GetIndex(R(i), C(j+1), D(m.m_stacks)));
                m_neighbors[m.GetIndex(i,j,m.m_stacks)] = neighbors;
            }
        break;
    case XLEFT:
        m_strips.resize(m.m_cols);
        for (int j = 0; j < m.m_cols+1; j++)
            for (int k = 0; k < m.m_stacks+1; k++)
            {
                if (j < m.m_cols)
                {
                    m_strips[j].push_back(m.GetIndex(0,j+1,k));
                    m_strips[j].push_back(m.GetIndex(0,j,k));
                }

                std::vector<int> neighbors;
                neighbors.push_back(m.GetIndex(R(0), C(j), D(k+1)));
                neighbors.push_back(m.GetIndex(R(0), C(j), D(k-1)));
                neighbors.push_back(m.GetIndex(R(0), C(j-1), D(k)));
                neighbors.push_back(m.GetIndex(R(0), C(j+1), D(k)));
                m_neighbors[m.GetIndex(0,j,k)] = neighbors;
            }
        break;
    case XRIGHT:
        m_strips.resize(m.m_cols);
        for (int j = 0; j < m.m_cols+1; j++)
            for (int k = 0; k < m.m_stacks+1; k++)
            {
                if (j < m.m_cols)
                {
                    m_strips[j].push_back(m.GetIndex(m.m_rows,j+1,k));
                    m_strips[j].push_back(m.GetIndex(m.m_rows,j,k));
                }

                std::vector<int> neighbors;
                neighbors.push_back(m.GetIndex(R(m.m_rows), C(j+1), D(k)));
                neighbors.push_back(m.GetIndex(R(m.m_rows), C(j-1), D(k)));
                neighbors.push_back(m.GetIndex(R(m.m_rows), C(j), D(k-1)));
                neighbors.push_back(m.GetIndex(R(m.m_rows), C(j), D(k+1)));
                m_neighbors[m.GetIndex(m.m_rows,j,k)] = neighbors;
            }
        break;
    case YBOTTOM:
        m_strips.resize(m.m_rows);
        for (int i = 0; i < m.m_rows+1; i++)
            for (int k = 0; k < m.m_stacks+1; k++)
            {
                if (i < m.m_rows)
                {
                    m_strips[i].push_back(m.GetIndex(i+1,0,k));
                    m_strips[i].push_back(m.GetIndex(i,0,k));
                }

                std::vector<int> neighbors;
                neighbors.push_back(m.GetIndex(R(i+1), C(0), D(k)));
                neighbors.push_back(m.GetIndex(R(i-1), C(0), D(k)));
                neighbors.push_back(m.GetIndex(R(i), C(0), D(k-1)));
                neighbors.push_back(m.GetIndex(R(i), C(0), D(k+1)));
                m_neighbors[m.GetIndex(i,0,k)] = neighbors;
            }
        break;
    case YTOP:
        m_strips.resize(m.m_rows);
        for (int i = 0; i < m.m_rows+1; i++)
            for (int k = 0; k< m.m_stacks+1; k++)
            {
                if (i < m.m_rows)
                {
                    m_strips[i].push_back(m.GetIndex(i+1,m.m_cols,k));
                    m_strips[i].push_back(m.GetIndex(i,m.m_cols,k));
                }

                std::vector<int> neighbors;
                neighbors.push_back(m.GetIndex(R(i), C(m.m_cols), D(k+1)));
                neighbors.push_back(m.GetIndex(R(i), C(m.m_cols), D(k-1)));
                neighbors.push_back(m.GetIndex(R(i-1), C(m.m_cols), D(k)));
                neighbors.push_back(m.GetIndex(R(i+1), C(m.m_cols), D(k)));
                m_neighbors[m.GetIndex(i,m.m_cols,k)] = neighbors;
            }
        break;
    }
}

void JelloMesh::FaceMesh::CalcDistToEye(const JelloMesh& m, const vec3& eyePos)
{
    std::vector<int> points = m_strips[(int) (m_strips.size()*0.5)];
    int idx = points[(int) (points.size()*0.5)];
    vec3 pos = m.GetParticle(m.m_vparticles, idx).position;
    distToEye = (pos - eyePos).Length();
}

bool JelloMesh::FaceMesh::compare(const FaceMesh& one, const FaceMesh& other)
{
    return one.distToEye > other.distToEye;
}



