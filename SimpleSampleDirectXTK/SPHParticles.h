#pragma once

#include <Eigen/Dense>
#include <vector>
#include <map>
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include "WavePacket.h"
using namespace std;
using namespace Eigen;

struct Particle {
	//constructor for particles so we can initialize them with velocities
	Particle(float _x, float _z, float _vx, float _vz) : x(_x, 0.05f,_z), v(_vx,0.f, _vz), f(0.f, 0.f,0.f), rho(0), p(0.f) {}
	Vector3f x, v, f;
	Vector2f x2D() {
		return Vector2f(x.x(), x.z());
	}
	float rho, p;
};
// solver parameters
const static Vector3f G(0.f, -9.8f ,0.f); // external (gravitational) forces
const static float REST_DENS = 1.f; // rest density
const static float GAS_CONST = 2000.f; // const for equation of state
const static float H = 16.f; // kernel radius
const static float HSQ = H*H; // radius^2 for optimization
const static float MASS = 90.f; // assume all particles have the same mass
const static float VISC = 250.f; // viscosity constant
const static float DT = 0.016f; // integration timestep

const static float SPLASH_LIFETIME = 1.f;

const static float surfaceTensionConstant = 32.f / (M_PI * std::pow(H, 9.f));
const static float surfaceTensionOffset = -std::pow(H, 6.f) / 64.f;

								 // smoothing kernels defined in Müller and their gradients
const static float POLY6 = 315.f / (65.f*M_PI*pow(H, 9.f));
const static float SPIKY_GRAD = -45.f / (M_PI*pow(H, 6.f));
const static float VISC_LAP = 45.f / (M_PI*pow(H, 6.f));

// simulation parameters
const static float EPS = H; // boundary epsilon
const static float BOUND_DAMPING = -0.6f;

const static int DAM_PARTICLES = 500;

// rendering projection parameters
const static int WINDOW_WIDTH = 800;
const static int WINDOW_HEIGHT = 600;
const static double VIEW_WIDTH = 1.5*800.f;
const static double VIEW_HEIGHT = 1.5*600.f;
struct SplashContainer {
	int val = 0;
	float time = 0.f;
	vector<Particle> particles;
};
class Particles {
public:
	void AddSplash(WAVE_PACKET* packet);
	void RemoveSplash(WAVE_PACKET* packet);
	void Update(int p_groundSizeY, int p_groundSizeX, float* p_distMap);
	bool IsSplashDone(WAVE_PACKET* packet);
	map<WAVE_PACKET*, SplashContainer> splashes;
	int			m_groundSizeX, m_groundSizeY;	// pixel size of the ground texture
private:
	float *m_distMap;
	float GetBoundaryDist(Vector2f &p);
	void Integrate();
	void ComputeDensityPressure(void);
	void ComputeForces(void);
	void CheckForDoneSplashes();
};
