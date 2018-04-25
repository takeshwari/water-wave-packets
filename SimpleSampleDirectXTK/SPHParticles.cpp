#include "SPHParticles.h"
#include "GlobalDefs.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
//This class keeps track of all SPH-simulated water particles in a scene and simulates them each update step
void Particles::Update(int p_groundSizeY, int p_groundSizeX, float* p_distMap) {
	m_groundSizeX = p_groundSizeX;
	m_groundSizeY = p_groundSizeY;
	m_distMap = p_distMap;
	ComputeDensityPressure();
	ComputeForces();
	//Simulation disabled until red dot graphics are debugged
	Integrate();
	CheckForDoneSplashes();
}

inline float Particles::GetBoundaryDist(Vector2f &p)
{
	Vector2f pTex = Vector2f(p.x() / SCENE_EXTENT + 0.5f, p.y() / SCENE_EXTENT + 0.5f);		// convert from world space to texture space
	float val1 = m_distMap[(int)(max(0, min(m_groundSizeY - 1, (int) pTex.y()*m_groundSizeY)))*m_groundSizeX + (int)(max(0, min(m_groundSizeX - 1,(int) pTex.x()*m_groundSizeX)))];
	float val2 = m_distMap[(int)(max(0, min(m_groundSizeY - 1, (int)pTex.y()*m_groundSizeY)))*m_groundSizeX + (int)(max(0, min(m_groundSizeX - 1, 1 + (int)pTex.x()*m_groundSizeX)))];
	float val3 = m_distMap[(int)(max(0, min(m_groundSizeY - 1, 1 + (int)pTex.y()*m_groundSizeY)))*m_groundSizeX + (int)(max(0, min(m_groundSizeX - 1, (int)pTex.x()*m_groundSizeX)))];
	float val4 = m_distMap[(int)(max(0, min(m_groundSizeY - 1, 1 + (int)pTex.y()*m_groundSizeY)))*m_groundSizeX + (int)(max(0, min(m_groundSizeX - 1, 1 + (int) pTex.x()*m_groundSizeX)))];
	float xOffs = (pTex.x()*m_groundSizeX) - (int)(pTex.x()*m_groundSizeX);
	float yOffs = (pTex.y()*m_groundSizeY) - (int)(pTex.y()*m_groundSizeY);
	float valH1 = (1.0f - xOffs)*val1 + xOffs*val2;
	float valH2 = (1.0f - xOffs)*val3 + xOffs*val4;
	return((1.0f - yOffs)*valH1 + yOffs*valH2);
}

//iterates through all splashes and checks if that splash is done. If so, remove it!
void Particles::CheckForDoneSplashes() {
	auto iter = splashes.begin();
	while (iter != splashes.end()) {
		bool isDone = false;

		SplashContainer splash = iter->second;
		bool isFlat = true;
		for (auto &p : splash.particles) {
			if (p.v.norm() > 0.6f) {
				isFlat = false;
			}
		}

		isDone = isFlat || splash.time > SPLASH_LIFETIME;
		if (isDone) {
			iter = splashes.erase(iter);
		}
		else {
			iter++;
		}
	}
}

//checks if there is no longer a splash entry for the given wave packet

bool Particles::IsSplashDone(WAVE_PACKET* packet) {
	return splashes.find(packet) == splashes.end();
}

//given a wave packet, generate a splash container for it and add it to the simulation
void Particles::AddSplash(WAVE_PACKET* packet) {
	SplashContainer splash = SplashContainer();
	//Particle p = Particle(packet->pos1.x, packet->pos1.y(), 0.f, 0.f);
	//PLACEHOLDER - just spawns one particle per packet! Used to debug red dot graphics
	/*for (int x = 0; x < 10; x++) {
		for (int z = 0; z < 10; z++) {
			splash.particles.push_back(Particle(5.f + (float)x, 30.f+ (float)z, 0.f, 0.f));
		}
	}*/
	//splash.particles.push_back(Particle(0.f, 0.f, 0.f, 0.f));
	

	// Plan out algorithm for packet here:

	/*
	
		1. Find packet width (width)
		2. Find packet length (length)
		3. Find packet forward vector (vF)
		4. Find packet right vector (vR)
		5. Find front left corner point of packet 9	(flCornerPt)
		6. Set spacing const (p)
		7. Find numCols (nC = width/p)
		8. Find numRows (nR = length/p)
		9. Find packet speed (s)
		10. Find particle velocities (pVel)

		- frontEdge is defined by:	 flCornerPt + vR * p * i [i from 0 -> numCols]
		- trailingEdge(pt) is defined by: pt - vF * p * i [i from 0 -> numRows]

		6. FOREACH (point pt in frontEdge):
				FOR (point trailPt in trailingEdge(pt)):
					add particle (pos=trailPt, vel= vF * s
	*/
	float width = packet->envelope * 2;
	float length = packet->envelope;
	Vector2f vF = packet->travelDir;
	Vector2f vR = Vector2f(vF.y(), -1 * vF.x());
	Vector2f flCornerPt = packet->midPos - packet->envelope * vR;
	float p = 0.5f;//spacing
	int nC = (int) (width / p);
	int nR = (int) (length / p);
	float s = (packet->speed1 + packet->speed2)/2.f;
	Vector2f pVel = vF * s;

	for (int w = 0; w <= nC; w++) {
		Vector2f frontEdgePoint = flCornerPt + vR * p * w;
		for (int d = 0; d <= nR; d++) {
			Vector2f depthPoint = frontEdgePoint - vF * p * d;
			splash.particles.push_back(Particle(depthPoint.x(),depthPoint.y(), pVel.x(), pVel.y()));
		}
	}
	//particle generation algorithm for splash particles
	//should initialize a bunch of particles for the splash across the surface of the wave packet with proper velocities
	/*
	for (float y = EPS; y < VIEW_HEIGHT - EPS*2.f; y += H)
		for (float x = VIEW_WIDTH / 4; x <= VIEW_WIDTH / 2; x += H)
			if (splash.particles.size() < DAM_PARTICLES)
			{
				float jitter = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
				splash.particles.push_back(Particle(x + jitter, y, 0.f, 0.f));
			}
	*/
	splashes.insert(std::make_pair(packet, splash));
}

void Particles::RemoveSplash(WAVE_PACKET* packet) {
	splashes.erase(packet);
}

void Particles::ComputeDensityPressure(void) {
	for (auto &splashEntry : splashes) {
		SplashContainer* splash = &splashEntry.second;
		splash->val = 1;
		for (auto &pi : splash->particles)
		{
			pi.rho = 0.f;
			for (auto &pj : splash->particles)
			{
				Vector3f rij = pj.x - pi.x;
				float r2 = rij.squaredNorm();

				if (r2 < HSQ)
				{
					// this computation is symmetric
					pi.rho += MASS*POLY6*pow(HSQ - r2, 3.f);
				}
			}
			pi.p = GAS_CONST*(pi.rho - REST_DENS);
		}
	}
}

void Particles::ComputeForces(void) {
	for (auto &splashEntry : splashes) {
		SplashContainer* splash = &splashEntry.second;
		for (auto &pi : splash->particles)
		{
			Vector3f fpress(0.f, 0.f,0.f);
			Vector3f fvisc(0.f, 0.f,0.f);
			for (auto &pj : splash->particles)
			{
				if (&pi == &pj)
					continue;

				Vector3f rij = pj.x - pi.x;
				float r = rij.norm();

				if (r < H)
				{
					// compute pressure force contribution
					fpress += -rij.normalized()*MASS*(pi.p + pj.p) / (2.f * pj.rho) * SPIKY_GRAD*pow(H - r, 2.f);
					// compute viscosity force contribution
					fvisc += VISC*MASS*(pj.v - pi.v) / pj.rho * VISC_LAP*(H - r);
				}
			}
			Vector3f fgrav = G * pi.rho;
			pi.f = fpress + fvisc + fgrav;
		}
	}
}

void Particles::Integrate() {
	for (auto &splashEntry : splashes) {
		SplashContainer* splash = &splashEntry.second;
		splash->time += DT;
		for (auto &p : splash->particles)
		{
			// forward Euler integration
			p.v += DT*p.f / p.rho;
			p.x += DT*p.v;

			// enforce boundary conditions
			/* //Old boundary check for 2d square
			if (p.x(0) - EPS < 0.0f)
			{
				p.v(0) *= BOUND_DAMPING;
				p.x(0) = EPS;
			}
			if (p.x(0) + EPS > VIEW_WIDTH)
			{
				p.v(0) *= BOUND_DAMPING;
				p.x(0) = VIEW_WIDTH - EPS;
			}
			if (p.x(1) - EPS < 0.0f)
			{
				p.v(1) *= BOUND_DAMPING;
				p.x(1) = EPS;
			}
			if (p.x(1) + EPS > VIEW_HEIGHT)
			{
				p.v(1) *= BOUND_DAMPING;
				p.x(1) = VIEW_HEIGHT - EPS;
			}*/
			//New boundary check for ground geometry
			//Check if the height of the ground at this particles' 2D position is greater than the particle's height above the water.
			
			Vector2f pos2D = p.x2D();
			float groundHeight = GetBoundaryDist(pos2D) * -1.f;
			if (groundHeight > p.x.y()) {
				//TODO figure out calculations for a reset velocity so that splash particles behave properly
				p.v(0) *= -0.5f;
				p.v(2) *= -0.5f;
				if (p.v.y() < -1.f) {
					p.v(1) *= -1.f;
				}
				Vector2f horiz = Vector2f(p.v.x(), p.v.z());
				p.v(1) += horiz.norm();
				//What we want to do is set the particle's position to be at the ground point instead and multiply its velocity by damping
				p.x(1) = groundHeight;
			}
			else if(p.x.y() < 0.05f)
			{
				p.v(1) *= -0.5f;
				p.x(1) = 0.05f;
			}
			
		}
	}
}
