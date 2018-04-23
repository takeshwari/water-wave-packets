#include "SPHParticles.h"
//This class keeps track of all SPH-simulated water particles in a scene and simulates them each update step
void Particles::Update(Packets* packets) {
	ComputeDensityPressure();
	ComputeForces();
	Integrate(packets);
}

//given a wave packet, generate a splash container for it and add it to the simulation
void Particles::AddSplash(WAVE_PACKET* packet) {
	SplashContainer splash = SplashContainer();
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
		SplashContainer splash = splashEntry.second;
		for (auto &pi : splash.particles)
		{
			pi.rho = 0.f;
			for (auto &pj : splash.particles)
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
		SplashContainer splash = splashEntry.second;
		for (auto &pi : splash.particles)
		{
			Vector3f fpress(0.f, 0.f);
			Vector3f fvisc(0.f, 0.f);
			for (auto &pj : splash.particles)
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

void Particles::Integrate(Packets* packets) {
	for (auto &splashEntry : splashes) {
		SplashContainer splash = splashEntry.second;
		for (auto &p : splash.particles)
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
			float groundHeight = packets->GetBoundaryDist(p.x2D()) * -1.f;
			if (groundHeight < p.x.z()) {
				//TODO figure out calculations for a reset velocity so that splash particles behave properly
				p.v *= 0.5f;
				//What we want to do is set the particle's position to be at the ground point instead and multiply its velocity by damping
				p.x(2) = groundHeight;
			}
		}
	}
}
