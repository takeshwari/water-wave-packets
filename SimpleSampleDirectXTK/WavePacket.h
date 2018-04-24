#include <Eigen/Dense>
using namespace Eigen;
#pragma once
struct WAVE_PACKET
{
	// positions, directions, speed of the tracked vertices
	Vector2f	pos1, pos2, pos3;				// 2D position
	Vector2f	dir1, dir2, dir3;				// current movement direction
	float		speed1, speed2, speed3;		// speed of the particle
	Vector2f	pOld1, pOld2, pOld3;			// position in last timestep (needed to handle bouncing)
	Vector2f	dOld1, dOld2, dOld3;			// direction in last timestep (needed to handle bouncing)
	float		sOld1, sOld2, sOld3;			// speed in last timestep (needed to handle bouncing)
	Vector2f	midPos;						// middle position (tracked each timestep, used for rendering)
	Vector2f	travelDir;					// travel direction (tracked each timestep, used for rendering)
	float		bending;					// point used for circular arc bending of the wave function inside envelope

											// bouncing and sliding
	bool		bounced1, bounced2, bounced3;	// indicates if this vertex bounced in this timestep
	bool		sliding3;					// indicates if the 3rd vertex is "sliding" (used for diffraction)
	bool		use3rd;						// indicates if the third vertex is present (it marks a (potential) sliding point)
											// wave function related
	float		phase;						// phase of the representative wave inside the envelope, phase speed vs. group speed
	float		phOld;						// old phase
	float		E;							// wave energy flux for this packet (determines amplitude)
	float		envelope;					// envelope size for this packet
	float		k, w0;						// w0 = angular frequency, k = current wavenumber
	float		k_L, w0_L, k_H, w0_H;			// w0 = angular frequency, k = current wavenumber,  L/H are for lower/upper boundary
	float		d_L, d_H;					// d = travel distance to reference wave (gets accumulated over time),  L/H are for lower/upper boundary
	float		ampOld;						// amplitude from last timestep, will be smoothly adjusted in each timestep to meet current desired amplitude
	float		dAmp;						// amplitude change in each timestep (depends on desired waveheight so all waves (dis)appear with same speed)
											// serial deletion step variable
	bool		toDelete;					// used internally for parallel deletion criterion computation
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};