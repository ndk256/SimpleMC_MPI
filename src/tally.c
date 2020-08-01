#include "header.h"

// Simple flux tally
void score_tally(Parameters *parameters, Material *material, Tally *tally, Particle *p)
{
  int ix, iy, iz;
  double vol,px=p->x,py=p->y,pz=p->z;
  
  while(px>=tally->dx*tally->n) px -= tally->n*tally->dx;
  while(py>=tally->dy*tally->n) py -= tally->dy*tally->n;
  while(pz>=tally->dz*tally->n) pz -= tally->dz*tally->n;

  // Volume
  vol = tally->dx * tally->dy * tally->dz;

  // Find the indices of the grid box of the particle
  ix = px/tally->dx;
  iy = py/tally->dy;
  iz = pz/tally->dz;

  // Scalar flux
  tally->flux[ix + tally->n*iy + tally->n*tally->n*iz] += 1./(vol * material->xs_t * parameters->n_particles);

  return;
}

void reset_tally(Tally *tally)
{
  memset(tally->flux, 0, tally->n*tally->n*tally->n*sizeof(double));

  return;
}
