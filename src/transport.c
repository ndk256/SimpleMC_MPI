#include "header.h"

// Main logic to move particle
void transport(Parameters *parameters, Geometry *geometry, double local_bounds[], Material *material, Buffer *sendbuf, Bank *fission_bank, Tally *tally, Particle *p)
{

  double d_b;
  double d_c;
  double d_e;
  double d;

  while(p->alive){
    
    ///Find distance to boundary of process section
    d_e = dist_to_edge(p, local_bounds); ///    

    // Find distance to boundary
    d_b = distance_to_boundary(geometry, p);

    // Find distance to collision
    d_c = distance_to_collision(material);

    // Take smallest of the distances
   /// d = d_b < d_c ? d_b : d_c;
if(d_b <= d_e && d_b <= d_c) d = d_b;
else if (d_e <= d_b && d_e <= d_c) d = d_e;
else d=d_c;

    // Advance particle
    p->x = p->x + d*p->u;
    p->y = p->y + d*p->v;
    p->z = p->z + d*p->w;
    // Case where particle crosses boundary
    if(d_b <= d_c&& d_b<= d_e){
      cross_surface(geometry, p);
}

///case of crossing into another process's subdomain
else if(d_e <= d_c && d_e <= d_b)
{
cross_process(local_bounds, p, sendbuf);
}
else{
    // Case where particle has collision
      collision(material, fission_bank, parameters->nu, p);
     }
      // Score tallies
      if(tally->tallies_on == TRUE){
        score_tally(parameters, material, tally, p);
      } 
  }
  return;
}

///Returns distance to nearest boundary of the process subdivision
double dist_to_edge(Particle *p, double s_coords[])
{
///mostly just a copy-paste of the distance_to_boundary function.
int i; double dist;
   double d = D_INF;
   double p_angles[6] = {p->u, p->u, p->v, p->v, p->w, p->w};
   double p_coords[6] = {p->x, p->x, p->y, p->y, p->z, p->z};
  
   for(i=0; i<6; i++){
     if(p_angles[i] == 0){
        dist = D_INF;
      }
      else{
        dist = (s_coords[i] - p_coords[i])/p_angles[i];
      if(dist < 0 || dist/p_angles[i]<=0){
         dist = D_INF;
}
      
      if(dist < d){
        d = dist;
      }
    }
 } 
 return d;
}

// Returns the distance to the nearest boundary for a particle traveling in a
// certain direction
double distance_to_boundary(Geometry *geometry, Particle *p)
{
  int i;
  double dist;
  double d = D_INF;
  int    surfaces[6] = {X0, X1, Y0, Y1, Z0, Z1};
  double p_angles[6] = {p->u, p->u, p->v, p->v, p->w, p->w};
  double p_coords[6] = {p->x, p->x, p->y, p->y, p->z, p->z};
  double s_coords[6] = {0, geometry->x, 0, geometry->y, 0, geometry->z};
  
  for(i=0; i<6; i++){
    if(p_angles[i] == 0){
      dist = D_INF;
    }
    else{
      dist = (s_coords[i] - p_coords[i])/p_angles[i];
      if(dist <= 0){
        dist = D_INF;
      }
    }
    if(dist < d){
      d = dist;
      p->surface_crossed = surfaces[i];
    }
  }

  return d;
}

// Returns the distance to the next collision for a particle
double distance_to_collision(Material *material)
{
  double d;

  if(material->xs_t == 0){
    d = D_INF;
  }
  else{
    d = -log(rn())/material->xs_t;
  }

  return d;
}

///handles particle crossing into another process's subdomain
///simply sends it to a buffer which will be passed out to the appropriate process later on
void cross_process(double localbounds[], Particle *p, Buffer *sendbuf)
 {
 if(p->x<=localbounds[0]) {
  if(sendbuf->n_banked[0]>=sendbuf->banksz[0])
   {sendbuf->tox0=realloc(sendbuf->tox0, sizeof(Particle)*2*sendbuf->banksz[0    ]);
    sendbuf->banksz[0] *= 2;}
  sendbuf->tox0[sendbuf->n_banked[0]] = *p;
  sendbuf->n_banked[0]++;}
 else if(p->x>=localbounds[1]) {
  if(sendbuf->n_banked[1]>=sendbuf->banksz[1])
   {sendbuf->tox1=realloc(sendbuf->tox1, sizeof(Particle)*2*sendbuf->banksz[1]);
    sendbuf->banksz[1] *= 2;}
  sendbuf->tox1[sendbuf->n_banked[1]] = *p;
  sendbuf->n_banked[1]++;}
 else if(p->y<=localbounds[2]) {
  if(sendbuf->n_banked[2]>=sendbuf->banksz[2])
   {sendbuf->toy0=realloc(sendbuf->toy0, sizeof(Particle)*2*sendbuf->banksz[2]);
    sendbuf->banksz[2] *= 2;}
  sendbuf->toy0[sendbuf->n_banked[2]] = *p;
  sendbuf->n_banked[2]++;}
 else if(p->y>=localbounds[3]) {
  if(sendbuf->n_banked[3]>=sendbuf->banksz[3])
   {sendbuf->toy1=realloc(sendbuf->toy1, sizeof(Particle)*2*sendbuf->banksz[3]);
    sendbuf->banksz[3] *= 2;}
  sendbuf->toy1[sendbuf->n_banked[3]] = *p;
  sendbuf->n_banked[3]++;}
 else if(p->z<=localbounds[4]) {
  if(sendbuf->n_banked[4]>=sendbuf->banksz[4])
   {sendbuf->toz0=realloc(sendbuf->toz0, sizeof(Particle)*2*sendbuf->banksz[4]);
    sendbuf->banksz[4] *= 2;}
  sendbuf->toz0[sendbuf->n_banked[4]] = *p;
  sendbuf->n_banked[4]++;}
 else if(p->z>=localbounds[5]) {
  if(sendbuf->n_banked[5]>=sendbuf->banksz[5])
   {sendbuf->toz1=realloc(sendbuf->toz1, sizeof(Particle)*2*sendbuf->banksz[5]);
    sendbuf->banksz[5] *= 2;}
  sendbuf->toz1[sendbuf->n_banked[5]] = *p;
  sendbuf->n_banked[5]++;}
 
 p->alive=FALSE; ///this is merely temporary
 
 return;
 }

// Handles a particle crossing a surface in the geometry
void cross_surface(Geometry *geometry, Particle *p)
{
  // Handle vacuum boundary conditions (particle leaks out)
  if(geometry->bc == VACUUM){
    p->alive = FALSE;
  }

  // Handle reflective boundary conditions
  else if(geometry->bc == REFLECT){
    if(p->x < 0.0){
      p->u = -p->u;
      p->x = 0.0;
    }
    else if(p->x > geometry->x){
      p->u = -p->u;
      p->x = geometry->x;
    }
    else if(p->y < 0.0){
      p->v = -p->v;
      p->y = 0.0;
    }
    else if(p->y > geometry->y){
      p->v = -p->v;
      p->y = geometry->y;
    }
    else if(p->z < 0.0){
      p->w = -p->w;
      p->z = 0.0;
    }
    else if(p->z > geometry->z){
      p->w = -p->w;
      p->z = geometry->z;
    }
  }
  
  // Handle periodic boundary conditions
  else if(geometry->bc == PERIODIC){
    if(p->surface_crossed == X0){
      p->x = geometry->x;
      if(sendbuf->n_banked[0]>=sendbuf->banksz[0])
   {sendbuf->tox0=realloc(sendbuf->tox0, sizeof(Particle)*2*sendbuf->banksz[0]);
    sendbuf->banksz[0] *= 2;}
  sendbuf->tox0[sendbuf->n_banked[0]] = *p;
  sendbuf->n_banked[0]++;
    }
    else if(p->surface_crossed == X1){
      p->x = 0;
      if(sendbuf->n_banked[1]>=sendbuf->banksz[1])
   {sendbuf->tox1=realloc(sendbuf->tox1, sizeof(Particle)*2*sendbuf->banksz[1]);
    sendbuf->banksz[1] *= 2;}
  sendbuf->tox1[sendbuf->n_banked[1]] = *p;
  sendbuf->n_banked[1]++;
    }
    else if(p->surface_crossed == Y0){
      p->y = geometry->y;
      if(sendbuf->n_banked[2]>=sendbuf->banksz[2])
   {sendbuf->toy0=realloc(sendbuf->toy0, sizeof(Particle)*2*sendbuf->banksz[2]);
    sendbuf->banksz[2] *= 2;}
  sendbuf->toy0[sendbuf->n_banked[2]] = *p;
  sendbuf->n_banked[2]++;
    }
    else if(p->surface_crossed == Y1){
      p->y = 0;
      if(sendbuf->n_banked[3]>=sendbuf->banksz[3])
   {sendbuf->toy1=realloc(sendbuf->toy1, sizeof(Particle)*2*sendbuf->banksz[3]);
    sendbuf->banksz[3] *= 2;}
  sendbuf->toy1[sendbuf->n_banked[3]] = *p;
  sendbuf->n_banked[3]++;
    }
    else if(p->surface_crossed == Z0){
      p->z = geometry->z;
      if(sendbuf->n_banked[4]>=sendbuf->banksz[4])
   {sendbuf->toz0=realloc(sendbuf->toz0, sizeof(Particle)*2*sendbuf->banksz[4]);
    sendbuf->banksz[4] *= 2;}
  sendbuf->toz0[sendbuf->n_banked[4]] = *p;
  sendbuf->n_banked[4]++;
    }
    else if(p->surface_crossed == Z1){
      p->z = 0;
      if(sendbuf->n_banked[5]>=sendbuf->banksz[5])
   {sendbuf->toz1=realloc(sendbuf->toz1, sizeof(Particle)*2*sendbuf->banksz[5]);
    sendbuf->banksz[5] *= 2;}
  sendbuf->toz1[sendbuf->n_banked[5]] = *p;
  sendbuf->n_banked[5]++;
    }
  }

  return;
}

void collision(Material *material, Bank *fission_bank, double nu, Particle *p)
{
  int n;
  int i = 0;
  double prob = 0.0;
  double cutoff;
  Nuclide nuc = {0, 0, 0, 0, 0};

  // Cutoff for sampling nuclide
  cutoff = rn()*material->xs_t;

  // Sample which nuclide particle has collision with
  while(prob < cutoff){
    nuc = material->nuclides[i];
    prob += nuc.atom_density*nuc.xs_t;
    i++;
  }

  // Cutoff for sampling reaction
  cutoff = rn()*nuc.xs_t;

  // Sample fission
  if(nuc.xs_f > cutoff){

    // Sample number of fission neutrons produced
    if(rn() > nu - (int)nu){
      n = nu;
    }
    else{
      n = nu + 1;
    }

    // Sample n new particles from the source distribution but at the current
    // particle's location
    if(fission_bank->n+n >= fission_bank->sz){
      fission_bank->resize(fission_bank);
    }
    for(i=0; i<n; i++){
      sample_fission_particle(&(fission_bank->p[fission_bank->n]), p);
      fission_bank->n++;
    }
    p->alive = FALSE;
  }

  // Sample absorption (disappearance)
  else if(nuc.xs_a > cutoff){
    p->alive = FALSE;
  }

  // Sample scattering
  else{
    p->mu = rn()*2 - 1;
    p->phi = rn()*2*PI;
    p->u = p->mu;
    p->v = sqrt(1 - p->mu*p->mu) * cos(p->phi);
    p->w = sqrt(1 - p->mu*p->mu) * sin(p->phi);
  }

  return;
}
