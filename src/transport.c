#include "header.h"

// Main logic to move particle
void transport(Parameters *parameters, Geometry *geometry, double local_bounds[], Material *material, Particle tox0[], Particle tox1[], Particle toy0[], Particle toy1[], Particle toz0[], Particle toz1[], int send_indices[], Bank *fission_bank, Tally *tally, Particle *p)
{

  double d_b;
  double d_c;
  double d_e;
  double d;

  while(p->alive){
   //if(parameters->local_rank==1) printf("loop\n");

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
//if(parameters->local_rank==1) printf("d is %.2f", d);
    // Advance particle
    p->x = p->x + d*p->u;
    p->y = p->y + d*p->v;
    p->z = p->z + d*p->w;
    // Case where particle crosses boundary
    if(d_b <= d_c&& d_b<= d_e){
      cross_surface(geometry, p);
//if(parameters->local_rank==1) printf("b");
}

///case of crossing into another process's subdomain
else if(d_e <= d_c && d_e <= d_b)
{
cross_process(local_bounds, p, tox0, tox1, toy0, toy1, toz0, toz1, send_indices);
//if(parameters->local_rank==1) printf("p");
}
else{
    // Case where particle has collision
      collision(material, fission_bank, parameters->nu, p);
//if(parameters->local_rank==1) printf("f");
}
//if(parameters->local_rank==1) printf("\n now to tally");
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
 ///not used:
//  int surfaces[6] = {X0, X1, Y0, Y1, Z0, Z1};
   double p_angles[6] = {p->u, p->u, p->v, p->v, p->w, p->w};
   double p_coords[6] = {p->x, p->x, p->y, p->y, p->z, p->z};
  
   for(i=0; i<6; i++){
     //printf("for dim%d\t", i); 
     if(p_angles[i] == 0){
       // printf("no dim");
        dist = D_INF;
      }
      else{
        dist = (s_coords[i] - p_coords[i])/p_angles[i];
//printf("dist %.6f", dist);       
if(dist <= 0 || dist/p_angles[i]<=0){
         dist = D_INF;
//printf(" which is a no");    
}
      
      if(dist < d){
        d = dist;
//printf("\ncurrent dist is:%.1f\n", d);      
}
    }
 } 
//printf("\tFINAL:%.2f\t", d);   
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
void cross_process(double localbounds[], Particle *p, Particle tox0[], Particle tox1[], Particle toy0[], Particle toy1[], Particle toz0[], Particle toz1[], int indices[])
{
if(p->x==localbounds[0]) {tox0[indices[0]] = *p; indices[0]++;}
 if(p->x==localbounds[1]) {tox1[indices[1]] = *p; indices[1]++;}
 if(p->y==localbounds[2]) {toy0[indices[2]] = *p; indices[2]++;}
 if(p->y==localbounds[3]) {toy1[indices[3]] = *p; indices[3]++;}
 if(p->z==localbounds[4]) {toz0[indices[4]] = *p; indices[4]++;}
 if(p->z==localbounds[5]) {toz1[indices[5]] = *p; indices[5]++;}

p->alive=FALSE; ///this is merely temporary

return;
/*
///locates which process the particle will cross into
///there's probably a better/cleaner way but hopefully this works
if(p->x==localbounds[0]) MPI_Cart_shift(comm, 0, -1, sendto, &sendto);
if(p->x==localbounds[1]) MPI_Cart_shift(comm, 0, 1, sendto, &sendto);
if(p->y==localbounds[2]) MPI_Cart_shift(comm, 1, -1, sendto, &sendto);
if(p->y==localbounds[3]) MPI_Cart_shift(comm, 1, 1, sendto, &sendto);
if(p->z==localbounds[4]) MPI_Cart_shift(comm, 2, -1, sendto, &sendto);
if(p->z==localbounds[5]) MPI_cart_shift(comm, 2, 1, sendto, &sendto);
///passing sendto as the "source process" argument allows for handling of corner cases

if(sendto==MPI_PROC_NULL) {cross_surface(geometry, p); return;} ///reread to check if this case will ever occur/be needed

MPI_Send(p, 1, PARTICLE, sendto, 0, comm);
 ///check which variant to use

return;
*/
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
    if(p->surface_crossed == X0){
      p->u = -p->u;
      p->x = 0.0;
    }
    else if(p->surface_crossed == X1){
      p->u = -p->u;
      p->x = geometry->x;
    }
    else if(p->surface_crossed == Y0){
      p->v = -p->v;
      p->y = 0.0;
    }
    else if(p->surface_crossed == Y1){
      p->v = -p->v;
      p->y = geometry->y;
    }
    else if(p->surface_crossed == Z0){
      p->w = -p->w;
      p->z = 0.0;
    }
    else if(p->surface_crossed == Z1){
      p->w = -p->w;
      p->z = geometry->z;
    }
  }
  
  // Handle periodic boundary conditions
  else if(geometry->bc == PERIODIC){
    if(p->surface_crossed == X0){
      p->x = geometry->x;
    }
    else if(p->surface_crossed == X1){
      p->x = 0;
    }
    else if(p->surface_crossed == Y0){
      p->y = geometry->y;
    }
    else if(p->surface_crossed == Y1){
      p->y = 0;
    }
    else if(p->surface_crossed == Z0){
      p->z = geometry->z;
    }
    else if(p->surface_crossed == Z1){
      p->z = 0;
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
