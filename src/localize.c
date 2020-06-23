#include "header.h"

Parameters *localize_parameters(Parameters *par,int dims[])
   {
   Parameters *p = malloc(sizeof(Parameters));
   *p = *par; //i forget how pointers work but hopefully this does
   
   p->gx = par->gx/dims[0];
   p->gy = par->gy/dims[1];
   p->gz = par->gz/dims[2];
   p->n_bins = par->n_bins/((dims[0]+dims[1]+dims[2])/3); ///hopefully this is a good enough approximation
   return p;
}

//Sends particles to processes along the axis indicated by dim
void distrib_particle(int prcoords[], int dim, MPI_Comm comm, Bank *sb, int myrank)
{
int d = prcoords[dim]; //number of processes on the given axis
int banksz = 0;
Particle *send; send = (Particle*) malloc((sb->n)*sizeof(Particle)); /// attempting to be dynamic
int destrank;
   
for(int i=1; i<d; i++) // move through the processes in the given dimension from nearest to farthest
{
   for(int i_p=0; i_p<sb->n; i_p++) // iterate through each particle to find if they belong on that process
   {
      /// this could be made more efficient; currently it will do a lot of double-checking
      double pcoords[3] = {sb->p[i_p].x, sb->p[i_p].y, sb->p[i_p].z}; // make p's coordinates into a generic/numeric reference
      /// modify ^ if repeated redeclarations isn't good
      if(pcoords[dim]/d == i)
      {
         send[banksz[i]] = sb->p[i_p];
         banksz++;
         sb->p[i_p]->alive = 0; ///double-check if correct, but basically kill the version that's not on the right processes
       //experimental alternative:
         /*
         sb->p[i_p] = sb->p[sb->n];
         sb->n--;
         */
      }
   }
   MPI_Cart_shift(comm, dim, i, &myrank, &destrank);
   MPI_Isend(send, banksz, PARTICLE, destrank, dim, comm, MPI_REQUEST_NULL); /// unsure about the last argument
   banksz=0;
}   

return;
}

void distribute_sb(int mycoords[], int prcoords[], int myrank, MPI_Comm comm, Bank *sb, Bank *mysb)
{
MPI_Status status;
int msg_size, index=0;

MPI_Barrier(comm);

//for the process representing (0,0,0)
if(mycoords[0]==0 && mycoords[1]==0 && mycoords[2] == 0)
{
distrib_particle(prcoords, 0, comm, sb, myrank);
}

MPI_Barrier(comm);

//for the processes handling (x,0,0)
if(mycoords[1]==0 && mycoords[2]==0){
if(mycoords[0]!=0)
{MPI_Probe(MPI_ANY_SOURCE, 0, comm, &status);
MPI_Get_count(&status, PARTICLE, &msg_size);
MPI_Recv(mysb->p, msg_size, PARTICLE, MPI_ANY_SOURCE, 0, comm, MPI_STATUS_IGNORE);
index+=msg_size;
mysb->n += msg_size;}

distrib_particle(prcoords, 1, comm, mysb, myrank);
}

MPI_Barrier(comm);

//for the processes handling (x,y,0)
if(mycoords[2]==0){
MPI_Probe(MPI_ANY_SOURCE, 1, comm, &status);
MPI_Get_count(&status, PARTICLE, &msg_size);
MPI_Recv(mysb->p+index, msg_size, PARTICLE, MPI_ANY_SOURCE, 1, comm, &status);
distrib_particle(prcoords, 2, comm, mysb, myrank);
index+=msg_size;
   mysb->n += msg_size;
}

MPI_Barrier(comm);

if(mycoords[2]!=0)
{MPI_Probe(MPI_ANY_SOURCE, 2, comm, &status);
MPI_Get_count(&status, PARTICLE, &msg_size);
MPI_Recv(mysb->p+index, msg_size, PARTICLE, status.MPI_SOURCE, 2, comm, &status);
mysb->n += msg_size;}

mysb->p = realloc(mysb->p, sizeof(Particle)*index); //no idea if this is right
return;
}
