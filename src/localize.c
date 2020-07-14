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
void distrib_particle(int prcoords[], int mycoords[], Parameters *par, int dim, Bank **sb)
{
int d = prcoords[dim]; //number of processes on the given axis
int banksz = 0;
Particle *send; send = (Particle*) malloc((sb->n)*sizeof(Particle)); /// attempting to be dynamic
int destrank;
   
for(int i=1; i<d; i++) // move through the processes in the given dimension from nearest to farthest
{
   for(int i_p=0; i_p<sb->n; i_p++) // iterate through each particle to find if they belong on that process
   {
      int pcoords[3] = {(*sb)->p[i_p].x/par->gx, (*sb)->p[i_p].y/par->gy, (*sb)->    p[i_p].z/par->gz};         
     //using integer division to get process-level coords
        /// (modify ^ if repeated redeclarations isn't good)
        if(pcoords[dim] != mycoords[dim] && pcoords[dim] == i)
      /// this could be made more efficient; currently it will do a lot of double-checking
      {
         send[banksz] = sb->p[i_p];
         banksz++;
         sb->p[i_p]->alive = FALSE; ///double-check if correct, but basically kill the version that's not on the right processes
       //experimental alternative:
         /*
         sb->p[i_p] = sb->p[sb->n];
         sb->n--;
         */
      }
   }
   MPI_Cart_shift(p->comm, dim, i, &(p->local_rank), &destrank);
    MPI_Request mpir = MPI_REQUEST_NULL;
   MPI_Isend(send, banksz, p->type, destrank, dim, p->comm, &mpir); /// unsure about the last argument
   banksz=0;
}   

return;
}

void distribute_sb(int mycoords[], Parameters *p, int prcoords[], Bank *sb, Bank **mysb)
{
MPI_Status status;
int msg_size, index=0;

MPI_Barrier(p->comm);

//for the process representing (0,0,0)
if(mycoords[0]==0 && mycoords[1]==0 && mycoords[2] == 0)
{
   (*mysb) = sb;
distrib_particle(prcoords, mycoords, p, 0, mysb);
}

MPI_Barrier(p->comm);

//for the processes handling (x,0,0)
if(mycoords[1]==0 && mycoords[2]==0){
if(mycoords[0]!=0)
{MPI_Probe(MPI_ANY_SOURCE, 0, parameters->comm, &status);
MPI_Get_count(&status, p->type, &msg_size);
MPI_Recv(mysb->p, msg_size, p->type, MPI_ANY_SOURCE, 0, p->comm, MPI_STATUS_IGNORE);
index+=msg_size;
(*mysb)->n += msg_size;}

distrib_particle(prcoords, mycoords, p, 1, mysb);
}

MPI_Barrier(p->comm);

//for the processes handling (x,y,0)
if(mycoords[2]==0 && mycoords[0]!=0) {
MPI_Probe(MPI_ANY_SOURCE, 1, p->comm, &status);
MPI_Get_count(&status, p->type, &msg_size);
MPI_Recv(mysb->p+index, msg_size, p->type, MPI_ANY_SOURCE, 1, p->comm, &status);
distrib_particle(prcoords, mycoords, p, 2, mysb);
index+=msg_size;
   (*mysb)->n += msg_size;
}

MPI_Barrier(p->comm);

if(mycoords[2]!=0)
{MPI_Probe(MPI_ANY_SOURCE, 2, p->comm, &status);
MPI_Get_count(&status, p->type, &msg_size);
MPI_Recv((*mysb)->p+index, msg_size, p->type, status.MPI_SOURCE, 2, p->comm, &status);
(*mysb)->n += msg_size;}

(*mysb)->p = realloc((*mysb)->p, sizeof(Particle)*index); //no idea if this is right
return;
}
