#include "header.h"

void localize_parameters(Parameters *par, int dims[])
{
par->gx = par->gx/dims[0];
par->gy = par->gy/dims[1];
par->gz = par->gz/dims[2];
par->n_bins = par->n_bins/((dims[0]+dims[1]+dims[2])/3); ///hopefully this is a good enough approximation

return;
}

void distrib_particle(int d, int mycoords[], Parameters *par, int dim, Bank **sb)
{
int banksz = 0;
Particle *send; 
send = (Particle*) malloc(((*sb)->n)*sizeof(Particle)); /// first dynamic line i wrote
 
for(int i=1; i<d; i++) // move through the processes in the given dimension from nearest to farthest
{

   for(int i_p=0; i_p<(*sb)->n; i_p++) // iterate through each particle to find if they belong on that process
   {    
      if((*sb)->p[i_p].alive){
      int pcoords[3] = {(*sb)->p[i_p].x/par->gx, (*sb)->p[i_p].y/par->gy, (*sb)->p[i_p].z/par->gz}; 
      //using integer division to get process-level coords
      /// modify ^ if repeated redeclarations isn't good
      
     if(pcoords[dim]!=mycoords[dim] && pcoords[dim]==i &&i!=mycoords[dim])
      {
         send[banksz] = (*sb)->p[i_p];
	 banksz++;
         (*sb)->p[i_p].alive = FALSE;         
	}
   }
   }
MPI_Request mpir = MPI_REQUEST_NULL;
   MPI_Isend(send, banksz, par->type, par->neighb[1+(dim*2)], dim, par->comm, &mpir); /// unsure about the last argument
   
banksz=0;

}   
free(send);

return;
}

void distribute_sb(int mycoords[], Parameters* p, int nprc[], Bank *sb, Bank **mysb)
{
MPI_Status status;
int msg_size, index=0;

MPI_Barrier(p->comm);

///for the process representing (0,0,0)
if(mycoords[0]==0&&mycoords[1]==0&&mycoords[2]==0)
{(*mysb) = sb;
distrib_particle(nprc[0], mycoords,p, 0, mysb);
}
MPI_Barrier(p->comm);

///for all the processes handling (x, 0, 0)
if(mycoords[1]==0 && mycoords[2]==0) {
     if(mycoords[0]!=0 && nprc[0]>1){
	MPI_Probe(MPI_ANY_SOURCE, 0, p->comm, &status);
      MPI_Get_count(&status, p->type, &msg_size);
      MPI_Recv((*mysb)->p, msg_size, p->type, MPI_ANY_SOURCE, 0, p->comm, MPI_STATUS_IGNORE);
///not sure if ^ should have status or MPI_STATUS_IGNORE
index+=msg_size;
(*mysb)->n += msg_size;

for(int i=0; i<index; i++)
 (*mysb)->p[i].alive=TRUE;
     }
if(nprc[1]>1)
distrib_particle(nprc[1], mycoords,p, 1, mysb);
}
MPI_Barrier(p->comm); //very necc

//for all the processes handling (x, y, 0)
if(mycoords[2]==0) {
   if(nprc[1]>1 && mycoords[1]!=0){
	MPI_Probe(MPI_ANY_SOURCE, 1, p->comm, &status);
MPI_Get_count(&status, p->type, &msg_size);
MPI_Recv((*mysb)->p+index, msg_size, p->type, MPI_ANY_SOURCE, 1, p->comm, &status);

for(int i=index; i<index+msg_size; i++)
 (*mysb)->p[i].alive=TRUE;

index+=msg_size;
(*mysb)->n += msg_size;}
	   	   
if(nprc[2]>1)
distrib_particle(nprc[2], mycoords,p, 2, mysb);
}
MPI_Barrier(p->comm);

if(nprc[2] > 1 && mycoords[2]!=0) {
MPI_Probe(MPI_ANY_SOURCE, 2, p->comm, &status);
MPI_Get_count(&status, p->type, &msg_size);

MPI_Recv((*mysb)->p+index, msg_size, p->type, status.MPI_SOURCE, 2, p->comm, &status); 
(*mysb)->n += msg_size;

for(int i=index; i<index+msg_size; i++)
(*mysb)->p[i].alive=TRUE;
} 

MPI_Barrier(p->comm);

return;}
