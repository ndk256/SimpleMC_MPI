#include "header.h"

///probably not necessary?
void localize_parameters(Parameters *par, int dims[])
{
//Parameters *p = malloc(sizeof(Parameters));
//*p = *par; ///i forget how pointers work but hopefully this does??

par->gx = par->gx/dims[0];
par->gy = par->gy/dims[1];
par->gz = par->gz/dims[2];
par->n_bins = par->n_bins/((dims[0]+dims[1]+dims[2])/3); ///hopefully this is a good enough approximation

return;
}

void distrib_particle(int d, int mycoords[], Parameters *par, int dim, Bank **sb)
{
int banksz = 0;
Particle *send; send = (Particle*) malloc(((*sb)->n)*sizeof(Particle)); /// attempting to be dynamic 
 
for(int i=1; i<d; i++) // move through the processes in the given dimension from nearest to farthest
{

   for(int i_p=0; i_p<(*sb)->n; i_p++) // iterate through each particle to find if they belong on that process
   {      
/// this could be made more efficient; currently it will do a lot of double-checking
 int pcoords[3] = {(*sb)->p[i_p].x/par->gx, (*sb)->p[i_p].y/par->gy, (*sb)->p[i_p].z/par->gz}; 
//if(i_p==0) printf("\n%.2f,%.2f,%.2f", (*sb)->p[i_p].x,(*sb)->p[i_p].y,(*sb)->p[i_p].z);
//using integer division to get process-level coords
      /// modify ^ if repeated redeclarations isn't good
      if(pcoords[dim]!=mycoords[dim] && pcoords[dim]==i &&i!=mycoords[dim])
      {
         send[banksz] = (*sb)->p[i_p];
       //  if(i_p==0) printf("\nsent as #%d", banksz);
	 banksz++;
         (*sb)->p[i_p].alive = FALSE;         
}
   }
   
MPI_Request mpir = MPI_REQUEST_NULL;
   MPI_Isend(send, banksz, par->type, par->neighb[1+(dim*2)], dim, par->comm, &mpir); /// unsure about the last argument
printf("\n %d sent %d\n\n", par->local_rank, banksz);
   
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
//printf("origin distributed");
}
//printf("%d has n=%d\n", p->local_rank, (*mysb)->n);
MPI_Barrier(p->comm);

///for all the processes handling (x, 0, 0)
if(nprc[0]>1 &&mycoords[0]!=0 && mycoords[1]==0&&mycoords[2]==0) {      
//printf("no");
MPI_Probe(MPI_ANY_SOURCE, 0, p->comm, &status);
      MPI_Get_count(&status, p->type, &msg_size);
      MPI_Recv((*mysb)->p, msg_size, p->type, MPI_ANY_SOURCE, 0, p->comm, MPI_STATUS_IGNORE);
///not sure if ^ should have status or MPI_STATUS_IGNORE
index+=msg_size;
(*mysb)->n += msg_size;

for(int i=0; i<index; i++)
 (*mysb)->p[i].alive=TRUE;

//printf("%d RECVD\n", p->local_rank); 

if(nprc[1]>1)
distrib_particle(nprc[1], mycoords,p, 1, mysb);
}
printf("\nn=%d", (*mysb)->n);
MPI_Barrier(p->comm); //might be unnecc

//for all the processes handling (x, y, 0)
if(nprc[1]>1 && mycoords[2]==0 && mycoords[1]!=0) {
printf("NO");
MPI_Probe(MPI_ANY_SOURCE, 1, p->comm, &status);
MPI_Get_count(&status, p->type, &msg_size);
MPI_Recv((*mysb)->p+index, msg_size, p->type, MPI_ANY_SOURCE, 1, p->comm, &status);

for(int i=index; i<index+msg_size; i++)
 (*mysb)->p[i].alive=TRUE;

if(nprc[2]>1)
distrib_particle(nprc[2], mycoords,p, 2, mysb);

index+=msg_size;
(*mysb)->n += msg_size;
}
//printf("\nn is now %d", (*mysb)->n);
MPI_Barrier(p->comm);

if(nprc[2] >1 && mycoords[2]!=0) {
printf("nou");
MPI_Probe(MPI_ANY_SOURCE, 2, p->comm, &status);
MPI_Get_count(&status, p->type, &msg_size);

//printf("%d received %d\n", p->local_rank, msg_size);
MPI_Recv((*mysb)->p+index, msg_size, p->type, status.MPI_SOURCE, 2, p->comm, &status); 
(*mysb)->n += msg_size;

for(int i=index; i<index+msg_size; i++)
(*mysb)->p[i].alive=TRUE;
} 

//mysb->p = realloc(mysb->p, sizeof(Particle)*index); ///i am flounderings
//printf("\n\nfinally, it's n=%d", (*mysb)->n);
MPI_Barrier(p->comm);

return;}
