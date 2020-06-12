#include <header.h>

Parameters *localize_parameters(Parameters *par,int dims[])
   {
   Parameters *p = malloc(sizeof(Parameters));
   *p = *par; //i forget how pointers work but hopefully this does
   
   p->gx = par->gx/dims[0];
   p->gy = par->gy/dims[1];
   p->gz = par->gz/dims[2];
   return p;
}

//Sends particles to processes along the axis indicated by dim
void distrib_particle(int prcoords[], int dim, MPI_Comm comm, Bank *sb, int myrank)
{
int d = prcoords[dim]; //number of processes on the given axis
Particle sendbank[sb->n][d];
int banksz[d];

for(int i_p=0; i_p<sb->n; i_p++)
{
double pcoords[3] = {sb->p[i_p].x, sb->p[i_p].y, sb->p[i_p].z};
int i= pcoords[dim]/d; //due to integer division, this should give the coordinates to send to
sendbank[banksz[i]][i] = sb->p[i_p];
banksz[i]++;}

int destrank;

for(int j=0; j<prcoords[dim]; j++)
{
MPI_Cart_shift(comm, dim, j, &myrank, &destrank);

//MPI_Isend(sendbank[][j], banksz[j], PARTICLE, destrank, dim, comm, MPI_REQUEST_NULL);
//the above line is not correct yet
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
if(mycoords[0]==0 && mycoords[2]==0){
if(mycoords[1]!=0)
{MPI_Probe(MPI_ANY_SOURCE, 0, comm, &status);
MPI_Get_count(&status, PARTICLE, &msg_size);
MPI_Recv(mysb, msg_size, PARTICLE, MPI_ANY_SOURCE, 0, comm, MPI_STATUS_IGNORE);
index+=msg_size;}

distrib_particle(prcoords, 1, comm, mysb, myrank);
}

MPI_Barrier(comm);

//for the processes handling (x,y,0)
if(mycoords[2]==0){
MPI_Probe(MPI_ANY_SOURCE, 1, comm, &status);
MPI_Get_count(&status, PARTICLE, &msg_size);
MPI_Recv(mysb+index, msg_size, PARTICLE, MPI_ANY_SOURCE, 1, comm, &status);
distrib_particle(prcoords, 2, comm, mysb, myrank);
index+=msg_size;
}

MPI_Barrier(comm);

if(mycoords[2]!=0)
{MPI_Probe(MPI_ANY_SOURCE, 2, comm, &status);
MPI_Get_count(&status, PARTICLE, &msg_size);
MPI_Recv(mysb+index, msg_size, PARTICLE, status.MPI_SOURCE, 2, comm, &status);
}

mysb->p = realloc(mysb->p, sizeof(Particle)*index); //no idea if this is right
return;
}
