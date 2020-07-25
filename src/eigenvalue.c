#include "header.h"

void run_eigenvalue(double localbounds[6], Parameters *parameters, Geometry *geometry, Material *material, Bank *source_bank, Bank *fission_bank, Tally *tally, double *keff)
{
  int i_b; // index over batches
  int i_a = -1; // index over active batches
  int i_g; // index over generations
  unsigned long i_p; // index over particles
  double keff_gen = 1; // keff of generation
  double keff_batch=0; // keff of batch
  double keff_mean=0; // keff mean over active batches
  double keff_std=0; // keff standard deviation over active batches

  // Loop over batches
  for(i_b=0; i_b<parameters->n_batches; i_b++){

    keff_batch = 0;

    // Turn on tallying and increment index in active batches
    if(i_b >= parameters->n_batches - parameters->n_active){
      i_a++;
      if(parameters->tally == TRUE){
        tally->tallies_on = TRUE;
      }
    }

    // Loop over generations
    for(i_g=0; i_g<parameters->n_generations; i_g++){

      // Set RNG stream for tracking
      set_stream(STREAM_TRACK);

      //set up buffer arrays for cross-process transport
      Particle tox0[parameters->n_particles]; ///possibly overdoing it on size but idk
Particle tox1[parameters->n_particles];
Particle toy0[parameters->n_particles];
Particle toy1[parameters->n_particles];
Particle toz0[parameters->n_particles];
Particle toz1[parameters->n_particles];
      int send_indices[6] = { }; ///initializes empty array
      // Loop over particles
      for(i_p=0; i_p<source_bank->n; i_p++){ //
	// Set seed for particle i_p by skipping ahead in the random number
	// sequence stride*(total particles simulated) numbers from the initial
	// seed. This allows for reproducibility of the particle history.
        rn_skip((i_b*parameters->n_generations + i_g)*parameters->n_particles + parameters->local_rank*i_p);
//^this will give some repetition still, but not as much as without multiplying by the local rank
	      
        // Transport the next particle
        transport(parameters, geometry, localbounds, material, tox0, tox1, toy0, toy1, toz0, toz1, send_indices, fission_bank, tally, &(source_bank->p[i_p])); 
}
int tosend=0, stillsend=0;
	
MPI_Barrier(parameters->comm);
	    
do{
stillsend=0;

if(send_indices[0]>0 || send_indices[1]>0 || send_indices[2]>0 || send_indices[3]>0 || send_indices[4]>0 || send_indices[5]>0) 
{tosend=1;}
else tosend=0;

MPI_Barrier(parameters->comm);
MPI_Allreduce(&tosend, &stillsend, 1, MPI_INT, MPI_SUM, parameters->comm);

if(stillsend==0) break;

///send particles to the instance of source_bank on the appropriate process
sendrecv_particles(parameters, source_bank, tox0, tox1, toy0, toy1, toz0, toz1, send_indices, localbounds);

  ///transport those particles
for(i_p =0; i_p<source_bank->n; i_p++){ 
	rn_skip((i_b*parameters->n_generations+i_g)*parameters->n_particles+i_p*parameters->local_rank*parameters->local_rank);

	///transport particle(s)
	transport(parameters, geometry, localbounds, material, tox0, tox1, toy0, toy1, toz0, toz1, send_indices, fission_bank, tally, &(source_bank->p[i_p]));
  }
}while(stillsend>0);

      // Switch RNG stream off tracking
      set_stream(STREAM_OTHER);
      rn_skip(i_b*parameters->n_generations + i_g);
////^ change?

MPI_Barrier(parameters->comm);

// Calculate generation k_effective and accumulate batch k_effective
      int total_fish;
//printf("%d's fissions:%d\n", parameters->local_rank, fission_bank->n);
      MPI_Reduce(&fission_bank->n, &total_fish, 1, MPI_INT, MPI_SUM, 0, parameters->comm);
      keff_gen = (double) total_fish / parameters->n_particles; ///ok to use parameters for that value?
      keff_batch += keff_gen;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank(source_bank, fission_bank);    
}

if(parameters->local_rank==0){   
 // Calculate k effective
    keff_batch /= parameters->n_generations;
    if(i_a >= 0){
      keff[i_a] = keff_batch;
    }
    calculate_keff(keff, &keff_mean, &keff_std, i_a+1);
}
    // Tallies for this realization 
   if(tally->tallies_on == TRUE){
      if(parameters->write_tally == TRUE)
       {
        Tally *overall_tally = init_tally(parameters);///
	MPI_Reduce(tally->flux, overall_tally->flux, tally->n*tally->n*tally->n, MPI_DOUBLE, MPI_SUM, 0, parameters->comm); ///check what the size of tally->flux[] is
 if(parameters->local_rank==0) write_tally(overall_tally, parameters->tally_file);///
     
      reset_tally(overall_tally);
    }}
	if(parameters->local_rank==0){
// Status text
    print_status(i_a, i_b, keff_batch, keff_mean, keff_std);
	}
}
reset_tally(tally);
   
  // Write out keff
  if(parameters->local_rank==0 && parameters->write_keff == TRUE){
    write_keff(keff, parameters->n_active, parameters->keff_file);
  }
  return;
}

void sendrecv_particles(Parameters *p, Bank *bank, Particle tox0[], Particle tox1[], Particle toy0[], Particle toy1[], Particle toz0[], Particle toz1[], int send_indices[6], double mybounds[6])
{
/// To do: re-allocate bank->p memory? check if necessary
	
MPI_Status status; int n_to_recv, recv_count=0;
MPI_Request reqs[6];

///sending along x direction///
MPI_Isend(tox0, send_indices[0]+1, p->type, p->neighb[0], 0, p->comm,&reqs[0]);
MPI_Isend(tox1, send_indices[1]+1, p->type, p->neighb[1], 1, p->comm, &reqs[1]);

send_indices[0]=0; send_indices[1]=0; ///resetting

///get any particles that have crossed in from other processes via x
MPI_Probe(p->neighb[0], 1, p->comm, &status); //this is blocking--is that desired?
MPI_Get_count(&status, p->type, &n_to_recv);

MPI_Recv(bank->p, n_to_recv, p->type, p->neighb[0], 1, p->comm, MPI_STATUS_IGNORE);
recv_count+=n_to_recv;

MPI_Probe(p->neighb[1], 0, p->comm, &status);
MPI_Get_count(&status, p->type, &n_to_recv);

MPI_Recv(bank->p+recv_count, n_to_recv, p->type, p->neighb[1], 0, p->comm, MPI_STATUS_IGNORE);
recv_count+=n_to_recv;

int b_ind=recv_count;
 
// sorting through the received particles to see which have/n't reached their final destination
for(int i=0; i<recv_count; i++)
{if(bank->p[i].y < mybounds[3] && bank->p[i].y >= mybounds[2] && bank->p[i].z < mybounds[5] && bank->p[i].z >= mybounds[4])
{bank->p[i].alive = TRUE;
} ///if it's in the right place it can be deemed alive again
else if(bank->p[i].y < mybounds[2]) {
toy0[send_indices[2]] = bank->p[i]; 
send_indices[2]++;
}
else if(bank->p[i].y >= mybounds[3]) {
toy1[send_indices[3]] = bank->p[i];
 send_indices[3]++;}
else if(bank->p[i].z < mybounds[4]) {
toz0[send_indices[4]] = bank->p[i]; 
send_indices[4]++;}
else if(bank->p[i].z >= mybounds[5]) {
toz1[send_indices[5]] = bank->p[i]; 
send_indices[5]++;}
}

MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);

///here is where y sending should begin///
MPI_Isend(toy0, send_indices[2]+1, p->type, p->neighb[2], 2, p->comm, &reqs[2]);
MPI_Isend(toy1, send_indices[3]+1, p->type, p->neighb[3], 3, p->comm,&reqs[3]);

send_indices[2]=0; send_indices[3]=0; ///resetting

MPI_Probe(p->neighb[2], 3, p->comm, &status);
MPI_Get_count(&status, p->type, &n_to_recv);

MPI_Recv(bank->p+recv_count, n_to_recv, p->type, p->neighb[2], 3, p->comm, MPI_STATUS_IGNORE);
recv_count+= n_to_recv;

MPI_Probe(p->neighb[3], 2, p->comm, &status);
  MPI_Get_count(&status, p->type, &n_to_recv);

  MPI_Recv(bank->p+recv_count, n_to_recv, p->type, p->neighb[3], 2, p->comm,  MPI_STATUS_IGNORE);
  recv_count+= n_to_recv;

for(int i=b_ind; i<recv_count; i++) 
  {if(bank->p[i].z < mybounds[5] && bank->p[i].z >= mybounds[4])
  {bank->p[i].alive = TRUE; } ///if it's in the right place it can be deemed alive agai    n
  else if(bank->p[i].z < mybounds[4]) {
toz0[send_indices[4]] = bank->p[i];
 send_indices[4]++;}
  else if(bank->p[i].z >= mybounds[5]) {
toz1[send_indices[5]] = bank->p[i]; 
send_indices[5]++;}
  }

MPI_Wait(&reqs[2], MPI_STATUS_IGNORE);
MPI_Wait(&reqs[3], MPI_STATUS_IGNORE);

///z sending////
MPI_Isend(toz0, send_indices[4], p->type, p->neighb[4], 4, p->comm, &reqs[4]);
MPI_Isend(toz1, send_indices[5], p->type, p->neighb[5], 5, p->comm, &reqs[5]);

send_indices[4]=0; send_indices[5]=0;

////z receiving
MPI_Probe(p->neighb[4], 5, p->comm, &status);
  MPI_Get_count(&status, p->type, &n_to_recv);
  
  MPI_Recv(bank->p+recv_count, n_to_recv, p->type, p->neighb[4], 5, p->comm,     MPI_STATUS_IGNORE);
  recv_count+=n_to_recv;

MPI_Probe(p->neighb[5], 4, p->comm, &status);
  MPI_Get_count(&status, p->type, &n_to_recv);
 
 MPI_Recv(bank->p+recv_count, n_to_recv, p->type, p->neighb[5], 4, p->comm, MPI_STATUS_IGNORE);
 recv_count+=n_to_recv;

MPI_Wait(&reqs[4], MPI_STATUS_IGNORE);
MPI_Wait(&reqs[5], MPI_STATUS_IGNORE);

///reinstate particles to alive
for(int i=0; i<recv_count; i++)
{bank->p[i].alive=TRUE;}

///resize bank
////sourcebank->resize=resize_particles???
bank->n = recv_count;
//bank->sz = recv_count; //?
	
return;
}

void synchronize_bank(Bank *source_bank, Bank *fission_bank)
{
  unsigned long i, j;
  unsigned long n_s = source_bank->n;
  unsigned long n_f = fission_bank->n;

  // If the fission bank is larger than the source bank, randomly select
  // n_particles sites from the fission bank to create the new source bank
  if(n_f >= n_s){

    // Copy first n_particles sites from fission bank to source bank
    memcpy(source_bank->p, fission_bank->p, n_s*sizeof(Particle));

    // Replace elements with decreasing probability, such that after final
    // iteration each particle in fission bank will have equal probability of
    // being selected for source bank
    for(i=n_s; i<n_f; i++){
      j = rni(0, i+1);
      if(j<n_s){
        memcpy(&(source_bank->p[j]), &(fission_bank->p[i]), sizeof(Particle));
      }
    }
  }

  // If the fission bank is smaller than the source bank, use all fission bank
  // sites for the source bank and randomly sample remaining particles from
  // fission bank
  else{

    // First randomly sample particles from fission bank
    for(i=0; i<(n_s-n_f); i++){
      j = rni(0, n_f);
      memcpy(&(source_bank->p[i]), &(fission_bank->p[j]), sizeof(Particle));
    }

    // Fill remaining source bank sites with fission bank
    memcpy(&(source_bank->p[n_s-n_f]), fission_bank->p, n_f*sizeof(Particle));
  }

  fission_bank->n = 0;

  return;
}

void calculate_keff(double *keff, double *mean, double *std, int n)
{
  int i;

  *mean = 0;
  *std = 0;

  // Calculate mean
  for(i=0; i<n; i++){
    *mean += keff[i];
  }
  *mean /= n;

  // Calculate standard deviation
  for(i=0; i<n; i++){
    *std += pow(keff[i] - *mean, 2);
  }
  *std = sqrt(*std/(n-1));

  return;
}
