#include "header.h"

int main(int argc, char *argv[])
{
    Parameters *parameters; // user defined parameters
    Geometry *geometry; // homogenous cube geometry
    Material *material; // problem material
    Bank *source_bank; // array for particle source sites
    Bank *fission_bank; // array for particle fission sites
    Tally *global_tally, *mytally; // scalar flux tally
    double *keff, *mykeff; // effective multiplication factor
    double t1=12, t2; // timers

MPI_Init(&argc, &argv); ///
 
  // Get inputs: set parameters to default values, parse parameter file,
    // override with any command line inputs, and print parameters
    parameters = init_parameters();
    parse_parameters(parameters);
    read_CLI(argc, argv, parameters);

int prcperdim[3]={0}, periodicity[3]={1,1,1}; ///{0,0,0};
 ///possibly reconsider these namings for greater clarity

  // Set up geometry
   geometry = init_geometry(parameters);
//if(geometry->bc==PERIODIC){periodicity[0]=1;periodicity[1]=1;periodicity[2]=1;}

///Establishing basic MPI necessity values and creating a 3D cartography
MPI_Comm_size(MPI_COMM_WORLD, &(parameters->n_prc)); ///

if(parameters->n_prc_auto)
MPI_Dims_create(parameters->n_prc, 3, prcperdim); ///
else{
prcperdim[0]=parameters->n_prc_x; prcperdim[1]=parameters->n_prc_y; 
prcperdim[2]=parameters->n_prc_z;
}

MPI_Cart_create(MPI_COMM_WORLD, 3, prcperdim, periodicity, 1, &(parameters->comm));
MPI_Comm_rank(parameters->comm, &(parameters->local_rank)); //

MPI_Barrier(parameters->comm); //or commworld?
int mycoords[3];
MPI_Cart_coords(parameters->comm, parameters->local_rank, 3, mycoords);
if(parameters->local_rank==0)  print_parameters(parameters);

///Defining boundaries of this process
double mybounds[6] = {mycoords[0]*(geometry->x/prcperdim[0]), (mycoords[0]+1)*(geometry->x/prcperdim[0]), mycoords[1]*(geometry->y/prcperdim[1]), (mycoords[1]+1)*(geometry->y/prcperdim[1]),mycoords[2]*(geometry->z/prcperdim[2]), (mycoords[2]+1)*geometry->z/prcperdim[2]}; ///is this correct?

///Getting the adjacent processes
for(int i=0; i<6; i++) {parameters->neighb[i]=0;} ////TEMP////
MPI_Cart_shift(parameters->comm, 0, 1, &(parameters->neighb[0]), &(parameters->neighb[1]));
MPI_Cart_shift(parameters->comm, 1, 1, &(parameters->neighb[2]), &(parameters->neighb[3]));
MPI_Cart_shift(parameters->comm, 2, 1, &(parameters->neighb[4]), &(parameters->neighb[5]));

  // Set initial RNG seed
  set_initial_seed(parameters->seed); ///is it an issue to have this occurring more than once (eg, on each process)?
  set_stream(STREAM_INIT);

  // Create files for writing results to
if(parameters->local_rank==0){
  init_output(parameters);}

Bank *my_sourcebank;
///idk if the below is good
my_sourcebank = init_bank(parameters->n_particles); ///

  // Set up material
  material = init_material(parameters);

  // Set up overall tally
 global_tally = init_tally(parameters);

if(mycoords[0]==0&&mycoords[1]==0&&mycoords[2]==0){
  source_bank = init_source_bank(parameters, geometry);
}

  // Set up array for k effective
  keff = calloc(parameters->n_active, sizeof(double));

// Localize parameters (dimensions)
localize_parameters(parameters, prcperdim);

///set up local tallies///
mytally = init_tally(parameters);

  // Create source bank and initial source distribution
if(mycoords[0]==0&&mycoords[1]==0&&mycoords[2]==0){
distribute_sb(mycoords, parameters, prcperdim, source_bank, &my_sourcebank);
printf("distrib_sb\n");
//free_bank(source_bank);
}
else distribute_sb(mycoords, parameters, prcperdim, my_sourcebank, &my_sourcebank);
//printf("\nand... n=%d", my_sourcebank->n);
  // Create (local) fission bank
  fission_bank = init_fission_bank(parameters);

  ///set up a local keff array///
mykeff = calloc(parameters->n_active,sizeof(double)); ///can i just do this?

MPI_Barrier(parameters->comm);

if(parameters->local_rank==0){
  center_print("SIMULATION", 79);
  border_print();
  printf("%-15s %-15s %-15s\n", "BATCH", "KEFF", "MEAN KEFF");

 // Start time
 // t1 = timer();
 t1 = MPI_Wtime();
}

///here we go
  run_eigenvalue(mybounds, parameters, geometry, material, my_sourcebank, fission_bank, mytally, mykeff);
//NOTE: segmentation fault somewhere in the above function
 ///but only at i_p == 1176

//printf("%d:eig complete!\n",parameters->local_rank);
MPI_Barrier(parameters->comm);

if(parameters->local_rank==0){
  // Stop time
///^that's pretty dramatic
  t2 =MPI_Wtime();// timer();

  printf("Simulation time: %f secs\n", t2-t1);
}

//  MPI_Type_free(&parameters->type);

  // Free memory
  free(keff); free(mykeff); 
 free_tally(mytally); free_tally(global_tally);  
free_bank(fission_bank); 
if(!(mycoords[0]==0&&mycoords[1]==0&&mycoords[2]==0)) 
   free_bank(my_sourcebank); 
else free_bank(source_bank);

 free_material(material);
  free(geometry);

MPI_Barrier(parameters->comm);

MPI_Type_free(&(parameters->type));
MPI_Comm_free(&parameters->comm);

free(parameters); 

MPI_Finalize(); ///

  return 0;
}
