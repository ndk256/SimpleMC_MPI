#include "header.h"

int main(int argc, char *argv[])
{
  Parameters *parameters; // user defined parameters
  Geometry *geometry; // homogenous cube geometry
  Material *material; // problem material
  Bank *source_bank; // array for particle source sites
  Bank *fission_bank; // array for particle fission sites
  Tally *mytally, *global_tally; // scalar flux tallies
  double *keff, *mykeff; // effective multiplication factor
  double t1 =12, t2; // timers

    MPI_Init(&argc, &argv);
  // Establish MPI values
    int prcperdim[3]={0}, periodicity[3]={1,1,1}; ///I'm unsure how this will be good for non-periodic materials but I'll just trust that it is
  MPI_Comm_size(MPI_COMM_WORLD, &(parameters->n_prc));
  
  // Get inputs: set parameters to default values, parse parameter file,
  // override with any command line inputs, and print parameters
  parameters = init_parameters();
  parse_parameters(parameters); ///Are there any issues with reading the same CLI on multiple processes?
  read_CLI(argc, argv, parameters);
  if(parameters->local_rank==0) print_parameters(parameters);
  
  if(parameters->n_prc_auto == TRUE)
  MPI_Dims_create(parameters->n_prc, 3, prcperdim); ///
  else{
  prcperdim[0]=parameters->n_prc_x; 
  prcperdim[1]=parameters->n_prc_y;
  prcperdim[2]=parameters->n_prc_z;
  }
    
  // Set initial RNG seed
  set_initial_seed(parameters->seed);
  set_stream(STREAM_INIT);

  // Create files for writing results to
  if(parameters->local_rank==0) init_output(parameters);

  // Set up geometry
  geometry = init_geometry(parameters);
 /* if(geometry->bc==PERIODIC){
    periodicity[0]=1;
    periodicity[1]=1;
    periodicity[2]=1;}*/
 
  // Create 3D topography and orient self
  MPI_Cart_create(MPI_COMM_WORLD, 3, prcperdim, periodicity, 1, &parameters->comm);
  int mycoords[3];   MPI_Comm_rank(parameters->comm, &(parameters->local_rank));
  MPI_Cart_coords(parameters->comm, parameters->local_rank, 3, mycoords);
  double mybounds[6] = {mycoords[0]*(geometry->x/prcperdim[0]), (mycoords[0]+1)*(geometry->x/prcperdim[0]), mycoords[1]*(geometry->y/prcperdim[1]), (mycoords[1]+1)*(geometry->y/prcperdim[1]),mycoords[2]*(geometry->z/prcperdim[2]), (mycoords[2]+1)*geometry->z/prcperdim[2]};
// ^sorry that's a long line^
  int myneighb[6];
  MPI_Cart_shift(parameters->comm, 0, 1, &(parameters->neighb[0]), &(parameters->neighb[1]));
  MPI_Cart_shift(parameters->comm, 1, 1, &(parameters->neighb[2]), &(parameters->neighb[3]));
  MPI_Cart_shift(parameters->comm, 2, 1, &(parameters->neighb[4]), &(parameters->neighb[5]));
  
  // Set up material
  material = init_material(parameters);

  // Set up tallies
  global_tally = init_tally(parameters);
  
  if(mycoords[0]==0&&mycoords[1]==0&&mycoords[2]==0)
      // Create source bank and initial source distribution
  {source_bank = init_source_bank(parameters, geometry);}
  
  // Set up array for k effective
  keff = calloc(parameters->n_active, sizeof(double));
  
  /// Set up localization
  localize_parameters(&parameters, prcperdim);
  
  // Create (local) fission bank
  fission_bank = init_fission_bank(parameters);
  
  my_tally = init_tally(parameters);
  my_sourcebank = init_bank(parameters->n_particles);
  mykeff = calloc(parameters->n_active, sizeof(double)); ///unsure if correct
  
  if(mycoords[0]==0&&mycoords[1]==0&&mycoords[2]==0)
      // Create source bank and initial source distribution
  {
  distribute_sb(mycoords, parameters, prcperdim, source_bank, &my_sourcebank);
  free_bank(source_bank); ///is this appropriate placement?
  }
  else distribute_sb(mycoords, parameters, prcperdim, my_sourcebank, &my_sourcebank);
  
  //add Barrier here?
  
  if(parameters->local_rank==0){
  center_print("SIMULATION", 79);
  border_print();
  printf("%-15s %-15s %-15s\n", "BATCH", "KEFF", "MEAN KEFF");
  }
  
  // Start time
  if(MPI_WTIME_IS_GLOBAL) t1 = MPI_Wtime(); /// utilize the wall clock if it's implemented
  ///note: test at some point to see if it is implemented/functional, and if so potentially remove the else
  else {MPI_Barrier(parameters->comm); t1 = timer();}
  
  run_eigenvalue(mybounds, parameters, geometry, material, my_sourcebank, fission_bank, mytally, mykeff);    
  
  // Stop time
  if(MPI_WTIME_IS_GLOBAL) {t2 = MPI_Wtime();}
  else {MPI_Barrier(parameters->comm); t2 = timer();}
  
  if(parameters->local_rank==0) printf("Simulation time: %f secs\n", t2-t1);

  // Free memory
  free(keff); free(mykeff);
  free_tally(global_tally); free(mytally);
  free_bank(fission_bank);
  if(!(mycoords[0]==0&&mycoords[1]==0&&mycoords[2]==0)) free_bank(my_sourcebank); //causes segmentation fault otherwise 
  free_material(material);
  free(geometry);
  
  MPI_Barrier(parameters->comm);
    
  MPI_Type_free(&parameters->type);
    MPI_Comm_free(&parameters->comm);
free(parameters);
  
  MPI_Finalize();

  return 0;
}
