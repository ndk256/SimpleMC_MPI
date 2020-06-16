#include "header.h"

int main(int argc, char *argv[])
{
  Parameters *parameters, *local_par; // user defined parameters
  Geometry *geometry; // homogenous cube geometry
  Material *material; // problem material
  Bank *source_bank; // array for particle source sites
  Bank *fission_bank; // array for particle fission sites
  Tally *mytally, *global_tally; // scalar flux tallies
  double *keff, *mykeff; // effective multiplication factor
  double t1, t2; // timers

  // Establish MPI values
  MPI_Comm_size(MPI_COMM_WORLD, &prcsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm cube; //really ought to have a better name
  MPI_Dims_create(prcsize, 3, prcperdim);
  
  // Get inputs: set parameters to default values, parse parameter file,
  // override with any command line inputs, and print parameters
  parameters = init_parameters();
  parse_parameters(parameters); ///Are there any issues with reading the same CLI on multiple processes?
  read_CLI(argc, argv, parameters);
  if(myrank==0) print_parameters(parameters);
  local_par = localize_parameters(parameters, prcperdim);

  MPI_Init(&argc, &argv);
  int prcsize, myrank, prcperdim[3], periodicity[3]={0,0,0};
  
  // Define Particle in MPI-sendable format
  MPI_Datatype PARTICLE;
  MPI_Datatype base[2] = {MPI_INT, MPI_DOUBLE};
  int blocks[2] = {2,8};
  MPI_Aint offsets[2], extent;
  MPI_Type_extent(MPI_INT, &extent);
  offsets[0] =0; offsets[1]=blocks[0]*extent;
  MPI_Type_struct(2, blocks, offsets, base, &PARTICLE);
  MPI_Type_commit(&PARTICLE);
  
  // Set initial RNG seed
  set_initial_seed(parameters->seed);
  set_stream(STREAM_INIT);

  // Create files for writing results to
  if(myrank==0) init_output(parameters);

  // Set up geometry
  geometry = init_geometry(parameters);
  if(geometry->bc==PERIODIC){
    periodicity[0]=1;
    periodicity[1]=1;
    periodicity[2]=1;}
 
  // Create 3D topography and orient self
  MPI_Cart_create(MPI_COMM_WORLD, 3, prcperdim, periodicity, 1, &cube);
  int mycoords[3];
  MPI_Cart_coords(cube, myrank, 3, mycoords);
  double mybounds[6] = {mycoords[0]*(geometry->x/prcperdim[0]), (mycoords[0]+1)*(geometry->x/prcperdim[0]), mycoords[1]*(geometry->y/prcperdim[1]), (mycoords[1]+1)*(geometry->y/prcperdim[1]),mycoords[2]*(geometry->z/prcperdim[2]), (mycoords[2]+1)*geometry->z/prcperdim[2]};
// ^sorry that's a long line^
  int myneighb[6];
  MPI_Cart_shift(cube, 0, -1, &myrank, &myneighb[0]);
  MPI_Cart_shift(cube, 0, 1, &myrank, &myneighb[1]);
  MPI_Cart_shift(cube, 1, -1, &myrank, &myneighb[2]);
  MPI_Cart_shift(cube, 1, 1, &myrank, &myneighb[3]);
  MPI_Cart_shift(cube, 2, -1, &myrank, &myneighb[4]);
  MPI_Cart_shift(cube, 2, 1, &myrank, &myneighb[5]);
  
  // Set up material
  material = init_material(parameters);

  // Set up tallies
  global_tally = init_tally(parameters);
  my_tally = init_tally(local_par);

  // Create source bank and initial source distribution
  my_sourcebank = init_source_bank(parameters, geometry); ///unsure of arguments
  if(mycoords[0]==0&&mycoords[1]==0&&mycoords[2]==0)
  {source_bank = init_source_bank(parameters, geometry);
  distribute_sb(mycoords, prcperdim, myrank, spatialgrid, source_bank, my_sourcebank);
  }
  else distribute_sb(mycoords, prcperdim, myrank, spatialgrid, my_sourcebank, my_sourcebank);
  
  ///free source_bank right now??
  
  // Create (local) fission bank
  fission_bank = init_fission_bank(parameters);

  // Set up array for k effective
  keff = calloc(parameters->n_active, sizeof(double));
  mykeff = calloc(local_par->n_active,sizeof(double)); ///unsure if correct
  
  if(myrank==0){
  center_print("SIMULATION", 79);
  border_print();
  printf("%-15s %-15s %-15s\n", "BATCH", "KEFF", "MEAN KEFF");
  }
  
  // Start time
  if(MPI_WTIME_IS_GLOBAL) t1 = MPI_Wtime(); /// utilize the wall clock if it's implemented
  ///note: test at some point to see if it is implemented/functional, and if so potentially remove the else
  else {MPI_Barrier(cube); t1 = timer();}
  
  run_eigenvalue(cube, myrank, myneighb, mybounds, local_par, geometry, material, my_sourcebank, fission_bank, mytally, mykeff);

  ///aggregate the tally///
    MPI_Reduce(mytally->flux, global_tally->flux, parameters->n_particles, MPI_DOUBLE, MPI_SUM, 0, cube);
    
  
  // Stop time
  if(MPI_WTIME_IS_GLOBAL) {t2 = MPI_Wtime();}
  else {MPI_Barrier(cube); t2 = timer();}

  if(myrank==0) printf("Simulation time: %f secs\n", t2-t1);

  // Free memory
  free(keff); free(mykeff);
  free_tally(global_tally); free(mytally);
  free_bank(fission_bank);
  if(myrank==0) free_bank(source_bank); ///may change
  free_material(material);
  free(geometry);
  free(parameters); free(local_par);
  
  MPI_Type_free(&PARTICLE);
  MPI_Comm_free(&cube);
  MPI_Finalize();

  return 0;
}
