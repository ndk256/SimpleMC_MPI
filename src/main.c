#include "header.h"

int main(int argc, char *argv[])
{
  Parameters *parameters; // user defined parameters
  Geometry *geometry; // homogenous cube geometry
  Material *material; // problem material
  Bank *source_bank; // array for particle source sites
  Bank *fission_bank; // array for particle fission sites
  Tally *tally; // scalar flux tally
  double *keff; // effective multiplication factor
  double t1, t2; // timers

  // Get inputs: set parameters to default values, parse parameter file,
  // override with any command line inputs, and print parameters
  parameters = init_parameters();
  parse_parameters(parameters);
  read_CLI(argc, argv, parameters);
  print_parameters(parameters);
  
  MPI_Init(&argc, &argv);
  int prcsize, myrank, prcperdim[3], periodicity[3]={0,0,0};

  // Establish MPI values
  MPI_Comm_size(MPI_COMM_WORLD, &prcsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm cube; //really ought to have a better name
  MPI_Dims_create(prcsize, 3, prcperdim);
  
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
  init_output(parameters);

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
  tally = init_tally(parameters);

  // Create source bank and initial source distribution
  source_bank = init_source_bank(parameters, geometry);

  // Create fission bank
  fission_bank = init_fission_bank(parameters);

  // Set up array for k effective
  keff = calloc(parameters->n_active, sizeof(double));

  center_print("SIMULATION", 79);
  border_print();
  printf("%-15s %-15s %-15s\n", "BATCH", "KEFF", "MEAN KEFF");

  // Start time
  t1 = timer();

  run_eigenvalue(parameters, geometry, material, source_bank, fission_bank, tally, keff);

  // Stop time
  t2 = timer();

  printf("Simulation time: %f secs\n", t2-t1);

  // Free memory
  free(keff);
  free_tally(tally);
  free_bank(fission_bank);
  free_bank(source_bank);
  free_material(material);
  free(geometry);
  free(parameters);
  
  MPI_Type_free(&PARTICLE);
  MPI_Comm_free(&cube);
  MPI_Finalize();

  return 0;
}
