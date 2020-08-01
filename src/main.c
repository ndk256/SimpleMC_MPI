#include "header.h"

int main(int argc, char *argv[])
{
    Parameters *parameters; // user defined parameters
    Geometry *geometry; // homogenous cube geometry
    Material *material; // problem material
    Bank *fission_bank; // array for particle fission sites
    Tally *mytally, *globaltally; // scalar flux tally
    double *keff, *mykeff; // effective multiplication factor
    double t1=12, t2; // timers

MPI_Init(&argc, &argv); ///
 
  // Get inputs: set parameters to default values, parse parameter file,
    // override with any command line inputs, and print parameters
    parameters = init_parameters();
    parse_parameters(parameters);
    read_CLI(argc, argv, parameters);

int prcperdim[3]={0}, periodicity[3]={1,1,1}; ///{0,0,0};

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

MPI_Barrier(parameters->comm);
int mycoords[3];
MPI_Cart_coords(parameters->comm, parameters->local_rank, 3, mycoords);
if(parameters->local_rank==0)  print_parameters(parameters);

///Defining boundaries of this process
double mybounds[6] = {mycoords[0]*(geometry->x/prcperdim[0]), (mycoords[0]+1)*(geometry->x/prcperdim[0]), mycoords[1]*(geometry->y/prcperdim[1]), (mycoords[1]+1)*(geometry->y/prcperdim[1]),mycoords[2]*(geometry->z/prcperdim[2]), (mycoords[2]+1)*geometry->z/prcperdim[2]}; ///is this correct?

///Getting the adjacent processes
MPI_Cart_shift(parameters->comm, 0, 1, &(parameters->neighb[0]), &(parameters->neighb[1]));
MPI_Cart_shift(parameters->comm, 1, 1, &(parameters->neighb[2]), &(parameters->neighb[3]));
MPI_Cart_shift(parameters->comm, 2, 1, &(parameters->neighb[4]), &(parameters->neighb[5]));

  // Set initial RNG seed
  set_initial_seed(parameters->seed); ///is it an issue to have this occurring more than once (eg, on each process)?
  set_stream(STREAM_INIT);

  // Create files for writing results to
if(parameters->local_rank==0){
  init_output(parameters);}

Bank *my_sourcebank; /// (local) array for particle source sites

  // Set up material
  material = init_material(parameters);

if(mycoords[0]==0&&mycoords[1]==0&&mycoords[2]==0){
  my_sourcebank = init_source_bank(parameters, geometry);
}
else
    my_sourcebank = init_bank(parameters->n_particles); ///

  // Set up array for k effective
  keff = calloc(parameters->n_active, sizeof(double));

    ///set up tallies///
mytally = init_tally(parameters);
    globaltally = init_tally(parameters);
    
// Localize parameters (dimensions)
localize_parameters(parameters, prcperdim);

  // Create source bank and initial source distribution
distribute_sb(mycoords, parameters, prcperdim, my_sourcebank, &my_sourcebank);

    // Create (local) fission bank
  fission_bank = init_fission_bank(parameters);

  ///set up a local keff array///
mykeff = calloc(parameters->n_active,sizeof(double));

MPI_Barrier(parameters->comm);

if(parameters->local_rank==0){
  center_print("SIMULATION", 79);
  border_print();
  printf("%-15s %-15s %-15s\n", "BATCH", "KEFF", "MEAN KEFF");

 // Start timer
 if(MPI_WTIME_IS_GLOBAL)
    t1 = MPI_Wtime();
    else t1 = timer();
}

  run_eigenvalue(mybounds, parameters, geometry, material, my_sourcebank, fission_bank, globaltally, mytally, mykeff);

MPI_Barrier(parameters->comm);

if(parameters->local_rank==0){
  // Stop timer
  if(MPI_WTIME_IS_GLOBAL) t2 =MPI_Wtime();
    else t2 = timer();

  printf("Simulation time: %f secs\n", t2-t1);
}

  // Free memory
  free(keff); free(mykeff); 
 free_tally(mytally); free_tally(globaltally);
free_bank(fission_bank); 
   free_bank(my_sourcebank); 

 free_material(material);
  free(geometry);

MPI_Barrier(parameters->comm);

MPI_Type_free(&(parameters->type));
MPI_Comm_free(&parameters->comm);

free(parameters); 

MPI_Finalize(); ///

  return 0;
}
