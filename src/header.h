#ifndef HEADER
#define HEADER

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<sys/time.h>
#include<math.h>
#include<float.h>
#include<unistd.h>
#include<string.h>
#include<mpi.h>

#define TRUE 1
#define FALSE 0

// Constants
#define PI 3.1415926535898
#define D_INF DBL_MAX

// Geometry boundary conditions
#define VACUUM 0
#define REFLECT 1
#define PERIODIC 2

// Reaction types
#define TOTAL 0
#define ABSORPTION 1
#define SCATTER 2
#define FISSION 3

// Surfaces
#define X0 0
#define X1 1
#define Y0 2
#define Y1 3
#define Z0 4
#define Z1 5

// RNG streams
#define N_STREAMS 3
#define STREAM_INIT 0
#define STREAM_TRACK 1
#define STREAM_OTHER 2

typedef struct Parameters_{
  unsigned long long seed; // RNG seed
  unsigned long n_particles; // number of particles
  int n_batches; // number of batches
  int n_generations; // number of generations per batch
  int n_active; // number of active batches
  int bc; // boundary conditions
  int n_nuclides; // number of nuclides in material
  int tally; // whether to tally
  int n_bins; // number of bins in each dimension of mesh
  double nu; // average number of fission neutrons produced
  double xs_a; // absorption macro xs
  double xs_s; // scattering macro xs
  double xs_f; // fission macro xs
  double gx; // geometry size in x
  double gy; // geometry size in y
  double gz; // geometry size in z
  int neighb[6]; // ranks of Cartesian neighbours
  int write_tally; // whether to output tallies
  int write_keff; // whether to output keff
  char *tally_file; // path to write tallies to
  char *keff_file; // path to write keff to
MPI_Datatype type; /// type used to transmit particle data in MPI messages
MPI_Comm comm; ///
int local_rank; ///
int n_prc;
int n_prc_x, n_prc_y, n_prc_z; /// number of processes per x, y, z dimension
int n_prc_auto; /// whether to construct process-dimension setup automatically (with MPI_Dims_create) or use user/file inputs
} Parameters;

typedef struct Particle_{
  int alive;
  double mu; // cosine of polar angle
  double phi; // azimuthal angle
  double u; // direction
  double v;
  double w;
  double x; // position
  double y;
  double z;
  int surface_crossed;
} Particle;

typedef struct Geometry_{
  int bc;
  double x;
  double y;
  double z;
} Geometry;

typedef struct Nuclide_{
  double xs_f; // fission micro xs
  double xs_a; // absorption micro xs
  double xs_s; // scattering micro xs
  double xs_t; // total micro xs
  double atom_density; // atomic density of nuclide in material
} Nuclide;

typedef struct Material_{
  double xs_f; // fission macro xs
  double xs_a; // absorption macro xs
  double xs_s; // scattering macro xs
  double xs_t; // total macro xs
  int n_nuclides;
  Nuclide *nuclides;
} Material;

typedef struct Tally_{
  int tallies_on; // whether tallying is currently turned on
  int n; // mumber of grid boxes in each dimension 
  double dx; // grid spacing
  double dy;
  double dz;
  double *flux;
} Tally;

typedef struct Bank_{
  unsigned long n; // number of particles
  unsigned long sz; // size of bank
  Particle *p; // particle array
  void (*resize)(struct Bank_ *b);
} Bank;

typedef struct Buffer_{
  unsigned int banksz[6], n_banked[6];
  Particle *tox0, *tox1, *toy0, *toy1, *toz0, *toz1;
} Buffer;

// io.c function prototypes
void parse_parameters(Parameters *parameters);
void read_CLI(int argc, char *argv[], Parameters *parameters);
void print_error(char *message);
void print_parameters(Parameters *parameters);
void border_print(void);
void fancy_int(long a);
void center_print(const char *s, int width);
void print_status(int i_a, int i_b, double keff_batch, double keff_mean, double keff_std);
void init_output(Parameters *parameters);
void write_tally(Tally *t, char *filename);
void write_keff(double *keff, int n, char *filename);

// utils.c funtion prototypes
double timer(void);

// prng.c function prototypes
double rn(void);
int rni(int min, int max);
void set_stream(int rn_stream);
void set_initial_seed(unsigned long long rn_seed0);
void rn_skip(long long n);

// initialize.c function prototypes
Parameters *init_parameters(void);
Buffer *init_buff(Parameters *parameters);
Geometry *init_geometry(Parameters *parameters);
Tally *init_tally(Parameters *parameters);
Material *init_material(Parameters *parameters);
Bank *init_fission_bank(Parameters *parameters);
Bank *init_source_bank(Parameters *parameters, Geometry *geometry);
Bank *init_bank(unsigned long n_particles);
void sample_source_particle(Geometry *geometry, Particle *p);
void sample_fission_particle(Particle *p, Particle *p_old);
void resize_particles(Bank *b);
void free_bank(Bank *b);
void free_buf(Buffer *bf);
void free_material(Material *material);
void free_tally(Tally *tally);

///localize.c function prototypes
void localize_parameters(Parameters *par, int dims[6]);
void distrib_particle(int d, int mycoords[3], Parameters *par, int dim, Bank **sb);
void distribute_sb(int mycoords[3], Parameters *p, int nprc[3], Bank *sb, Bank **mysb);

// transport.c function prototypes
double distance_to_boundary(Geometry *geometry, Particle *p);
double distance_to_collision(Material *material);
double dist_to_edge(Particle *p, double s_coords[6]);
void cross_surface(Geometry *geometry, Particle *p);
void collision(Material *material, Bank *fission_bank, double nu, Particle *p);
void transport(Parameters *parameters, Geometry *geometry, double local_bounds[], Material *material, Buffer* sendbuf, Bank *fission_bank, Tally *tally, Particle *p);
void cross_process(double localbounds[], Particle *p, Buffer *sendbuf);

// eigenvalue.c function prototypes
void run_eigenvalue(double localbounds[6], Parameters *parameters, Geometry *geometry, Material *material, Bank *source_bank, Bank *fission_bank, Tally *overall_tally,Tally *tally, double *keff);
void sendrecv_particles(Parameters *p, Bank *bank, Buffer * sendbuf, double mybounds[6]);
void synchronize_bank(Bank *source_bank, Bank *fission_bank);
void calculate_keff(double *keff, double *mean, double *std, int n);

// tally.c function prototypes
void score_tally(Parameters *parameters, Material *material, Tally *tally, Particle *p);
void reset_tally(Tally *tally);

#endif
