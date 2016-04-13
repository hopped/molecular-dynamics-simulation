/*
 * Copyright (c) 2012-2014 Christoph Niethammer <christoph.niethammer@googlemail.com>
 * Copyright (c) 2014      Mohamad Amer Wafai <amerwafai@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

/**
 * @file
 *
 * The main program
 */

#include "config.h"
#include "typedefs.h"

#ifdef GPU
#include "moleculecontainer/basicN2/basicN2.h"
#endif

#include "simulation.h"
#include "phasespace.h"
#include "moleculecontainer/moleculecontainer.h"
#include "ensemble.h"
#include "io/generators/grid_generator.h"
#include "io/generators/dropOverBasin_generator.h"

#include "utils/timer.h"
#include "utils/logger.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
//#include <papi.h>
#ifdef MPI
#include <mpi.h>
#include "mpi_utils.h"
#include "moleculecontainer/basicN2/basicN2.h"
#endif
char generator_name[256];
void print_version_info() {
	LOG(Info, "Christoph's Test Kernel MD, version %s\n", VERSION_STRING);
}

void usage() {
    fprintf(stderr, "main [[OPTIONS]...] [[FILE]...] \n"
"   -N, --num-molecules=ULONG\n"
"       number of molecules\n"
"   -V, --volume=DOUBLE\n"
"   --generator-drop-radius=DOUBLE\n"
"   --generator=NAME\n"
"          allowed names: grid, drop, dropOverBasin\n"
"   --domain-type=TYPE\n"
"         allowed types: cube, rectangle\n"
"   -T, --temperature=DOUBLE\n"
"   -n, --rho=DOUBLE\n"
"   --simulation-time=DOUBLE\n"
"   -v, --verbose=INT\n"
"   --cutoff-radius=DOUBLE\n"
"   --simulation-end-time=DOUBLE\n"
"   --simulation-start-time=DOUBLE\n"
"   --simulation-equilibration-time=INTEGER\n"
"   --n-per-subdomain=INTEGER\n"
"   --cpu-part=INTEGER\n"
"   --thermostat=TYPE\n"
"          allowed types: velocity-scaling\n"
"   --povray-output\n"
"   --ascii-output\n"
"   --molecule-container=TYPE\n"
"          allowed types: MOLECULEBLOCKS, BASICN2\n"
"   --gridgenerator-lattice-system=SYSTEM\n"
"   --gridgenerator-lattice-centering=CENTERING\n"
"   --gridgenerator-random-seed=INT\n"
"   --gridgenerator-dist-boltzmann=(0|1)\n"
"   --gridgenerator-random-shift=DOUBLE\n"
"   --external-force-field=DOUBLE,DOUBLE,DOUBLE\n"
"   -h, --help\n"
            );
}

void print_help() {
    usage();
    exit(0);
}

int main(int argc, char *argv[]) {
    INIT_GLOBAL_LOG("main.log");
    print_version_info();
    int myrank = 0;
    int size = 1;
    phasespace_t *psp = NULL;
    ensemble_t target_ensemble;
    ensemble_init(&target_ensemble);
    domain_t domain;

//    int Events[2] = {PAPI_FP_OPS, PAPI_TOT_INS};
//    long long values[2];

    /* set up as in the NVE MD example for comparison */
    simulation.dt         = 0.005;
    simulation.start_time = 0;
    simulation.end_time   = 0.05; //150;
    simulation.equilibration_time = 50;
    config.density = 0;
    config.cutoff_radius = 1;
    config.cutoff_radius_sq = 1;
    config.externalForceField[0] = 0.;
    config.externalForceField[1] = 0.;
    config.externalForceField[2] = 0.;
    config.lattice_system = cubic;
    config.lattice_centering = primitive;
    config.gridgenerator_random_seed = 1; /* defalt seed, see man srand */
    config.gridgenerator_dist_boltzmann = 0;
    config.gridgenerator_random_shift = -1; /* disabled */
    strcpy(generator_name, "grid");
    config.generator_dropRadius = 0;

    int c;
    while (1) {
        int option_index = 0;
        static struct option long_options[] = {
            {"ascii-output",0,0,0},
            {"cutoff-radius",1,0,0},
            {"external-force-field",1,0,0},
            {"domain-type",1,0,0},
            {"simulation-equilibration-time",1,0,0},
            {"gridgenerator-dist-boltzmann",1,0,0},
            {"gridgenerator-lattice-centering",1,0,0},
            {"gridgenerator-lattice-system",1,0,0},
            {"gridgenerator-random-seed",1,0,0},
            {"gridgenerator-random-shift",1,0,0},
            {"generator-drop-radius",1,0,0},
            {"generator",1,0,0},
            {"help",0,0,'h'},
            {"mass",1,0,'m'},
            {"molecule-container",1,0,0},
            {"num-molecules",1,0,'N'},
            {"povray-output",0,0,0},
            {"rho",1,0,'n'},
            {"simulation-end-time",1,0,0},
            {"simulation-start-time",1,0,0},
            {"temperature",1,0,'T'},
            {"thermostat",1,0,0},
            {"timestep-length",1,0,0},
            {"verbose",2,0,'v'},
            {"volume",1,0,'V'},
            {"n-per-subdomain",1,0,0},
            {"cpu-part",1,0,0},
            {0,0,0,0}
        };

        c = getopt_long (argc, argv, "hN:m:n:V:T:v::",
                long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
            case 0:
                if( strcmp( "cutoff-radius", long_options[option_index].name ) == 0 ) {
                    real rc = atof(optarg);
                    LOG(Info, "Cutoff radius: %"PRIreal"\n", rc);
                    config.cutoff_radius = rc;
                    config.cutoff_radius_sq = rc * rc;
                }
                else if( strcmp( "domain-type", long_options[option_index].name ) == 0 ) {
                    domain_set_type(&domain, optarg);
                    LOG(Info, "Domain type: %s (%d)\n", optarg, domain.type);
                }
                else if( strcmp( "external-force-field", long_options[option_index].name ) == 0 ) {
                    char delimiter[] = ",;:";
                    char *ptr;
                    int d;
                    for( d = 0, ptr = strtok(optarg, delimiter);  d < 3 && ptr != NULL;  d++, ptr = strtok(NULL, delimiter) ) {
                        config.externalForceField[d] = atof(ptr);
                    }
                    LOG(Info, "External force field: [%lf, %lf, %lf]\n", config.externalForceField[0], config.externalForceField[1], config.externalForceField[2]);
                }
                else if( strcmp( "generator", long_options[option_index].name ) == 0 ) {
                    strncpy(generator_name, optarg, 255);
                    generator_name[255] = '\n';
                    LOG(Info, "Generator: %s\n", generator_name);
                }
                else if( strcmp( "generator-drop-radius", long_options[option_index].name ) == 0 ) {
                    config.generator_dropRadius = atof(optarg);
                    LOG(Info, "Generator drop radius: %lf\n", config.generator_dropRadius);
                }
                else if( strcmp( "simulation-start-time", long_options[option_index].name ) == 0 ) {
                    simulation.start_time = (real) atof(optarg);
                    LOG(Info, "Simulation start time: %"PRIreal"\n", simulation.start_time);
                }
                else if( strcmp( "simulation-end-time", long_options[option_index].name ) == 0 ) {
                    simulation.end_time = (real) atof(optarg);
                    LOG(Info, "Simulation end time: %"PRIreal"\n", simulation.end_time);
                }
                else if( strcmp( "simulation-equilibration-time", long_options[option_index].name ) == 0 ) {
                    simulation.equilibration_time = (real) atof(optarg);
                    LOG(Info, "Simulation equilibration time: %"PRIreal"\n", simulation.equilibration_time);
                }
                else if( strcmp( "timestep-length", long_options[option_index].name ) == 0 ) {
                    simulation.dt = (real) atof(optarg);
                    LOG(Info, "Timestep length: %"PRIreal"\n", simulation.dt);
                }
                else if( strcmp( "thermostat", long_options[option_index].name ) == 0 ) {
                    simulation_set_thermostat(&simulation, optarg);
                    LOG(Info, "Thermostat: %s\n", optarg);
                }
                else if( strcmp( "ascii-output", long_options[option_index].name ) == 0 ) {
                    config.ascii_output = 1;
                    LOG(Info, "ASCII output: yes\n");
                }
                else if( strcmp( "povray-output", long_options[option_index].name ) == 0 ) {
                    config.povray_output = 1;
                    LOG(Info, "Povray output: yes\n");
                }
                else if( strcmp( "molecule-container", long_options[option_index].name ) == 0 ) {
                    config.molecule_container_type = get_mc_type(optarg);
                    LOG(Info, "Molecule container type: %s (%d)\n", optarg, config.molecule_container_type);
                }
                else if( strcmp( "gridgenerator-lattice-system", long_options[option_index].name ) == 0 ) {
                    config.lattice_system = lattice_system(optarg);
                    if(config.lattice_system == -1) {
                        LOG(Error, "Unknown lattice system specified: %s\n", optarg);
                        exit(1);
                    }
                    LOG(Info, "Grid generator, lattice system: %s\n", optarg);
                }
                else if( strcmp( "gridgenerator-lattice-centering", long_options[option_index].name ) == 0 ) {
                    config.lattice_centering = lattice_centering(optarg);
                    if(config.lattice_centering == -1) {
                        LOG(Error, "Unknown lattice centering specified: %s\n", optarg);
                        exit(1);
                    }
                    LOG(Info, "Grid generator, lattice centering: %s\n", optarg);
                }
                else if( strcmp( "gridgenerator-random-seed", long_options[option_index].name ) == 0 ) {
                    config.gridgenerator_random_seed = atoi(optarg);
                    LOG(Info, "Grid genenerator random seed: %d\n", config.gridgenerator_random_seed);
                }
                else if( strcmp( "gridgenerator-random-shift", long_options[option_index].name ) == 0 ) {
                    config.gridgenerator_random_shift = atof(optarg);
                    LOG(Info, "Grid genenerator random shift: %d\n", config.gridgenerator_random_shift);
                }
                else if( strcmp( "gridgenerator-dist-boltzmann", long_options[option_index].name ) == 0 ) {
                    config.gridgenerator_dist_boltzmann = atoi(optarg);
                    LOG(Info, "Grid genenerator use Maxwell Boltzmann velocity distribution: %d\n", config.gridgenerator_dist_boltzmann);
                }
                else if( strcmp( "cpu-part", long_options[option_index].name ) == 0 ) {
                    target_ensemble.CPU_PART = atoi(optarg);
                    LOG(Info, "Number of Subdomains running on CPU is: %d\n", target_ensemble.CPU_PART);
                }
                else if( strcmp( "n-per-subdomain", long_options[option_index].name ) == 0 ) {
                    target_ensemble.NPS = atoi(optarg);
                    LOG(Info, "Number of Molecules Per Subdomain: %d\n", target_ensemble.NPS);
                }
                else {
                    LOG(Info, "Option %s", long_options[option_index].name);
                    if (optarg)
                        LOG(Info, " mit Argument %s", optarg);
                    LOG(Info, "\n");
                }
                break;
            case 'h':
                print_help();
                break;
            case 'v': 
                logger_set_output_level(GLOBAL_LOG, Info);
                if (optarg) {
                    LogLevel ll = atoi(optarg);
                    logger_set_output_level(GLOBAL_LOG, ll);
                }
                LOG(Info, "Verbose mode: %d\n", logger_get_output_level(GLOBAL_LOG));
                break;
                     
            case 'm':
                config.m = atof(optarg);
                LOG(Info, "Mass of molecule: %"PRIreal"\n", config.m);
                break;
            case 'N':
                target_ensemble.N = atol(optarg);
                LOG(Info, "Ensemble.N: %ld\n", target_ensemble.N);
                break;
            case 'n':
                config.density = atof(optarg);
                LOG(Info, "rho: %lf\n", config.density);
                break;
            case 'V':
                target_ensemble.V = atof(optarg);
                LOG(Info, "Ensemble.V: %lf\n", target_ensemble.V);
                break;
            case 'T':
                target_ensemble.T = atof(optarg);
                target_ensemble.type = target_ensemble.type | 0x001;
                LOG(Info, "Ensemble.T: %lf\n", target_ensemble.T);
                break;
            case '?':
                break;
            default:
                LOG(Warning, "?? getopt returned character 0%o ??\n", c);
                break;
        }
    }
    if (optind < argc)
    {
        LOG(Info, "Input files:\n");
        while (optind < argc)
            LOG(Info, "  %s\n", argv[optind++]);
    } 
#ifdef MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    target_ensemble.N /= size;
    if ( target_ensemble.N > 0 && config.density > 0) {
        target_ensemble.V = target_ensemble.N / config.density;
        LOG(Info, "Ensemble.V: %lf\n", target_ensemble.V);
    }
    target_ensemble.type = NVT;
    ensemble_print_config(&target_ensemble);


    simulation_init(&simulation);
    real box_length[3];
    box_length[0] = pow(target_ensemble.V, 1./3.);
    box_length[1] = box_length[0];
    box_length[2] = box_length[0];
    domain_set_size(&domain, box_length);
    int i;
    real min=domain.L[0];
    for (i=1;i<3;i++){
	  if(min>domain.L[i]) min = domain.L[i];	
    }
    if(config.cutoff_radius>min/2){
       printf("Cutoff Radius must be <= %"PRIreal"\n",min/2);
       exit(0);
    }

    psp = psp_alloc();
#ifdef MPI
    MPI_Create_dims(myrank);
    MPI_Map(psp, myrank);
//    printf("rank %d has positive neighbors in X: %d, Y: %d, Z: %d\n", ((basicN2 *)( psp->mc->container ))->myrank,
//           ((basicN2 *)( psp->mc->container ))->Nx[0], ((basicN2 *)( psp->mc->container ))->Ny[0], ((basicN2 *)( psp->mc->container ))->Nz[0]);
//    printf("rank %d has negative neighbors in X: %d, Y: %d, Z: %d\n", ((basicN2 *)( psp->mc->container ))->myrank,
//           ((basicN2 *)( psp->mc->container ))->Nx[1], ((basicN2 *)( psp->mc->container ))->Ny[1], ((basicN2 *)( psp->mc->container ))->Nz[1]);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    LOG(Info, "Process %d Generating phase space %s...\n", myrank, generator_name);
//    grid_generator(&target_ensemble, cubic, primitive, &domain, psp);
    if( strcmp(generator_name, "grid") == 0 ) {
        grid_generator(&target_ensemble, config.lattice_system, config.lattice_centering, &domain, psp);
    }
    else if( strcmp(generator_name, "dropOverBasin") == 0 ) {
        double l = domain.L[2];
        double basinHight = 2.*l/5.;
        double dropPos[3] = {l/2.,l/2.,7.*l/10.};
        double dropRadius = l/10.;
        if(config.generator_dropRadius > 0) {
            dropRadius = config.generator_dropRadius;
        }
        LOG(Info, "Basin hight: %lf\n", basinHight);
        LOG(Info, "Drop position: [%lf, %lf, %lf]\n", dropPos[0], dropPos[1], dropPos[2]);
        LOG(Info, "Drop radius: %lf\n", dropRadius);
        dropOverBasin_generator(&target_ensemble, basinHight, dropPos, dropRadius, &domain, psp);
    }
    else if( strcmp(generator_name, "drop") == 0 ) {
        double l = domain.L[2];
        double basinHight = 0.;
        double dropPos[3] = {l/2.,l/2.,l/2.};
        double dropRadius = l/4.;
        if(config.generator_dropRadius > 0) {
            dropRadius = config.generator_dropRadius;
        }
        LOG(Info, "Basin hight: %lf\n", basinHight);
        LOG(Info, "Drop position: [%lf, %lf, %lf]\n", dropPos[0], dropPos[1], dropPos[2]);
        LOG(Info, "Drop radius: %lf\n", dropRadius);
        dropOverBasin_generator(&target_ensemble, basinHight, dropPos, dropRadius, &domain, psp);
    }
    else {
        LOG(Error, "Unknown generator '%s'\n", generator_name);
        exit(1);
    }
    mc_print_stats(psp->mc, &target_ensemble);

    LOG(Info, "Process %d Initial ensemble info ...\n", myrank);
    ensemble_t initial_ensemble;
    ensemble_init(&initial_ensemble);
    initial_ensemble.NPS=target_ensemble.NPS;
    initial_ensemble.CPU_PART=target_ensemble.CPU_PART;
    psp_calc_ensemble_values(psp, &domain, &initial_ensemble);
//    ensemble_print_info(&initial_ensemble);
    if (ENABLED == config.ascii_output) {  
        char filename[256];
        snprintf(filename, 255, "psp-%d_initial.dat", myrank);
        mc_print_ascii(psp->mc, filename);
    }
    if (ENABLED == config.povray_output) {
        char filename[256];
        snprintf(filename, 255, "psp-%d_initial.pov", myrank);
        mc_print_pov(psp->mc, filename);
    }

    real stime;

    size_t timesteps = 0;
    double start_timer, end_timer;
    LOG(Info, "Process %d Calculating inital forces.\n", myrank);
    ensemble_t ensemble;
    ensemble_init(&ensemble);
    ensemble.NPS=target_ensemble.NPS;
    ensemble.CPU_PART=target_ensemble.CPU_PART;
    LOG(Info, "Process %d Update Halos to calculate initial forces.\n", myrank);
    psp_update(psp, &domain);
    psp_update_halo(psp, &domain);
    psp_estimate_halos(psp,&domain);
    psp_decompose(psp, &domain);
    LOG(Info, "Process %d resets forces.\n", myrank);
    psp_reset_forces_and_momenta(psp);
    psp_calc_forces(psp,&ensemble.U_pot);
    psp_integrate_postf(psp, &simulation, &domain);
    LOG(Info, "Process %d Starting main simulation loop.\n", myrank);
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    start_timer = timer();
//    if (PAPI_start_counters (Events, 2) != PAPI_OK) {
//        printf("Error starting counters\n");
//        exit(1);
//    }

    for (stime = simulation.start_time; stime < simulation.end_time; stime += simulation.dt) {
        LOG(Info, "Process %d Current simulation time: %" PRIreal "\n", myrank, stime );

        psp_integrate_pref(psp, &simulation, &domain);
        psp_update(psp, &domain);
        psp_update_halo(psp, &domain);
        if (ENABLED == config.ascii_output) {
            char filename[256];
            snprintf(filename, 255, "psp-%" PRIreal "_%d.dat", stime, myrank);
            mc_print_ascii(psp->mc, filename);
        }

        if (ENABLED == config.povray_output) {
            char filename[256];
            snprintf(filename, 255, "psp-%" PRIreal "_%d.pov", stime, myrank);
            mc_print_pov(psp->mc, filename);
        }
        psp_decompose(psp, &domain);
        ensemble_init(&ensemble);
        ensemble.NPS=target_ensemble.NPS;
        ensemble.CPU_PART=target_ensemble.CPU_PART;
        psp_reset_forces_and_momenta(psp);
        psp_calc_forces(psp, &ensemble.U_pot);
        psp_integrate_postf(psp, &simulation, &domain);
        psp_calc_ensemble_values(psp, &domain, &ensemble);
        if(stime < simulation.equilibration_time && simulation.thermostat_type == VELOCITY_SCALING ) {
            LOG(Info, "Process %d Applying velocity scaling thermostat T=%"PRIreal", T_target =%"PRIreal"\n", myrank ,ensemble.T, target_ensemble.T);
            psp_apply_thermostat(psp, ensemble.T, target_ensemble.T);
            /* recalculate ensemble values as they changed due to the termostat */
            psp_calc_ensemble_values(psp, &domain, &ensemble);
        }

        printf("Process %d %" PRIreal "\t", myrank, stime );
        ensemble_print_info(&ensemble);

        timesteps++;
    }
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
//    if (PAPI_stop_counters (values, 2) != PAPI_OK) {
//        printf("Error stopping counters\n");
//        exit(1);
//    }

    end_timer = timer();

//    printf("Process %d has PAPI_FP_OPS = %lld\n", myrank, values[0]);
//    printf("Process %d has PAPI_TOT_INS = %lld\n", myrank, values[1]);


    if (ENABLED == config.ascii_output) {
        char filename[256];
        snprintf(filename, 255, "psp-final_%d.dat", myrank);
        mc_print_ascii(psp->mc, filename);
    }
    if (ENABLED == config.povray_output) {
        char filename[256];
        snprintf(filename, 255, "psp-final_%d.pov", myrank);
        mc_print_pov(psp->mc, filename);
    }

    ensemble_t final_ensemble;
    ensemble_init(&final_ensemble);
    final_ensemble.NPS=target_ensemble.NPS;
    psp_calc_ensemble_values(psp, &domain, &final_ensemble);
    mc_print_stats(psp->mc, &final_ensemble);
    LOG(Info, "Process %d Computation in the main look took: %lf sec\n", myrank, end_timer - start_timer);
//    printf("Process %d spent %lf sec\n", myrank, end_timer - start_timer);
//    printf("Process %d performs %lg MFLOPS\n", myrank, values[0]*1e-6/(end_timer - start_timer));
//    printf("Process %d performs %lld FOP\n", myrank, values[0]);
//    printf(" %ld\t%lg\n", final_ensemble.N, ( end_timer - start_timer ) / (double)( final_ensemble.N * size ));
    if(myrank==0){
      printf(" %ld\t%lg\n", final_ensemble.N, ( end_timer - start_timer ));
    }
    LOG(Info, "Process %d Computational time per step and molecule: %le sec\n", myrank, ( end_timer - start_timer ) / ( timesteps * final_ensemble.N ));

    psp_destroy(psp);
    psp_free(psp);

    FINALIZE_GLOBAL_LOG();
#ifdef MPI
    MPI_Finalize();
#endif
    return 0;
}
