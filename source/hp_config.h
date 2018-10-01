#ifndef HP_CONFIG_H
#define HP_CONFIG_H

#include "config/config.h"

EMP_BUILD_CONFIG( HPConfig,
  GROUP(HARDWARE_GROUP, "Hardware settings"),
  VALUE(TAG_WIDTH, size_t, 16, "How many bits per hardware tag"),
  VALUE(UID,       int,     0, "Position that the UID will be in hw trait vector"),
  VALUE(VOTE,      size_t,  1, "Position that the Vote will be in hw trait vector"),
  VALUE(POSX,      size_t,  2, "Position that the Coordinate X will be in hw trait vector"),
  VALUE(POSY,      size_t,  3, "Position that the Coordinate Y will be in hw trait vector"),
  VALUE(FOE,       size_t,  4, "Friend or foe hardware?"),
  VALUE(MAX_CORES, size_t, 20, "Maximum number of cores a hardware can spawn."),
  GROUP(GRAPH_GROUP, "Graph settings"),
  VALUE(GRA_DIM,  size_t,       3, "Dimension of graph"),
  VALUE(NUM_ITER, size_t,     150, "Number of iterations per trial."),
  VALUE(NUM_FRI,  size_t,       8, "Number of friends in the graph."),
  VALUE(NUM_ENE,  size_t,       1, "Number of enemies in the graph."),
  VALUE(GRA_TYPE, size_t,       0, "Type of graph we are about to use."),
  VALUE(MIN_BND,  size_t,       1, "Lower bound on random numbers."),
  VALUE(MAX_BND,  size_t, 1000000, "Uper bound on the random numbers."),
  GROUP(MUTATION_GROUP, "Mutation settings"),
  VALUE(MIN_FUN_CNT, size_t,   3, "Minimum number of functions each hardware should have."),
  VALUE(MAX_FUN_CNT, size_t,  10, "Maximum number of functions each hardware should have."),
  VALUE(MIN_FUN_LEN, size_t,   10, "Minimum number of instructions each function will have."),
  VALUE(MAX_FUN_LEN, size_t,  64, "Maximum number of instructions each function will have."),
  VALUE(MAX_TOT_LEN,  size_t, 512, "Maximum size of hardware genome."),
  VALUE(MIN_BIN_THSH, size_t, 0.0, "Minimum Threshold for "),
  GROUP(EXPERIMENT_GROUP, "Experiment settings"),
  VALUE(POP_SIZE,  size_t, 200, "Population size."),
  VALUE(NUM_GENS,  size_t, 10000, "Number of generations per experiments."),
  VALUE(RNG_SEED,  size_t,    17, "Random number seed."),
  VALUE(EVAL_SIZE, size_t,     5, "Number of bad guys a good guy will face per run."),
  VALUE(TOURN_SIZE, size_t,    4, "Number or organims competing in tournament selection."),
  VALUE(SNAP_SHOT,  size_t,   50, "Time that we will take a snapshot of population"),
  VALUE(COHORT_SIZE, size_t, 20, "Number of agents in a cohort"),
  VALUE(COHORT_TOTAL, size_t, 10, "Number of cohorts that exist"),
  VALUE(RESULTS_DIR, std::string, "No_Adress", "Address of where the results will go to!")
)

#endif
