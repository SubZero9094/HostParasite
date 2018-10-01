#ifndef HP_EXPERIMENT_H
#define HP_EXPERIMENT_H

//TODO: 
// Write to all output files for each genome
// Write all stats in the csv
// Finish python program to do analysis

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "hp_config.h"
#include "ParaGraph.h"
#include "../../Empirical/source/tools/Random.h"
#include "../../Empirical/source/tools/random_utils.h"
#include "../../Empirical/source/hardware/EventDrivenGP.h"
#include "../../Empirical/source/Evolve/World.h"
#include "../../Alex/signalgp-adventures/adventures/utility_belt/source/utilities.h"
#include "../../Empirical/source/Evolve/Lex_Cohort_Selection.h" //Include my selection function.


/* CONSTEXPR FOR HARDWARE */

//Number of bits in the tag per hardware
constexpr size_t TAG_WIDTH = 16;
constexpr double VALUE = .0000001;

const std::string NOP_PATH = "genome1.txt";
const std::string BASIC_PATH = "genome.txt";
const std::string ADVANCE_PATH = "genome2.txt";

/* NEW TYPE DECLARATIONS FOR HARDWARE*/

//Type for a actural hardware
using hardware_t = emp::EventDrivenGP_AW<TAG_WIDTH>;
//Type for the hardware genome
using program_t = hardware_t::Program;
//Type for the hardwares state 
using state_t = hardware_t::State;
//Instruction object for hardware
using inst_t = hardware_t::inst_t;
//Instruction library for hardware
using inst_lib_t = hardware_t::inst_lib_t;
//Type for Events for each hardware
using event_t = hardware_t::event_t;
//Event Libarary for hardware
using event_lib_t = hardware_t::event_lib_t;
//Function type for the hardware
using function_t = hardware_t::Function;
//Memory type for hardware
using memory_t = hardware_t::memory_t;
//Mutator
using mutant_t = toolbelt::SignalGPMutator<hardware_t>;

/* NEW TYPE DECLARATIONS FOR SIMPLICITY*/

using coor_t = std::pair<size_t, size_t>;


class Experiment
{
  struct Agent;

  //World of Agents
  using world_t = emp::World<Agent>;  

  struct Agent
  {
    //Without Parasite
    double mWithout = -100;

    //With Parasite
    double mBest = -100;
    double mWorst = 1000000000;
    double mAverage = -100;
    std::vector<double> mScores;
    double mMedian = -100;

    //DNA
    program_t mGenome;

    Agent(const program_t & p) : mGenome(p) {;}
    program_t & GetGenome() {return mGenome;}

    void PrintInfo()
    {
      std::cout << "Agent Vector: ";
      for(auto x : mScores)
        std::cout << x << " ";
      std::cout << std::endl;

      std::cout << "Best Without: " << mWithout << std::endl;
      std::cout << "Best Score: " << mBest << std::endl;
      std::cout << "Worst Score: " << mWorst << std::endl;
      std::cout << "Average: " << mAverage << std::endl;
      std::cout << "Median: " << mMedian << std::endl;
    }

    //mWithout
    void Without(double x)
    {
      mWithout = x;
    }

    //mAverage
    //mMedian
    //mBest
    //mWorst
    void CalcMetrics(size_t size)
    {
      if(mScores.size() != size)
      {
        std::cout << "Agent scores not correct size!" << std::endl;
        std::cout << "Size: " << size << " mScores.size(): " << mScores.size() << std::endl;
        exit(-1);
      }

      //Calculating average
      double sum = 0.0;
      double denom = (double) size;
      for(auto x : mScores)
        sum += x;
      mAverage = sum / denom;

      //Calculating median
      auto temp = mScores;
      std::sort(temp.begin(), temp.end());
      mWorst = temp[0];
      mBest = temp[temp.size() - 1];
      //even
      if(size%2 == 0)
      {
        mMedian = (temp[size/2] + temp[(size/2) - 1]) / 2.0;
      }
      else
      {
        mMedian = temp[size/2];
      }

    }

    void Insert(double x)
    {
      mScores.push_back(x);
    }

    double GetMean()
    {
      if(mAverage == -100)
      {
        std::cout << "mAverage not updated" << std::endl;
        exit(-1);
      }
      return mAverage;
    }

    double GetMedian()
    {
      if(mMedian == -100)
      {
        std::cout << "mMedian not updated" << std::endl;
        exit(-1);
      }
      return mMedian;
    }

    double GetBest()
    {
      if(mBest == -100)
      {
        std::cout << "mBest not updated" << std::endl;
        exit(-1);
      }
      return mBest;
    }

    double GetWorst()
    {
      if(mWorst == 3.14)
      {
        std::cout << "mWorst not updated" << std::endl;
        exit(-1);
      }
      return mWorst;
    }

    double GetWithout()
    {
      if(mWithout == -100)
      {
        std::cout << "mWithout not updated" << std::endl;
        exit(-1);
      }
      return mWithout;
    }

    double GetScore(size_t i)
    {
      if(i >= mScores.size())
      {
        std::cout << "Agent GetScore: " << i << " > mScores.size(): " <<  mScores.size() << std::endl;
        std:: 
        exit(-1);
      }
      return mScores[i];
    }

  };

  public:
    Experiment(const HPConfig & config) :
    POP_SIZE(config.POP_SIZE()), NUM_GENS(config.NUM_GENS()),
    RNG_SEED(config.RNG_SEED()), EVAL_SIZE(config.EVAL_SIZE()),
    TOURN_SIZE(config.TOURN_SIZE()), GRA_DIM(config.GRA_DIM()),
    GRA_TYPE(config.GRA_TYPE()), SNAP_SHOT(config.SNAP_SHOT()),
    NUM_ITER(config.NUM_ITER()), NUM_ENE(config.NUM_ENE()), 
    NUM_FRI(config.NUM_FRI()), COHORT_SIZE(config.COHORT_SIZE()),  
    COHORT_TOTAL(config.COHORT_TOTAL()), MIN_FUN_CNT(config.MIN_FUN_CNT()),
    MAX_FUN_CNT(config.MAX_FUN_CNT()), MIN_FUN_LEN(config.MIN_FUN_LEN()), 
    MAX_FUN_LEN(config.MAX_FUN_LEN()), MAX_TOT_LEN(config.MAX_TOT_LEN()),
    RESULTS_DIR(config.RESULTS_DIR())
    {
      mRng = emp::NewPtr<emp::Random>(RNG_SEED);
      inst_lib = emp::NewPtr<inst_lib_t>();
      event_lib = emp::NewPtr<event_lib_t>();
      mGraph = emp::NewPtr<ParaGraph>(config, mRng);
      mGoodWorld = emp::NewPtr<world_t>(*mRng, "GoodWorld");
      mBadWorld = emp::NewPtr<world_t>(*mRng, "BadWorld");
      mMutant = emp::NewPtr<mutant_t>(MIN_FUN_CNT, MAX_FUN_CNT, MIN_FUN_LEN, MAX_FUN_LEN, MAX_TOT_LEN);
      THEORY_MAX_NOPARA = NUM_ITER * GRA_DIM * GRA_DIM;
      THEORY_MAX_PARA = NUM_ITER * (NUM_FRI);

      mGoodCohort.resize(COHORT_TOTAL);
      mBadCohort.resize(COHORT_TOTAL);

      for(auto & co : mGoodCohort)
        co.resize(COHORT_SIZE);
      
      for(auto & co : mBadCohort)
        co.resize(COHORT_SIZE);

      for(size_t i = 0; i < POP_SIZE; ++i)
      {
        mGoodPos.push_back(i);
        mBadPos.push_back(i);
      }

      SetVote(config.VOTE()); SetUID(config.UID());
      SetPOSX(config.POSX()); SetPOSY(config.POSY());

      mBestPos.resize(TOTAL_POS);

      // create and open the .csv file
      //mCSV.open(FILE_NAME);

      //mCSV << "Generation" << "," << "Best Host" << "," << "Best Host -> Parasite" << "," << "Best Parasite" << "," << "Best Parasite -> Host" << std::endl;
    }

    ~Experiment()
    {
      mGraph.Delete();
      mGoodWorld.Delete();
      mBadWorld.Delete();
      mRng.Delete();
      mMutant.Delete();
      inst_lib.Delete();
      event_lib.Delete();
    }

    /* FUNCTIONS DEDICATED TO THE EXPERIMENT */

    //Run the experiment
    void Run();

    //Confiugre all the neccesary things
    //Configure instructions, events, create world of pop size with nop genome
    void Config_All();

    //Evalute each agent for 
    std::vector<double> Evaluation_step();

    //Selection
    void Selection_step();

    //Update
    void Update_step();

    //Return an genome full of nops
    program_t Genome_NOP();

    //Return basic program
    program_t Genome_BASIC();

    //Return advance program
    program_t Genome_ADVANCE();

    //Calculate the averages of the population
    std::vector<double> CalcMetrics(size_t size, double without);

    std::vector<double> StoreResults(double & max_mean, double & min_mean, double & max_med, double & min_med, double & best, double & worst, double without)
    {
      std::vector<double> results; results.resize(6);
      results[RMAX_MEAN] = max_mean;
      results[RMIN_MEAN] = min_mean;
      results[RMAX_MEDIAN] = max_med;
      results[RMIN_MEDIAN] = min_med;
      results[RMAX_BEST] = best;
      results[RMAX_WORST] = worst;
      results[RMAX_WITHOUT] = without;
      return results;
    }

    /* FUNCTIONS DEDICATED TO THE CONFIGURATIONS, INSTRUCTIONS, EVENTS*/

    //Will make the instruction library
    void Config_Inst();

    //Will set the output memory of a hardware to all of its neighboors
    static void Inst_BroadcastMail(hardware_t & hw, const inst_t & inst);

    //Will send the vote to all of a hardwares neighboors
    static void Inst_BroadcastVote(hardware_t & hw, const inst_t & inst);

    //Will load the UID of a hardware into its working buffer
    static void Inst_GetUID(hardware_t & hw, const inst_t & inst);
    
    //Will get the vote of a hardware
    static void Inst_GetVote(hardware_t & hw, const inst_t & inst);

    //Will set the vote of the hardware
    static void Inst_SetVote(hardware_t & hw, const inst_t & inst);

    //Will make the event library
    void Config_Events();

    //Will actually do the event
    void Dispatch_Broadcast(hardware_t & hw, const event_t & e);

    //Will do the event to send vote out
    void Dispatch_BroadcastVote(hardware_t & hw, const event_t & e);

    //Will spawn a core for the event
    static void Handle_Broadcast(hardware_t & hw, const event_t & e);

    //Will create the hardware and add it to the world
    void Config_HW(program_t p);

    //Will configure the world
    void Config_World();

    //Will create the cohort and shuffle them
    void Create_Cohort();

    //Will create the fitness function vector
    void Create_FitFun();


    /* FUNCTIONS DEDICATED TO SETTERS */
    static void SetVote(size_t x) {VOTE = x;}
    static void SetUID(size_t x) {UID = x;}
    static void SetPOSX(size_t x) {POSX = x;}
    static void SetPOSY(size_t x) {POSY = x;}


    /* FUNCTIONS DEDICATED TO TEST ParaGraph */
    
    //Will test if the graph is actually GRA_DIM X GRA_DIM
    //PASS
    void GraphTest1();

    //Will test the adjacency list and create enemeies
    //PASS
    void GraphTest2();

    //Will test setting up initial traits
    //PASS
    void GraphTest3();

    //Will test setting an enemy node
    //PASS
    void GraphTest4();

    //Will test to see if the FinalVotes vector is correct
    //PASS
    void GraphTest5();    

    //Will test to se if scoring functions are correct
    //PASS
    void GraphTest6();

    //Will test loading GENOME
    //PASS
    void GraphTest7();

    //Will test step by step process of the hardware
    //PASS
    void GraphTest8();

    //Will test to see if system without parasite will work
    //PASS
    void GraphTest9();

    //Will test the cohort initialization and testing
    //PASS
    void GraphTest10();

    //Will test the CalcMetrics function
    void GraphTest11();

    //Will test the Single Selection function created
    void GraphTest12();

    //Will test the overall selection function
    void GraphTest13();

  private:

    /* EXPERIMENT SPECIFIC PARAMATERS */

    //Size of the population we are evaluating.
    size_t POP_SIZE;
    //Number of genearations evaluating during experiments
    size_t NUM_GENS;
    //Random Number Seed Number
    size_t RNG_SEED;
    //Number of bad guys a good guy will face per Run
    size_t EVAL_SIZE;
    //Number of organisms competing during tournament do_selection_step
    size_t TOURN_SIZE;
    //Dimension of the graph
    size_t GRA_DIM;
    //Pointer for random number generator
    emp::Ptr<emp::Random> mRng;
    //Pointer for Graph
    emp::Ptr<ParaGraph> mGraph;
    //Graph tyep
    size_t GRA_TYPE;
    //vector to hold positions of hardware
    std::vector<size_t> mHwPos;
    //World to hold good organisms
    emp::Ptr<world_t> mGoodWorld;
    //World to hold bad organisms
    emp::Ptr<world_t> mBadWorld;
    //Mutator
    emp::Ptr<mutant_t> mMutant;
    //Snapshot index
    size_t SNAP_SHOT;
    //All the iterations
    size_t NUM_ITER;
    //Count of how many time broadcast vote is called
    double mCount;
    //Theoretical max with parasite
    size_t THEORY_MAX_PARA;
    //Theoretical max without parasite
    size_t THEORY_MAX_NOPARA;
    //Number of enemeies in the system
    size_t NUM_ENE;
    //Number of friends in the system
    size_t NUM_FRI;
    //Container of indexes for good world
    emp::vector<size_t> mGoodPos;
    //Container of indexes for bad world
    emp::vector<size_t> mBadPos;
    //Container of Good Cohort
    emp::vector<emp::vector<size_t>> mGoodCohort;
    //Container of Bad Cohort
    emp::vector<emp::vector<size_t>> mBadCohort;
    //Container for all the fitness functions for the good cohort
    emp::vector<std::function<double(Agent&)>> mGFitFuns;
    //Container for all the fitness functions for the bad cohort
    emp::vector<std::function<double(Agent&)>> mBFitFuns;
    //Size of each cohort
    size_t COHORT_SIZE;
    //Total number of cohorts
    size_t COHORT_TOTAL;

    /* HARDWARE SPECIFIC PARAMATERS */

    //Position of UID within hw trait vector
    static size_t UID;
    //Position of VOTE within hw trait vector
    static size_t VOTE;
    //Position of X coordinate within hw trait vector
    static size_t POSX;
    //Position of Y coordinate within hw trait vector
    static size_t POSY;
    //Instruction Library for hardware
    emp::Ptr<inst_lib_t> inst_lib;
    //Event Library for hardware
    emp::Ptr<event_lib_t> event_lib;
    //Directory of BASIC genome
    std::string GENOME_BASIC = "BASIC.txt";
    //Directory of NOP genome
    std::string GENOME_NOP = "NOP.txt";
    //Directory of ADVANCE genome
    std::string GENOME_ADVANCE = "ADVANCE.txt";

    /* MUTATION SPECIFIC PARAMATERS */
    size_t MIN_FUN_CNT;
    size_t MAX_FUN_CNT;
    size_t MIN_FUN_LEN;
    size_t MAX_FUN_LEN;
    size_t MAX_TOT_LEN;

    /* POSITIONS FOR RESULTS */
    size_t RMAX_MEAN = 0;
    size_t RMIN_MEAN = 1;
    size_t RMAX_MEDIAN = 2;
    size_t RMIN_MEDIAN = 3;
    size_t RMAX_BEST = 4;
    size_t RMAX_WORST = 5;
    size_t RMAX_WITHOUT = 6;

    /* POSITIONS FOR BEST PERFORMERS */
    size_t POS_WITHOUT = 0;
    size_t POS_WITH_HOST = 1;
    size_t POS_WITH_PARA = 2;

    /* UPDATE ACCORDING TO POSITIONS FOR RESULTS AND PERFORMERS */
    size_t TOTAL_RESULTS = 7;
    size_t TOTAL_POS = 3;

    /* VARIABLES TO HELP RECORD THE DATA*/

    //Will store all the positions of the best performers
    emp::vector<size_t> mBestPos;

    //File path for the folder we are storing data in
    std::string RESULTS_DIR;

    //Actual file we are adding too
    std::ofstream mGoodDNAs_p;
    std::ofstream mBadDNAs_p;
    std::ofstream mGoodWitoutDNA;
};

size_t Experiment::UID = 0;
size_t Experiment::VOTE = 0;
size_t Experiment::POSX = 0;
size_t Experiment::POSY = 0;

/* FUNCTIONS DEDICATED TO THE EXPERIMENT */

void Experiment::Run()
{
  Experiment::Config_All();

  for(size_t i = 0; i < NUM_GENS; ++i)
  {
    std::cout << "GEN: " << i << std::endl;
    auto v = Experiment::Evaluation_step();

    std::cout << "POPMAX_Mean: " << v[RMAX_MEAN] << " POPMIN_Mean: " << v[RMIN_MEAN] << " POPMAX_Median: "\
    << v[RMAX_MEDIAN] << " POPMIN_Median: " << v[RMIN_MEDIAN] << " POPMAX_Best: " << v[RMAX_BEST] << " POPMAX_Worst: " << v[RMAX_WORST] << " POPMAX_Host_WO: " << v[RMAX_WITHOUT] << std::endl;
    std::cout << std::endl;

    if((i%50) == 0)
    {
      mGoodWorld->GetOrg(0).GetGenome().PrintProgramFull();
    }


    Experiment::Selection_step();
    Experiment::Update_step();
  }
}

//Evalute each agent for 
std::vector<double> Experiment::Evaluation_step()
{
  //first will be the best score with no parasite
  //second will the the best score with parasite
  double best_without = -100;
  double best_with = -100;
  double best_para = -100;
  Create_Cohort();
  //Select on of the cohorts
  for(size_t i = 0; i < COHORT_TOTAL; ++i)
  {
    //Select a specific host within the cohort
    for(size_t j = 0; j < COHORT_SIZE; ++j)
    {
      //Grab good organism and put in the grid and then run it.
      mGraph->Reset();
      Agent & good_agent = mGoodWorld->GetOrg(mGoodCohort[i][j]);
      program_t & good_pro = good_agent.GetGenome();
      mGraph->SetPara(false);
      mGraph->GenomeNoPara(good_pro);
      double without = mGraph->RunNoPara(NUM_ITER);
      good_agent.Without(without);

      if(without > best_without)
      {
        best_without = without;
        mBestPos[POS_WITHOUT] = mGoodCohort[i][j];
      }

      //Face all parasites within that cohort
      for(size_t x = 0; x < COHORT_SIZE; ++x)
      {
        //Grab bad agent and get its program
        Agent & bad_agent = mBadWorld->GetOrg(mBadCohort[i][x]);
        program_t & bad_pro = bad_agent.GetGenome();

        //Reset the graph and load both good and bad organism
        mGraph->Reset();
        mGraph->SetPara(true);
        mGraph->CreatePara(NUM_ENE);
        mGraph->GenomePara(good_pro, bad_pro);
        double score = mGraph->RunPara(NUM_ITER);

        //update positions
        if(score > best_with)
        {
          best_with = score;
          mBestPos[POS_WITH_HOST] = mGoodCohort[i][j];
        }
        if((THEORY_MAX_PARA - score) > best_para)
        {
          best_para = THEORY_MAX_PARA - score;
          mBestPos[POS_WITH_PARA] = mBadCohort[i][x];
        }

        //Store data within the organisms
        good_agent.Insert(score);
        bad_agent.Insert(THEORY_MAX_PARA - score);
      }
    }
  }

  return CalcMetrics(COHORT_SIZE, best_without);
}

std::vector<double> Experiment::CalcMetrics(size_t size, double without)
{
  auto pop = mGoodWorld->GetFullPop();
  for(auto agent : pop )
  {
    agent->CalcMetrics(size);
  }

  double best_med = pop[0]->GetMedian();
  double worst_med = pop[0]->GetMedian();
  double best_avg = pop[0]->GetMean();
  double worst_avg = pop[0]->GetMean();
  double best_score = pop[0]->GetBest();
  double bestw_score = pop[0]->GetWorst();

  if(pop.size() == 1)
  {
    return StoreResults(best_avg, worst_avg, best_med, worst_med, best_score, bestw_score, without);
  }

  else
  {
    for(size_t i = 1; i < pop.size(); ++i)
    {
      if(best_med < pop[i]->GetMedian())
        best_med = pop[i]->GetMedian();

      if(worst_med > pop[i]->GetMedian())
        worst_med = pop[i]->GetMedian();

      if(best_avg < pop[i]->GetMean())
        best_avg = pop[i]->GetMean();

      if(worst_avg > pop[i]->GetMean())
        worst_avg = pop[i]->GetMean();

      if(best_score < pop[i]->GetBest())
        best_score = pop[i]->GetBest();

      if(bestw_score < pop[i]->GetWorst())
        bestw_score = pop[i]->GetWorst();
    }

    return StoreResults(best_avg, worst_avg, best_med, worst_med, best_score, bestw_score, without);
  }
}


//Selection
//TODO
void Experiment::Selection_step()
{
  emp::CohortHostSelect(*mGoodWorld, mGoodCohort, mGFitFuns, COHORT_SIZE, COHORT_TOTAL, mRng);
  emp::CohortHostSelect(*mBadWorld, mBadCohort, mBFitFuns, COHORT_SIZE, COHORT_TOTAL, mRng);
}

//Update
void Experiment::Update_step()
{
  mGoodWorld->Update();
  mGoodWorld->DoMutations();

  mBadWorld->Update();
  mBadWorld->DoMutations();
}

//Return an genome full of nops
program_t Experiment::Genome_NOP()
{
  program_t ancestor_prog(inst_lib);
  std::ifstream ancestor_fstream(GENOME_NOP);

  if(!ancestor_fstream.is_open())
  {
    std::cout << "Failed to open program file: " << GENOME_NOP << " ... Exiting..." << std::endl;
    exit(-1);
  }

  ancestor_prog.Load(ancestor_fstream);
  return ancestor_prog;
}

//Return basic program
program_t Experiment::Genome_BASIC()
{
  program_t ancestor_prog(inst_lib);
  std::ifstream ancestor_fstream(GENOME_BASIC);

  if(!ancestor_fstream.is_open())
  {
    std::cout << "Failed to open program file: " << GENOME_BASIC << " ... Exiting..." << std::endl;
    exit(-1);
  }

  ancestor_prog.Load(ancestor_fstream);
  return ancestor_prog;
}

//Return advance program
program_t Experiment::Genome_ADVANCE()
{
  program_t ancestor_prog(inst_lib);
  std::ifstream ancestor_fstream(GENOME_ADVANCE);

  if(!ancestor_fstream.is_open())
  {
    std::cout << "Failed to open program file: " << GENOME_ADVANCE << " ... Exiting..." << std::endl;
    exit(-1);
  }

  ancestor_prog.Load(ancestor_fstream);
  return ancestor_prog;
}

/* FUNCTIONS DEDICATED TO THE CONFIGURATIONS, INSTRUCTIONS, EVENTS */

void Experiment::Config_All()
{
  Experiment::Config_Events();
  Experiment::Config_Inst();
  Experiment::Config_World();
  Experiment::Create_FitFun();

  program_t nop = Experiment::Genome_NOP();

  mGoodWorld->Inject(nop, POP_SIZE);
  mBadWorld->Inject(nop, POP_SIZE);

  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);
  mGraph->ConfigureTraits();
  mGraph->CreateAdj(GRA_DIM, GRA_TYPE);
}

//Will make the instruction library
void Experiment::Config_Inst()
{
  //Experiment specific instructions
  inst_lib->AddInst("GetUID", Inst_GetUID, 1, "UID => Local Memory Arg1");
  inst_lib->AddInst("Broadcast", Inst_BroadcastMail, 1, "Output Memory => hw.neighboors");

  //inst_lib->AddInst("BroadcastVote", Inst_BroadcastVote, 1, "Output Memory[Arg1] = Vote => hw.neighboors");
  inst_lib->AddInst("GetVote", Inst_GetVote, 1, "Vote => Local Memory Arg1");
  inst_lib->AddInst("SetVote", Inst_SetVote, 1, "Local Memory Arg1 => Hw.Trait[VOTE]");

  // - Setup the instruction set. -
  // Standard instructions:
  inst_lib->AddInst("Inc", hardware_t::Inst_Inc, 1, "Increment value in local memory Arg1");
  inst_lib->AddInst("Dec", hardware_t::Inst_Dec, 1, "Decrement value in local memory Arg1");
  inst_lib->AddInst("Not", hardware_t::Inst_Not, 1, "Logically toggle value in local memory Arg1");
  inst_lib->AddInst("Add", hardware_t::Inst_Add, 3, "Local memory: Arg3 = Arg1 + Arg2");
  inst_lib->AddInst("Sub", hardware_t::Inst_Sub, 3, "Local memory: Arg3 = Arg1 - Arg2");
  inst_lib->AddInst("Mult", hardware_t::Inst_Mult, 3, "Local memory: Arg3 = Arg1 * Arg2");
  inst_lib->AddInst("Div", hardware_t::Inst_Div, 3, "Local memory: Arg3 = Arg1 / Arg2");
  inst_lib->AddInst("Mod", hardware_t::Inst_Mod, 3, "Local memory: Arg3 = Arg1 % Arg2");
  inst_lib->AddInst("TestEqu", hardware_t::Inst_TestEqu, 3, "Local memory: Arg3 = (Arg1 == Arg2)");
  inst_lib->AddInst("TestNEqu", hardware_t::Inst_TestNEqu, 3, "Local memory: Arg3 = (Arg1 != Arg2)");
  inst_lib->AddInst("TestLess", hardware_t::Inst_TestLess, 3, "Local memory: Arg3 = (Arg1 < Arg2)");
  inst_lib->AddInst("If", hardware_t::Inst_If, 1, "Local memory: If Arg1 != 0, proceed; else, skip block.", emp::ScopeType::BASIC, 0, {"block_def"});
  inst_lib->AddInst("While", hardware_t::Inst_While, 1, "Local memory: If Arg1 != 0, loop; else, skip block.", emp::ScopeType::BASIC, 0, {"block_def"});
  inst_lib->AddInst("Countdown", hardware_t::Inst_Countdown, 1, "Local memory: Countdown Arg1 to zero.", emp::ScopeType::BASIC, 0, {"block_def"});
  inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "Close current block if there is a block to close.", emp::ScopeType::BASIC, 0, {"block_close"});
  inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "Break out of current block.");
  inst_lib->AddInst("Call", hardware_t::Inst_Call, 0, "Call function that best matches call affinity.", emp::ScopeType::BASIC, 0, {"affinity"});
  inst_lib->AddInst("Return", hardware_t::Inst_Return, 0, "Return from current function if possible.");
  inst_lib->AddInst("SetMem", hardware_t::Inst_SetMem, 2, "Local memory: Arg1 = numerical value of Arg2");
  inst_lib->AddInst("CopyMem", hardware_t::Inst_CopyMem, 2, "Local memory: Arg1 = Arg2");
  inst_lib->AddInst("SwapMem", hardware_t::Inst_SwapMem, 2, "Local memory: Swap values of Arg1 and Arg2.");
  inst_lib->AddInst("Input", hardware_t::Inst_Input, 2, "Input memory Arg1 => Local memory Arg2.");
  inst_lib->AddInst("Output", hardware_t::Inst_Output, 2, "Local memory Arg1 => Output memory Arg2.");
  inst_lib->AddInst("Commit", hardware_t::Inst_Commit, 2, "Local memory Arg1 => Shared memory Arg2.");
  inst_lib->AddInst("Pull", hardware_t::Inst_Pull, 2, "Shared memory Arg1 => Local memory Arg2.");
  inst_lib->AddInst("Nop", hardware_t::Inst_Nop, 0, "No operation.");
}

//Will set the output memory of a hardware to all of its neighboors
void Experiment::Inst_BroadcastMail(hardware_t & hw, const inst_t & inst)
{
  state_t & state = hw.GetCurState();
  hw.TriggerEvent("BroadcastMail", inst.affinity, state.output_mem);
}

//Will set the output memory of a hardware to all of its neighboors
void Experiment::Inst_BroadcastVote(hardware_t & hw, const inst_t & inst)
{
  double vote = hw.GetTrait(VOTE);
  state_t & state = hw.GetCurState();
  state.SetOutput(inst.args[0], vote);
  hw.TriggerEvent("BroadcastVote", inst.affinity, state.output_mem);
}

//Will load the UID of a hardware into its working buffer
void Experiment::Inst_GetUID(hardware_t & hw, const inst_t & inst)
{
  state_t & state = hw.GetCurState();
  state.SetLocal(inst.args[0], hw.GetTrait(UID));
}

//Will get the vote of a hardware
void Experiment::Inst_GetVote(hardware_t & hw, const inst_t & inst)
{
  state_t & state = hw.GetCurState();
  state.SetLocal(inst.args[0], hw.GetTrait(VOTE));
}

//Will set the vote of the hardware
void Experiment::Inst_SetVote(hardware_t & hw, const inst_t & inst)
{
  state_t & state = hw.GetCurState();
  double vote = state.GetLocal(inst.args[0]);
  hw.SetTrait(VOTE, vote);
}

//Will make the event library
void Experiment::Config_Events()
{
  event_lib->AddEvent("BroadcastMail", Experiment::Handle_Broadcast, "Send output memory to all neighbors.");
  event_lib->RegisterDispatchFun("BroadcastMail", [this](hardware_t & hw, const event_t & e)
  {
    this->Dispatch_Broadcast(hw, e);
  });


  event_lib->AddEvent("BroadcastVote", Experiment::Handle_Broadcast, "Send output memory to all neighbors.");
  event_lib->RegisterDispatchFun("BroadcastVote", [this](hardware_t & hw, const event_t & e)
  {
    auto graph = this->mGraph;
    double vote = hw.GetTrait(VOTE);

    if(graph->Find(vote))
    {
      this->mCount+= 1;
    }


    this->Dispatch_Broadcast(hw, e);
  });
}

//Will spawn a core for the event
void Experiment::Handle_Broadcast(hardware_t & hw, const event_t & e)
{
  hw.SpawnCore(e.affinity, hw.GetMinBindThresh(), e.msg);
}

//Will actually do the event
void Experiment::Dispatch_Broadcast(hardware_t & hw, const event_t & e)
{
  coor_t p = std::make_pair(hw.GetTrait(POSX), hw.GetTrait(POSY));
  auto team = mGraph->GetNodeNeig(p.first, p.second);

  for(auto pair : team)
  {
    auto node = mGraph->GetNode(pair.first, pair.second);
    node->mHW->QueueEvent(e);
  }
}

//Will configure the world
//TODO
void Experiment::Config_World()
{
  mGoodWorld->Reset();
  mGoodWorld->SetPopStruct_Mixed(true);
  mGoodWorld->SetMutFun([this](Agent & agent, emp::Random & rnd)
  {
    program_t & pro = agent.mGenome;
    return this->mMutant->ApplyMutations(pro, rnd);
  });

  mBadWorld->Reset();
  mBadWorld->SetPopStruct_Mixed(true);
  mBadWorld->SetMutFun([this](Agent & agent, emp::Random & rnd)
  {
    program_t & pro = agent.mGenome;
    return this->mMutant->ApplyMutations(pro, rnd);
  });
}

void Experiment::Create_Cohort()
{
  emp::Shuffle(*mRng, mBadPos);
  emp::Shuffle(*mRng, mGoodPos);

  for(size_t i = 0; i < COHORT_TOTAL; ++i)
  {
    for(size_t j = 0; j < COHORT_SIZE; ++j)
    {
      mGoodCohort[i][j] = mGoodPos[(i*20) + j];
      mBadCohort[i][j] = mBadPos[(i*20) + j];
    }
  }
}

//Will create the fitness function vector
void Experiment::Create_FitFun()
{
  //Will fill the vector with COHORT_SIZE many functions and the function without any parasites
  mGFitFuns.resize(COHORT_SIZE + 1);
  mBFitFuns.resize(COHORT_SIZE);

  for(size_t i = 0; i < COHORT_SIZE; ++i)
  {
    mGFitFuns[i] = [i](Agent & agent)
    {
      return agent.GetScore(i);
    };

    mBFitFuns[i] = [i](Agent & agent)
    {
      return agent.GetScore(i);
    };
  }

  mGFitFuns[COHORT_SIZE] = [](Agent & agent)
  {
    return agent.GetWithout();
  };

}


/* FUNCTIONS DEDICATED TO TEST ParaGraph */
  
//Will test if the graph is actually GRA_DIM X GRA_DIM
void Experiment::GraphTest1()
{
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);
  auto g = mGraph->GetGraph();

  if(g.size() != GRA_DIM)
    std::cout << "WRONG" << std::endl;

  for(auto col : g)
  {
    if(col.size() != GRA_DIM)
      std::cout << "WRONG" << std::endl;
  }
  std::cout << "Done with GraphTest1!" << std::endl;

}

//Will test the adjacency list
void Experiment::GraphTest2()
{
  //Create adj list
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);
  mGraph->CreateAdj(GRA_DIM, GRA_TYPE);
  mGraph->PrintFriends();
}

//Print out the mRandomNums and hardware traits
void Experiment::GraphTest3()
{
  //Create adj list
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);
  mGraph->CreateAdj(GRA_DIM, GRA_TYPE);
  mGraph->ConfigureTraits();
  mGraph->PrintRanNums();
  mGraph->PrintUV();
}

//Will test setting an parasite node
void Experiment::GraphTest4()
{
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);
  mGraph->CreateAdj(GRA_DIM, GRA_TYPE);
  mGraph->ConfigureTraits();

  //Create the enemies
  mGraph->CreatePara(NUM_ENE);
  mGraph->PrintFOE();
}

//Will test to see if the FinalVotes vector is correct
void Experiment::GraphTest5()
{
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);
  mGraph->CreateAdj(GRA_DIM, GRA_TYPE);
  mGraph->ConfigureTraits();

  std::cout << "RUN1->All nodes are in concensus" << std::endl;
  auto uid = mGraph->GetRandNums();
  mGraph->SetVotes(uid[0]);
  mGraph->MakeFinalVotes();
  mGraph->PrintRanNums();
  mGraph->PrintFinVote();
  mGraph->PrintUV();

  std::cout << "RUN2->All nodes except one illegal vote are in concensus" << std::endl;
  mGraph->ConfigureTraits();
  uid = mGraph->GetRandNums();
  mGraph->SetVotes(uid[0]);
  mGraph->SetVote(0,0,0);
  mGraph->MakeFinalVotes();
  mGraph->PrintRanNums();
  mGraph->PrintFinVote();
  mGraph->PrintUV();

  std::cout << "RUN3->All nodes except one leagal vote are in concensus" << std::endl;
  mGraph->ConfigureTraits();
  uid = mGraph->GetRandNums();
  mGraph->SetVotes(uid[0]);
  mGraph->SetVote(1,1,uid[1]);
  mGraph->MakeFinalVotes();
  mGraph->PrintRanNums();
  mGraph->PrintFinVote();
  mGraph->PrintUV();
}

//Will test to see if scoring functions are correct
void Experiment::GraphTest6()
{
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);
  mGraph->CreateAdj(GRA_DIM, GRA_TYPE);
  mGraph->ConfigureTraits();

  std::cout << "RUN1_NOPARA-> All nodes are in concensus" << std::endl;
  mGraph->SetPara(false);
  auto uid = mGraph->GetRandNums();
  mGraph->SetVotes(uid[0]);
  mGraph->MakeFinalVotes();
  mGraph->PrintUV();
  std::cout << "CONSENSUS-> SHOULD BE: " << NUM_ENE + NUM_FRI << " GUESS: " << mGraph->Consensus() << std::endl;
  std::cout << "LARGEST-V-> SHOULD BE: " << NUM_ENE + NUM_FRI << " GUESS: " << mGraph->LargestLegal() << std::endl;
  std::cout << "LEGALVOTE-> SHOULD BE: " << NUM_ENE + NUM_FRI << " GUESS: " << mGraph->LegalVotes() << std::endl;
  std::cout << "-----------------------------------------------\n" << std::endl;

  mGraph->ConfigureTraits();
  mGraph->CreatePara(NUM_ENE);
  mGraph->SetPara(true);
  std::cout << "RUN1-PARA> All nodes are in concensus, NO PARA" << std::endl;
  uid = mGraph->GetRandNums();
  mGraph->SetVotes(uid[0]);
  mGraph->MakeFinalVotes();
  mGraph->PrintUV();
  std::cout << "CONSENSUS-> SHOULD BE: " << NUM_FRI << " GUESS: " << mGraph->Consensus() << std::endl;
  std::cout << "LARGEST-V-> SHOULD BE: " << NUM_FRI << " GUESS: " << mGraph->LargestLegal() << std::endl;
  std::cout << "LEGALVOTE-> SHOULD BE: " << NUM_FRI << " GUESS: " << mGraph->LegalVotes() << std::endl;
  std::cout << "-----------------------------------------------\n" << std::endl;

  mGraph->ConfigureTraits();
  std::cout << "RUN2_NOPARA-> All nodes are in concensus except one legal vote" << std::endl;
  mGraph->SetPara(false);
  uid = mGraph->GetRandNums();
  mGraph->SetVotes(uid[0]);
  mGraph->SetVote(0,0,uid[1]);
  mGraph->MakeFinalVotes();
  mGraph->PrintUV();
  std::cout << "CONSENSUS-> SHOULD BE: " << 0 << " GUESS: " << mGraph->Consensus() << std::endl;
  std::cout << "LARGEST-V-> SHOULD BE: " << NUM_ENE + NUM_FRI - 1 << " GUESS: " << mGraph->LargestLegal() << std::endl;
  std::cout << "LEGALVOTE-> SHOULD BE: " << NUM_ENE + NUM_FRI << " GUESS: " << mGraph->LegalVotes() << std::endl;
  std::cout << "-----------------------------------------------\n" << std::endl;

  mGraph->ConfigureTraits();
  mGraph->CreatePara(NUM_ENE);
  mGraph->SetPara(true);
  std::cout << "RUN2-PARA> All nodes are in concensus except one legal vote" << std::endl;
  uid = mGraph->GetRandNums();
  mGraph->SetVotes(uid[0]);
  mGraph->SetVote(1,1,uid[1]);
  mGraph->MakeFinalVotes();
  mGraph->PrintUV();
  std::cout << "CONSENSUS-> SHOULD BE: " << 0 << " GUESS: " << mGraph->Consensus() << std::endl;
  std::cout << "LARGEST-V-> SHOULD BE: " << NUM_FRI - 1 << " GUESS: " << mGraph->LargestLegal() << std::endl;
  std::cout << "LEGALVOTE-> SHOULD BE: " << NUM_FRI << " GUESS: " << mGraph->LegalVotes() << std::endl;
  std::cout << "-----------------------------------------------\n" << std::endl;

  mGraph->ConfigureTraits();
  std::cout << "RUN3_NOPARA-> All nodes are in concensus except one illegal vote" << std::endl;
  mGraph->SetPara(false);
  uid = mGraph->GetRandNums();
  mGraph->SetVotes(uid[0]);
  mGraph->SetVote(0,0,0);
  mGraph->MakeFinalVotes();
  mGraph->PrintUV();
  std::cout << "CONSENSUS-> SHOULD BE: " << 0 << " GUESS: " << mGraph->Consensus() << std::endl;
  std::cout << "LARGEST-V-> SHOULD BE: " << NUM_ENE + NUM_FRI - 1 << " GUESS: " << mGraph->LargestLegal() << std::endl;
  std::cout << "LEGALVOTE-> SHOULD BE: " << NUM_ENE + NUM_FRI - 1 << " GUESS: " << mGraph->LegalVotes() << std::endl;
  std::cout << "-----------------------------------------------\n" << std::endl;

  mGraph->ConfigureTraits();
  mGraph->CreatePara(NUM_ENE);
  mGraph->SetPara(true);
  std::cout << "RUN3-PARA> All nodes are in concensus except one legal vote" << std::endl;
  uid = mGraph->GetRandNums();
  mGraph->SetVotes(uid[0]);
  mGraph->SetVote(1,1,0);
  mGraph->MakeFinalVotes();
  mGraph->PrintUV();
  std::cout << "CONSENSUS-> SHOULD BE: " << 0 << " GUESS: " << mGraph->Consensus() << std::endl;
  std::cout << "LARGEST-V-> SHOULD BE: " << NUM_FRI - 1 << " GUESS: " << mGraph->LargestLegal() << std::endl;
  std::cout << "LEGALVOTE-> SHOULD BE: " << NUM_FRI-1 << " GUESS: " << mGraph->LegalVotes() << std::endl;
  std::cout << "-----------------------------------------------\n" << std::endl;
}

//Will test loading GENOME
void Experiment::GraphTest7()
{
  Experiment::Config_Inst();
  Experiment::Config_Events();

  std::cout << "NO PARASITE IN SYSTEM" << std::endl;
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);
  mGraph->ConfigureTraits();
  mGraph->CreateAdj(GRA_DIM, GRA_TYPE);

  program_t good = Experiment::Genome_NOP();
  std::cout << std::endl;

  mGraph->SetPara(false);
  mGraph->GenomeNoPara(good);
  mGraph->PrintGenomes();

  std::cout << "\nPARASITE IN THE SYSTEM" << std::endl;
  mGraph->Reset();
  mGraph->ConfigureTraits();
  program_t bad = Experiment::Genome_BASIC();
  mGraph->SetPara(true);
  mGraph->CreatePara(NUM_ENE);
  mGraph->GenomePara(good, bad);
  mGraph->PrintGenomes();

  std::cout << "\n OTHER PARASITE IN THE SYSTEM" << std::endl;
  mGraph->Reset();
  mGraph->ConfigureTraits();
  program_t meh = Experiment::Genome_ADVANCE();
  mGraph->SetPara(true);
  mGraph->CreatePara(NUM_ENE);
  mGraph->GenomePara(meh, good);
  mGraph->PrintGenomes();
}

//Will test step by step process of the hardware
void Experiment::GraphTest8()
{
  std::cout << "GRAPHTEST8: \n" << std::endl;

  Experiment::Config_Inst();
  Experiment::Config_Events();
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);
  
  std::cout << "NO PARASITE IN SYSTEM" << std::endl;  
  mGraph->ConfigureTraits();
  mGraph->CreateAdj(GRA_DIM, GRA_TYPE);
  mGraph->SetPara(false);
  program_t advance = Experiment::Genome_ADVANCE();
  mGraph->Reset();
  mGraph->GenomeNoPara(advance);

  mGraph->PrintFriends();

  double score = mGraph->RunNoPara(NUM_ITER);

  std::cout << "NO PARASITE IN THE SYSTEM SCORE : " << score << " THEO MAX: " << THEORY_MAX_NOPARA << std::endl;

  std::cout << "PARASITE IN THE SYSTEM!!!!" << std::endl;

  mGraph->Reset();
  mGraph->SetPara(true);
  mGraph->CreatePara(NUM_ENE);
  program_t bad = Experiment::Genome_NOP();
  mGraph->GenomePara(advance, bad);

  double s = mGraph->RunPara(NUM_ITER);
  std::cout << "PARASITE IN THE SYSTEM SCORE: " << s << " THEO MAX: " << THEORY_MAX_PARA << std::endl;

}

//Will test the system without parasites
void Experiment::GraphTest9()
{
  std::cout << "GRAPHTEST8: \n" << std::endl;

  Experiment::Config_Inst();
  Experiment::Config_Events();
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);
  
  std::cout << "NO PARASITE IN SYSTEM" << std::endl;  
  mGraph->ConfigureTraits();
  mGraph->CreateAdj(GRA_DIM, GRA_TYPE);
  mGraph->SetPara(false);
  program_t advance = Experiment::Genome_ADVANCE();
  mGraph->Reset();
  mGraph->GenomeNoPara(advance);

  std::cout << "Score: " << mGraph->RunNoPara(NUM_ITER);
}

//Test to see if cohorts have correct things in them
void Experiment::GraphTest10()
{
  Experiment::Create_Cohort();

  if(mGoodCohort.size() != COHORT_TOTAL)
  {
    std::cout << "Good cohort does not have correct number of cohorts!" << std::endl;
    exit(0);
  }

  if(mBadCohort.size() != COHORT_TOTAL)
  {
    std::cout << "Bad cohort does not have correct number of cohorts!" << std::endl;
    exit(0);
  }

  for(auto v : mGoodCohort)
  {
    if(v.size() != COHORT_SIZE)
    {
      std::cout << "Good cohort not correct size!" << std::endl;
      exit(0);
    }
  }

  for(auto v : mBadCohort)
  {
    if(v.size() != COHORT_SIZE)
    {
      std::cout << "Bad cohort not correct size!" << std::endl;
      exit(0);
    }
  }

  std::cout << "GOOD COHORT" << std::endl;
  for(size_t i = 0; i < COHORT_TOTAL; ++i)
  {
    std::cout << mGoodCohort[i].size() << std::endl;
    std::cout << i << ":";
    for(size_t j = 0; j < COHORT_SIZE; ++j)
    {
      std::cout << " " << mGoodCohort[i][j];
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;

  std::cout << "BAD COHORT" << std::endl;
  for(size_t i = 0; i < COHORT_TOTAL; ++i)
  {
    std::cout << mBadCohort[i].size() << std::endl;
    std::cout << i << ":";
    for(size_t j = 0; j < COHORT_SIZE; ++j)
    {
      std::cout << " " << mBadCohort[i][j];
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;

  Experiment::Create_Cohort();

  std::cout << "GOOD COHORT" << std::endl;
  for(size_t i = 0; i < COHORT_TOTAL; ++i)
  {
    std::cout << mGoodCohort[i].size() << std::endl;
    std::cout << i << ":";
    for(size_t j = 0; j < COHORT_SIZE; ++j)
    {
      std::cout << " " << mGoodCohort[i][j];
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;

  std::cout << "BAD COHORT" << std::endl;
  for(size_t i = 0; i < COHORT_TOTAL; ++i)
  {
    std::cout << mBadCohort[i].size() << std::endl;
    std::cout << i << ":";
    for(size_t j = 0; j < COHORT_SIZE; ++j)
    {
      std::cout << " " << mBadCohort[i][j];
    }
    std::cout << std::endl;
  }

}

//Will test the CalcMetrics function
void Experiment::GraphTest11()
{
  Experiment::Config_Inst();
  Experiment::Config_Events();
  auto p = Genome_NOP();
  size_t size = 10;
  mGoodWorld->Inject(p, 3);

  for(size_t i = 0; i < 3; ++i)
  {
    Agent & temp = mGoodWorld->GetOrg(i);
    for(size_t j = 0; j < size; ++j)
    {
      temp.Insert(mRng->GetUInt(100));
    }
  }

  auto v = CalcMetrics(size, 10.0);
 
  for(size_t i = 0; i < mGoodWorld->size(); ++i)
  {
    std::cout << "ORG " << i << ": " << std::endl;
    mGoodWorld->GetOrg(i).PrintInfo();
    std::cout << std::endl;
  }

  std::cout << "Best Mean: " << v[RMAX_MEAN] << " Worst Mean " << v[RMIN_MEAN] << " Best Median: " << v[RMAX_MEDIAN] << " Worst Median: " << v[RMIN_MEDIAN] << " Best best: " << v[RMAX_BEST] << " Best Worst: " << v[RMAX_WORST] << std::endl;  

}

//Will test the single host selection algorithm
void Experiment::GraphTest12()
{
  Experiment::Config_Inst();
  Experiment::Config_Events();
  auto p = Genome_NOP();
  std::vector<size_t> vec;
  mGoodWorld->Inject(p, 10);

  for(size_t i = 0; i < mGoodWorld->GetSize(); ++i)
  {
    vec.push_back(i);
    auto & org = mGoodWorld->GetOrg(i);
    org.mBest = i;
  }
  //mGoodWorld->GetOrg(8).mBest = 9;
  std::function<double(Agent&)> f = [](Agent & agent)
    {
      return agent.mBest;
    };
  auto v = emp::SingleHostSelection(*mGoodWorld, vec, f);

  for(size_t i = 0; i < mGoodWorld->GetSize(); ++i)
  {
    std::cout << i << ": " << mGoodWorld->GetOrg(i).mBest << std::endl;
  }

  for(auto i : v)
  {
    std::cout << i << std::endl;
  }

  
}

//Will test overall selection
void Experiment::GraphTest13()
{
  Experiment::Config_Inst();
  Experiment::Config_Events();
  auto p = Genome_NOP();
  emp::vector<size_t> a{0,1,2}; emp::vector<size_t> b{3,4,5};emp::vector<size_t> c{6,7,8};
  emp::vector<emp::vector<size_t>> cohorts;
  cohorts.push_back(a); cohorts.push_back(b); cohorts.push_back(c); 
  mGoodWorld->Inject(p, 9);

  for(size_t i = 0; i < mGoodWorld->GetSize(); ++i)
  {
    auto & org = mGoodWorld->GetOrg(i);
    org.mBest = mRng->GetUInt(50, 100);
    org.mWorst = mRng->GetUInt(0, 49);
  }

  mGoodWorld->GetOrg(0).mBest = 100;
  mGoodWorld->GetOrg(1).mBest = 100;

  mGoodWorld->GetOrg(0).mWorst = 50;
  mGoodWorld->GetOrg(1).mWorst = 50;

  for(size_t i = 0; i < mGoodWorld->GetSize(); ++i)
  {
    std::cout << "ID: " << i << std::endl;
    std::cout << "mGood: " << mGoodWorld->GetOrg(i).mBest << std::endl;
    std::cout << "mWorst: " << mGoodWorld->GetOrg(i).mWorst << std::endl;
    std::cout << std::endl;
  }

  std::function<double(Agent&)> f = [](Agent & agent)
    {
      std::cout << "mBest" << std::endl;
      return agent.mBest;
    };
  std::function<double(Agent&)> g = [](Agent & agent)
    {
      std::cout << "mWorst" << std::endl;
      return agent.mWorst;
    };
  mGFitFuns.clear(); mGFitFuns.push_back(f); mGFitFuns.push_back(g);

  emp::CohortHostSelect(*mGoodWorld, cohorts, mGFitFuns, 3, 3, mRng);
}

#endif