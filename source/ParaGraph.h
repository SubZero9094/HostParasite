#ifndef HP_PARAGRAPH_H
#define HP_PARAGRAPH_H

#include <iostream>
#include <set>


#include "hp_config.h"
#include "../../Empirical/source/tools/Random.h"
#include "../../Empirical/source/tools/random_utils.h"
#include "../../Empirical/source/hardware/EventDrivenGP.h"
#include "../../Empirical/source/base/assert.h"

struct Node;

/* CONSTEXPR FOR HARDWARE */

//Number of bits in the tag per hardware
constexpr size_t _TAG_WIDTH_ = 16; 

/* NEW TYPE DECLARATIONS FOR HARDWARE*/

//Type for a actural hardware
using hardware_t = emp::EventDrivenGP_AW<_TAG_WIDTH_>;
//Type for the hardware genome
using program_t = hardware_t::Program;
//Instruction object for hardware
using ins_t = hardware_t::inst_t;
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

/* NEW TYPE DECLARATIONS FOR SIMPLICITY*/
using graph_t = std::vector<std::vector<emp::Ptr< Node >>>;
using coor_t = std::pair<size_t, size_t>;
using vote_t = std::unordered_map<double, double>;
using nodes_t = std::vector<emp::Ptr<Node>>;
using randnum_t = std::vector<size_t>;

/* STRUCT USED TO IMITATE A Node WITHIN THE SYSTEM */
struct Node 
{
  //Hardware within the Nodes
  emp::Ptr<hardware_t> mHW;
  //All of the Nodes that this Nodes has access too
  std::set<coor_t> mFriends;
  //Friend or Foe
  //0 -> friend
  //1 -> foe
  //mHW->TraitVector[FOE] = bool

  Node(emp::Ptr<inst_lib_t> ilib, emp::Ptr<event_lib_t> elib, emp::Ptr<emp::Random> rng)
  {
    mHW = emp::NewPtr<hardware_t>(ilib, elib, rng);
  }
};

class ParaGraph
{
    public:
        ParaGraph(const HPConfig & config, emp::Ptr<emp::Random> rng = nullptr) :
        GRA_DIM(config.GRA_DIM()), NUM_ITER(config.NUM_ITER()), NUM_FRI(config.NUM_FRI()),
        NUM_ENE(config.NUM_ENE()), mRng(rng), RNG_SEED(config.RNG_SEED()), UID(config.UID()),
        VOTE(config.VOTE()), POSX(config.POSX()), POSY(config.POSY()), FOE(config.FOE()),
        MAX_BND(config.MAX_BND()), MIN_BND(config.MIN_BND()), MIN_BIN_THSH(config.MIN_BIN_THSH()),
        MAX_CORES(config.MAX_CORES())
        {
        }

        ~ParaGraph()
        {
            for(size_t i = 0; i < mGraph.size(); ++i)
            {
                for(size_t j = 0; j < mGraph[i].size(); ++j)
                {
                mGraph[i][j]->mHW.Delete();
                mGraph[i][j].Delete();
                }
            }
            
            mNodes.clear();
        }

        /* FUNCTIONS DEDICATED TO THE STRUCTURE */

        //Will create the graph in a matrix and linear format
        void CreateGraph(size_t dim, size_t type, emp::Ptr<inst_lib_t> lib, emp::Ptr<event_lib_t> eve);

        //Will create the adjacency list for the nodes
        void CreateAdj(size_t dim, size_t type);

        //Will create the enemeies in the system
        void CreatePara(size_t ene);

        //Will reset all nodes to be friens
        void AllFri();

        //Will look to see if number is inside the UIDS
        bool Find(size_t x);

        //Will load UID into the hardware
        void ConfigureTraits();

        //Will set votes for all hardware in the graph
        void SetVotes(size_t vote);

        //Will set a vote for a single hardware
        void SetVote(size_t x, size_t y, size_t vote);

        //Will make the dictonary of votes
        void MakeFinalVotes();

        //Will check to see if they are in consensus
        double Consensus();

        //Largest legal vote
        double LargestLegal();

        //All legal votes
        double LegalVotes();

        //Load Genome to all with no parasite
        void GenomeNoPara(program_t p);

        //Load Genome to all with parasites enabled
        void GenomePara(program_t good, program_t bad);

        
        /* FUNCTIONS DEDICATED TO RunNoParaNING EXPERIMENT */

        //Run the experiment with no parasite
        double RunNoPara(size_t iter);

        //Run the experiment with parasite
        double RunPara(size_t iter);

        //Reset the graph
        void Reset();

        /* FUNCTIONS DEDICATED TO BE SETTER */

        //Are there parasites this experiment
        void SetPara(bool x) {mPara = x;}

    

        /* FUNCTIONS DEDICATED TO BE GETTER */
        
        //Return the graph
        graph_t GetGraph() const {return mGraph;}

        //Return mFinalVotes
        vote_t GetFinalVotes() const {return mFinalVotes;}

        //Return mRandNums
        randnum_t GetRandNums() const {return mRandomNums;}

        //Return all of a nodes friends
        std::set<coor_t> GetNodeNeig(size_t x, size_t y) { return mGraph[x][y]->mFriends;}

        //Return a specific node
        Node* GetNode(size_t x, size_t y) {return mGraph[x][y];}

        /* FUNCTIONS DEDICATED TO PRINTING OUT CRAP */

        //Will print out the friends list per node
        void PrintFriends();

        //Will print out who is a friend and who is a foe
        void PrintFOE();

        //Will print uid and vote
        void PrintUV();

        //Will print out the traits of each hardware
        void PrintTraits();

        //Will print the random numbers generated
        void PrintRanNums();

        //Will print the vote in mFinalVote
        void PrintFinVote();

        //Will print out the genome and type of node
        void PrintGenomes();


    private:
        /* CLASS SPECIFIC VARIABLES */

        //Holder of Node objects in matrix form
        graph_t mGraph;
        //Will hold the random numbers generated 
        randnum_t mRandomNums;
        //Will hold all the Nodes in a linear fashion
        nodes_t mNodes;
        //Will hold the final votes 
        vote_t mFinalVotes;
        //The scheduler
        emp::vector<coor_t> mSchedule;
        //The holder of parasite positions
        std::vector<coor_t> mParaCoors;

        /* GRAPH SPECIFIC PARAMATERS */

        //Dimension of the structure
        size_t GRA_DIM;
        //Number of iterations
        size_t NUM_ITER;
        //Number of Friends
        size_t NUM_FRI;
        //Number of enemies
        size_t NUM_ENE;
        //Random number generator
        emp::Ptr<emp::Random> mRng;
        //Random Number seed
        size_t RNG_SEED;
        //Do we have a parasite in the system?
        bool mPara;

        /* HARDWARE SPECIFIC PARAMATERS */

        //Position of UID within hw trait vector
        size_t UID;
        //Position of VOTE within hw trait vector
        size_t VOTE;
        //Position of X coordinate within hw trait vector
        size_t POSX;
        //Position of Y coordinate within hw trait vector
        size_t POSY;
        //Position of kind of node it is
        size_t FOE;
        //Uperbound on random numbers
        size_t MAX_BND;
        //Lowerbound on ranodm numbers
        size_t MIN_BND;
        //Minimum threshold
        double MIN_BIN_THSH;
        //Max number of cores a hardware can spawn
        size_t MAX_CORES;

};

/* FUNCTIONS DEDICATED TO RUNNING EXPERIMENT */

//Run the experiment
double ParaGraph::RunNoPara(size_t iter)
{
    if(mPara == true)
    {
        std::cout << "mPara == true: RunNoPara()" << std::endl;
        exit(-1);
    }

    double score = 0;
    for(size_t i = 0; i < iter; ++i)
    {
        emp::Shuffle(*mRng, mSchedule);
        for(auto p : mSchedule)
        {
            mGraph[p.first][p.second]->mHW->SingleProcess();
        }
        ParaGraph::MakeFinalVotes();
        score += ParaGraph::Consensus();
    }
    ParaGraph::MakeFinalVotes();
    score += ParaGraph::LargestLegal();
    score += ParaGraph::LegalVotes();

    return score;
}

//Run the experiment with parasite
double ParaGraph::RunPara(size_t iter)
{
    if(mPara == false)
    {
        std::cout << "mPara == false: RunPara" << std::endl;
        exit(-1);
    }
    
    double score = 0;
    for(size_t i = 0; i < iter; ++i)
    {
        emp::Shuffle(*mRng, mSchedule);
        for(auto p : mSchedule)
        {
            mGraph[p.first][p.second]->mHW->SingleProcess();
        }
        ParaGraph::MakeFinalVotes();
        score += ParaGraph::Consensus();
    }
    ParaGraph::MakeFinalVotes();
    score += ParaGraph::LegalVotes();
    score += ParaGraph::LargestLegal();

    return score;
}

//Reset the graph
void ParaGraph::Reset()
{
    mRandomNums.clear();
    mFinalVotes.clear();
    mParaCoors.clear();
    ParaGraph::ConfigureTraits();

    for(auto n : mNodes)
    {
        n->mHW->ResetHardware();
        n->mHW->SpawnCore(0, memory_t(), true);
    }
}

/* FUNCTIONS DEDICATED TO THE STRUCTURE */

//Will create the graph and set the traits of the hardware
//While also setting some hardware hyperparamaters
void ParaGraph::CreateGraph(size_t dim, size_t type, emp::Ptr<inst_lib_t> lib, emp::Ptr<event_lib_t> eve)
{
    std::cout << "CREATING GRAPH OF TYPE: " << type << std::endl;
    if (type == 0)
    {
        for(size_t i = 0; i < dim; ++i)
        {
            std::vector<emp::Ptr< Node >> temp;
            mGraph.push_back(temp);
            for(size_t j = 0; j < dim; ++j)
            {
                emp::Ptr<Node> n = emp::NewPtr<Node>(lib, eve, mRng);
                n->mHW->SetMinBindThresh(MIN_BIN_THSH);
                n->mHW->SetTrait(POSX, i);
                n->mHW->SetTrait(POSY, j);
                n->mHW->SetTrait(FOE, 0);
                n->mHW->SetMaxCores(MAX_CORES);
                mNodes.push_back(n);
                mGraph[i].push_back(n);
                mSchedule.emplace_back(std::make_pair(i,j));
            }
        }

        if(mNodes.size() != (GRA_DIM * GRA_DIM))
        {
            std::cout << "mNodes not the correct size!" << std::endl;
            exit(0);
        }
    }

    std::cout << "GRAPH CREATED! \n" << std::endl;
}

//Will create the adjacency list for the nodes
void ParaGraph::CreateAdj(size_t dim, size_t type)
{
    std::cout << "CREATING ADJ LIST TYPE: " << type << std::endl;
    if(type == 0)
    {
        for(size_t i = 0; i < mGraph.size(); ++i)
        {
            for(size_t j = 0; j < mGraph[i].size(); ++j)
            {
                //Right
                mGraph[i][j]->mFriends.emplace(std::make_pair((((i+dim)+1) % dim), j));
                //Left
                mGraph[i][j]->mFriends.emplace(std::make_pair((((i+dim)-1) % dim), j));
                //Up
                mGraph[i][j]->mFriends.emplace(std::make_pair(i, (((j+dim)+1) % dim)));
                //Down
                mGraph[i][j]->mFriends.emplace(std::make_pair(i, (((j+dim)-1) % dim)));
            }
        }
    }
    std::cout << "ADJ LIST CREATED!\n" << std::endl;
}

//Will create the parasites in the system
void ParaGraph::CreatePara(size_t ene)
{
    mParaCoors.clear();
    // std::cout << "CREATING ENEMIES LIST" << std::endl;
    auto enemies = emp::Choose(*mRng, mSchedule.size(), ene);

    for(auto i : enemies)
    {
        auto p = mSchedule[i];
        mGraph[p.first][p.second]->mHW->SetTrait(FOE, 1);
        mParaCoors.push_back(p);
    }

    

    if(ene != mParaCoors.size())
    {
        std::cout << "NOT THE CORRECT NUMBER OF ENEMIES" << std::endl;
        exit(0);
    }

    // std::cout << "DONE CREATING ENEMIES LIST! \n" << std::endl;
}

//Will reset all nodes to be friens
void ParaGraph::AllFri()
{
    for(auto n : mNodes)
        n->mHW->SetTrait(FOE, false);
}

//Will look to see if number is inside the UIDS
bool ParaGraph::Find(size_t x)
{
    for(auto y : mRandomNums)
    {
        if(x == y)
        {
            return true;
        }
    }
    return false;
}

//Will load new UIDs into the hardware and there is a check to see if they are the 
//same ammount of UIDs as nodes. Also will set vote to -33 and foe to 0.
void ParaGraph::ConfigureTraits()
{
    mRandomNums.clear();

    while(mRandomNums.size() != (GRA_DIM * GRA_DIM))
    {
        size_t num = mRng->GetUInt(MIN_BND, MAX_BND);

        if(!ParaGraph::Find(num))
        {
            mRandomNums.push_back(num);
        }
    }

    if(mRandomNums.size() != mNodes.size())
    {
        std::cout << "mRandomNums and mNodes are not the same size!" << std::endl;
        exit(0);
    }

    for(size_t i = 0; i < mRandomNums.size(); ++i)
    {
        mNodes[i]->mHW->SetTrait(UID, mRandomNums[i]);
        mNodes[i]->mHW->SetTrait(VOTE, -333);
        mNodes[i]->mHW->SetTrait(FOE, 0);

    }
}

//Will set votes for all hardware in the graph
void ParaGraph::SetVotes(size_t vote)
{
    for(auto n : mNodes)
    {
        n->mHW->SetTrait(VOTE, vote);
    }
}

//Will set a vote for a single hardware
void ParaGraph::SetVote(size_t x, size_t y, size_t vote)
{
    mGraph[x][y]->mHW->SetTrait(VOTE, vote);
}

//Will make the dictonary of votes for only friendly nodes in the system
void ParaGraph::MakeFinalVotes()
{
    mFinalVotes.clear();

    for(auto n : mNodes)
    {
        size_t vote = n->mHW->GetTrait(VOTE);
        size_t foe =  n->mHW->GetTrait(FOE);
        if(ParaGraph::Find(vote) && (foe == 0))
        {
            if(mFinalVotes.find(vote) == mFinalVotes.end())
            {
                mFinalVotes[vote] = 1;
            }

            else
            {
                ++mFinalVotes[vote];
            }
        }
    }
}

//Will check to see if they are in consensus
//Special case because all non-parasite nodes must be in consensus
double ParaGraph::Consensus()
{
    size_t cnt = 0;
    for(auto p : mFinalVotes)
    {
        if(p.second > cnt)
        {
            cnt = p.second;
        }
    }

    //We have a parasite in the system and they have consesnsus
    if((cnt == NUM_FRI) && (mPara == true))
    {
        return NUM_FRI;
    }
    //We do not have a paraiste and we have reached a consensus
    else if ((cnt == (NUM_FRI + NUM_ENE)) && (mPara == false))
    {
        return (NUM_FRI + NUM_ENE);
    }
    //consensus has not been reached
    else
    {
        return 0.0;
    }
}

//Largest legal vote
double ParaGraph::LargestLegal()
{
    size_t cnt = 0;
    for(auto p : mFinalVotes)
    {
        if(p.second > cnt)
        {
            cnt = p.second;
        }
    }

    return cnt;
}

//All legal votes
double ParaGraph::LegalVotes()
{
    size_t sum = 0;
    for(auto p : mFinalVotes)
    {
        sum += p.second;
    }

    return sum;
}

//Load Genome to all with no parasite
void ParaGraph::GenomeNoPara(program_t p)
{
    if(mPara == true)
    {
        std::cout << "mPara -> TRUE" << std::endl;
        exit(-1);
    }

    if(mParaCoors.size() > 0)
    {
        std::cout << "MORE THAN ONE PARASITE -> GenomePara()" << std::endl;
        exit(-1);
    }

    for(auto n : mNodes)
    {
        n->mHW->SetProgram(p);
    }
}

//Load Genome to all with parasites enabled
void ParaGraph::GenomePara(program_t good, program_t bad)
{
    if(mPara == false)
    {
        std::cout << "mPara -> FALSE" << std::endl;
        exit(-1);
    }

    if(mParaCoors.size() == 0)
    {
        std::cout << "NO PARASITE -> GenomePara()" << std::endl;
        exit(-1);
    }

    for(auto n : mNodes)
    {
        auto type = n->mHW->GetTrait(FOE);
        if(type == 0)
        {
            n->mHW->SetProgram(good);
        }
        else
        {
            n->mHW->SetProgram(bad);
        }
    }
}


/* FUNCTIONS DEDICATED TO PRINTING OUT CRAP */

//Will print friends list per node
void ParaGraph::PrintFriends()
{
    for(size_t i = 0; i < mGraph.size(); ++i)
    {
        for(size_t j = 0; j < mGraph[i].size(); ++j)
        {
            std::cout << "(" << i << "," << j << "): ";
            for(auto p : mGraph[i][j]->mFriends)
            {
                std::cout << "(" << p.first << "," << p.second << "), ";
            }
            std::cout << std::endl;
        }
    }
}

//Will print out kind of node everyone is
void ParaGraph::PrintFOE()
{
    for(size_t i = 0; i < mGraph.size(); ++i)
    {
        for(size_t j = 0; j < mGraph[i].size(); ++j)
        {
            std::cout << "(" << i << ", " << j << "): " << mGraph[i][j]->mHW->GetTrait(FOE) << std::endl;
        }
    }
}

//Will print the traits
void ParaGraph::PrintUV()
{
    for(size_t i = 0; i < mGraph.size(); ++i)
    {
        for(size_t j = 0; j < mGraph[i].size(); ++j)
        {
            std::cout << "(" << i << ", " << j << "): " << std::endl;
            std::cout << "UID: " << mGraph[i][j]->mHW->GetTrait(UID) << std::endl;
            std::cout << "VOTE: " << mGraph[i][j]->mHW->GetTrait(VOTE) << std::endl;
            std::cout << "FOE: " << mGraph[i][j]->mHW->GetTrait(FOE) << std::endl;
        }
    }
    std::cout << std::endl;
}

//Will print the traits
void ParaGraph::PrintTraits()
{
    for(size_t i = 0; i < mGraph.size(); ++i)
    {
        for(size_t j = 0; j < mGraph[i].size(); ++j)
        {
            std::cout << "(" << i << ", " << j << "): " << std::endl;
            std::cout << "UID: " << mGraph[i][j]->mHW->GetTrait(UID) << std::endl;
            std::cout << "VOTE: " << mGraph[i][j]->mHW->GetTrait(VOTE) << std::endl;
            std::cout << "POSX: " <<  mGraph[i][j]->mHW->GetTrait(POSX) << std::endl;
            std::cout << "POSY: " << mGraph[i][j]->mHW->GetTrait(POSY) << std::endl;
            std::cout << std::endl;
        }
    }
}

//Will print the random numbers genrated
void ParaGraph::PrintRanNums()
{
    std::cout << "mRanNums: ";
    for(auto x : mRandomNums)
    {
        std::cout << x << ", ";
    }
    std::cout << std::endl;
}

void ParaGraph::PrintFinVote()
{
    for(auto p : mFinalVotes)
    {
        std::cout << "Vote: " << p.first << " Cnt: " << p.second << std::endl; 
    }
    std::cout << std::endl;
}

//Will print out kind of node everyone is
void ParaGraph::PrintGenomes()
{
    for(size_t i = 0; i < mGraph.size(); ++i)
    {
        for(size_t j = 0; j < mGraph[i].size(); ++j)
        {
            std::cout << "(" << i << ", " << j << "): " << mGraph[i][j]->mHW->GetTrait(FOE) << std::endl;
            mGraph[i][j]->mHW->PrintProgramFull();
            std::cout << std::endl;
        }
    }
}
#endif