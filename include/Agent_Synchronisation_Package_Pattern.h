/****************************************
***** Agent_Sychronisation_Pattern_ *****
****************************************/

// Contains everything required for syncrhonising/moving agents across processes. 
// That includes the agent packages implementation, the package provider and receiver

/****************
* Include Files
****************/
// Include agent related files
#include "Virus_Cell_Agent.h"



/**********************
* Serializable Agent Package - A package containing all potentailly required variables for creating/updating an agent across processes.
* Used to move agents between processes and update and agent copies which are found in the bufferzone.
**********************/
struct VirusCellInteractionAgentPackage {
	
public:

    /*** Variables applicable to all agent types ***/
    int id;
    int rank;
    int type;
    int currentRank;
    double lifespan;
    int age;
    int internalState;

    /*** Specific To Epithelial Cell Agents ***/
    int externalState = -1;     // The state of the cell which other agents perceive.
    double infectedLifespan = -1;      // The amount of time a cell can live after getting infected.
    int timeInfected = -1;      // The amount of time the epithelial cell has been infected for
    double divisionRate = -1;      // The time that needs to pass between two divisions of an epithelial cell agent.
    int timeSinceLastDivision = -1;     //The time which has passed since the last division of the cell.
    double releaseDelay = -1.0;     // Time that needs to pass after a cell gets infected before it starts releasing new viruses in the grid/infecting neighbouring cells.
    double displayVirProteinsDelay = -1.0;    // Time that needs to pass after a cell gets infected before it starts appearing infected to the other agents.
    double extracellularReleaseProb = -1.0;    // Probability of releasing a new virus particle in the extracellular space, when the cell starts producing the progeny virus
    double cellToCellTransmissionProb = -1.0;     // Probability of directly infecting a neighbouring cell as a form of new virus particle release
    double virionReleaseRate = -1.0;      // The virus count rate which a virus producing cell releases each hour.
    int countOfVirionsToRelease = -1;     // The count of viruses that a virus producing cell will release at the current step.
    
    // The remainder of virus particles from the release rate (any decimal - if the release rate is 1.3, we'll keep 0.3 as remainder, so when it sums over time
    // it will result in additional virus particles to be released)
    double virionReleaseRemainder = -1.0;

    /*** Specific to Virions ***/
    double penetrationProbability = -1.0;   // Probability of a virus to penetrate a cell (to infect it)
    double clearanceProbability = -1.0;     // Probability of a virus getting cleared through unmodelled immune response mechanisms.
    double clearanceProbScaler = -1.0;      // Parameter for scaling the probability of a virus getting cleared for each immune cell at the virus location.

    /*** Specific to the two immune cell types. Applicable to both innate and specialised immune cell agents ***/
    double infectedCellRecognitionProb = -1.0;
    double infectedCellEliminationProb = -1.0;
    // The count of new specialised immune cell agents which are recruited by an immune cell per discovered infection.
    double  specialisedImmuneCellRecruitRate = -1.0;
    // The count of new specialised cells that the immune cell needs to recruit at this timestep.
    int countOfSpecCellsToRecruit = -1;
    // The fractional remainder of specialised immune cells from the recruit rate (any decimal - if the recruit rate is 1.3, 
    // we'll keep 0.3 as remainder, so when it sums over time it will result in additional specialised immune cells to be recruited)
    double specCellsRecruitRemainder = -1.0;

    /*** Specific to the innate immune cell agents ***/
    // Tre probability of an innate immune cell to recruit specialised immune cells, when they discover an invection
    double specialisedImmuneCellRecruitProb = -1.0;
    // The count of new innate immune cell agents which are recruited by an innate immune cell per discovered infection.
    double innateImmuneCellRecruitRate = -1.0;
    // The count of new innate cells that this immune cell needs to recruit at this timestep.
    int countOfInnateCellsToRecruit = -1;
    // The fractional remainder of innate immune cells from the recruit rate (any decimal - if the recruit rate is 1.3, 
    // we'll keep 0.3 as remainder, so when it sums over time it will result in additional innate immune cells to be recruited)
    double innateCellsRecruitRemainder = -1.0;

   
    /*** Specific for handling buffer zone agent changes propagating back to original agents ***/
    int     typeOfChangeToMakeToAgent = -1;
    int     agentToEditID = -1;
    int     agentToEditStartRank = -1;
    int     agentToEditType = -1;
    int     agentToEditCurrRank = -1;
  
	
    /* Constructors */
    VirusCellInteractionAgentPackage(); // For serialization

    // Constructor of the serializable package to be used for Epithelial Cell Agents
    VirusCellInteractionAgentPackage(int _id, int _rank, int _type, int _currentRank, double _lifespan, int _age, int _internalState, int _externalState, double _infectedLifespan, 
                int _timeInfected, double _divisionRate, int _timeSinceLastDivision, double _releaseDelay, double _displayVirProteinsDelay, 
                double _extracellularReleaseProb, double _cellToCellTransmissionProb,
                double _virionReleaseRate, int _countOfVirionsToRelease, double _virionReleaseRemainder,
                int _typeOfChangeToMakeToAgent, int _agentToEditID, int _agentToEditStartRank, int _agentToEditType, int _agentToEditCurrRank);

    // Constructor of the serializable package to be used for Virion Agents
    VirusCellInteractionAgentPackage(int _id, int _rank, int _type, int _currentRank, double _lifespan, int _age, int _internalState, 
                                    double _penetrationProbability, double _clearanceProbability, double _clearanceProbScaler);

    // Constructor of the serializable package to be used for Innate Immune Cell Agents
    VirusCellInteractionAgentPackage(int _id, int _rank, int _type, int _currentRank, double _lifespan, int _age, int _internalState, 
                                    double _infectedCellRecognitionProb, double _infectedCellEliminationProb,
                                    double _specialisedImmuneCellRecruitProb, double _specialisedImmuneCellRecruitRate,
                                    int _countOfSpecCellsToRecruit, double _specCellsRecruitRemainder,
                                    double _innateImmuneCellRecruitRate, double _countOfInnateCellsToRecruit,
                                    double _innateCellsRecruitRemainder);

    // Constructor of the serializable package to be used for Specialised Immune Cell Agents
    VirusCellInteractionAgentPackage(int _id, int _rank, int _type, int _currentRank, double _lifespan, int _age, int _internalState, 
                                    double _infectedCellRecognitionProb, double _infectedCellEliminationProb,
                                    double _specialisedImmuneCellRecruitRate, int _countOfSpecCellsToRecruit, double _specCellsRecruitRemainder);
	
    /* For archive packaging */
    template<class Archive>

    // Serialises the passed variables from the archive and stores them in the package variables.
    void serialize(Archive &ar, const unsigned int version){
        ar & id;
        ar & rank;
        ar & type;
        ar & currentRank;
        ar & lifespan;
        ar & age;
        ar & internalState;

        // Specific to epithelial cell agents
        ar & externalState;
        ar & infectedLifespan;
        ar & timeInfected;
        ar & divisionRate;
        ar & timeSinceLastDivision;
        ar & releaseDelay;
        ar & displayVirProteinsDelay;
        ar & extracellularReleaseProb;
        ar & cellToCellTransmissionProb;
        ar & virionReleaseRate;
        ar & countOfVirionsToRelease;
        ar & virionReleaseRemainder;

        // Specific to virion agents
        ar & penetrationProbability;
        ar & clearanceProbability;
        ar & clearanceProbScaler;

        // Specific to the two immune cell types
        ar & infectedCellRecognitionProb;
        ar & infectedCellEliminationProb;
        ar & specialisedImmuneCellRecruitRate;
        ar & countOfSpecCellsToRecruit;
        ar & specCellsRecruitRemainder;

        // Specific to the innnate immune cell agents
        ar & specialisedImmuneCellRecruitProb;
        ar & innateImmuneCellRecruitRate;
        ar & countOfInnateCellsToRecruit;
        ar & innateCellsRecruitRemainder;

        // Specific for handling buffer zone agent changes propagating back to original agents
        ar & typeOfChangeToMakeToAgent;
        ar & agentToEditID;
        ar & agentToEditStartRank;
        ar & agentToEditType;
        ar & agentToEditCurrRank;
    }
};





/**********************
*   Agent Package Provider class. Builds and provides the packages which contain the agents' state variables and parameters
**********************/
class VirusCellInteractionAgentsPackageProvider {
	
private:
    repast::SharedContext<VirusCellInteractionAgents>* agentsContext;
	
public:
	
    VirusCellInteractionAgentsPackageProvider(repast::SharedContext<VirusCellInteractionAgents>* contextPtr);
	
    void providePackage(VirusCellInteractionAgents * agent, std::vector<VirusCellInteractionAgentPackage>& out);
	
    void provideContent(repast::AgentRequest req, std::vector<VirusCellInteractionAgentPackage>& out);
	
};





/**********************
* Agent Package Receiver class. Receives the packages which contain the agents' state variables and parameters and updates/creates the agents in the appropriate process.
**********************/
class VirusCellInteractionAgentsPackageReceiver {
	
private:
    repast::SharedContext<VirusCellInteractionAgents>* agentsContext;
	
public:
	
    VirusCellInteractionAgentsPackageReceiver(repast::SharedContext<VirusCellInteractionAgents>* agentPtr);
	
    VirusCellInteractionAgents * createAgent(VirusCellInteractionAgentPackage package);
	
    void updateAgent(VirusCellInteractionAgentPackage package);
};