/* Virus_Cell_Agent.h */
#ifndef VIRUS_CELL_AGENT
#define VIRUS_CELL_AGENT

/**********************
*   Include files
**********************/
#include "repast_hpc/AgentId.h"
#include "repast_hpc/SharedContext.h"

// Spatial Projection includes
#include "repast_hpc/SharedDiscreteSpace.h"
#include "repast_hpc/GridComponents.h"


/* Base Agent Class */
class VirusCellInteractionAgents
{

protected:
    repast::AgentId     agentId;
    int              agentLifespan;
    int                 agentAge;

public:
    // Constructors
    VirusCellInteractionAgents(repast::AgentId theId);
    VirusCellInteractionAgents(repast::AgentId theId, double theLifespan, int theAge);

    // Destructor
    virtual ~VirusCellInteractionAgents();

    /* Required Getters */
    virtual repast::AgentId& getId(){                   return agentId;    }
    virtual const repast::AgentId& getId() const {      return agentId;    }

    /* Getters specific to this these specific agents */
    double getLifespan(){                                      return agentLifespan;      }
    double getAge(){                                  return agentAge;  }

    /* Setter */
    void set(int currentRank, double newLifespan, double newAge);

    virtual void doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    virtual void move(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace) {}
};



/* Serializable Agent Package */
struct VirusCellInteractionAgentPackage {
	
public:
    int    id;
    int    rank;
    int    type;
    int    currentRank;
    int lifespan;
    int age;
    int privateState;

    // Specific To Epithelial Cells
    int publicState;
    int infectedLifespan;
    int timeInfected;
    int divisionRate;
    int timeSinceLastDivision;
    double releaseDelay;
    double displayVirProteinsDelay;
    double extracellularReleaseProb;    // Probability of releasing a new virus particle in the extracellular space, when the cell starts producing the progeny virus
    double cellToCellTransmissionProb;     // Probability of directly infecting a neighbouring cell as a form of new virus particle release
    double virionReleaseRate;      // The virus count rate which a virus producing cell releases each hour.
    int countOfVirionsToRelease;     // The count of viruses that a virus producing cell will release at the current step.
    // The remainder of virus particles from the release rate (any decimal - if the release rate is 1.3, we'll keep 0.3 as remainder, so when it sums over time
    // it will result in additional virus particles to be released)
    double virionReleaseRemainder;

    // Specific to Virions
    double penetrationProbability;
    double clearanceProbability;
    double clearanceProbScaler;

    // Specific to the two immune cell types
    double infectedCellRecognitionProb;
    double infectedCellEliminationProb;

    // Specific to the innate immune cell agents
    double  specialisedImmuneCellRecruitProb;

    // Specific for handling buffer zone agent changes propagating back to original agents
    int     typeOfChangeToMakeToAgent;
    int     agentToEditID;
    int     agentToEditStartRank;
    int     agentToEditType;
    int     agentToEditCurrRank;
  

	
    /* Constructors */
    VirusCellInteractionAgentPackage(); // For serialization
    VirusCellInteractionAgentPackage(int _id, int _rank, int _type, int _currentRank, int _lifespan, int _age, int _privateState, int _publicState, int _infectedLifespan, 
                int _timeInfected, int _divisionRate, int _timeSinceLastDivision, double _releaseDelay, double _displayVirProteinsDelay, 
                double _extracellularReleaseProb, double _cellToCellTransmissionProb,
                double _virionReleaseRate, int _countOfVirionsToRelease, double _virionReleaseRemainder,
                double _penetrationProbability, double _clearanceProbability, double _clearanceProbScaler,
                double _infectedCellRecognitionProb, double _infectedCellEliminationProb,
                double _specialisedImmuneCellRecruitProb,
                int _typeOfChangeToMakeToAgent, int _agentToEditID, int _agentToEditStartRank, int _agentToEditType, int _agentToEditCurrRank);
	
    /* For archive packaging */
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version){
        ar & id;
        ar & rank;
        ar & type;
        ar & currentRank;
        ar & lifespan;
        ar & age;
        ar & privateState;

        // Specific to epithelial cell agents
        ar & publicState;
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

        // Specific to the innnate immune cell agents
        ar & specialisedImmuneCellRecruitProb;

        // Specific for handling buffer zone agent changes propagating back to original agents
        ar & typeOfChangeToMakeToAgent;
        ar & agentToEditID;
        ar & agentToEditStartRank;
        ar & agentToEditType;
        ar & agentToEditCurrRank;
    }

};
#endif // VIRUS_CELL_AGENT