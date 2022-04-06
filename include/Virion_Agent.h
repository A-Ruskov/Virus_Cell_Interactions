/*  Virion_Agent.h  */
#ifndef VIRION_AGENT
#define VIRION_AGENT

/**********************
*   Include files
**********************/
#include "repast_hpc/AgentId.h"
#include "repast_hpc/SharedContext.h"

// Spatial Projection includes
#include "repast_hpc/SharedDiscreteSpace.h"

// Include the file containing the parent class.
#include "Virus_Cell_Agent.h"


/**********************
*   Forward class declarations
**********************/
class EpithelialCellAgent;


/**********************
*   The Virion Agenrt Class
**********************/
class VirionAgent: public VirusCellInteractionAgents
{
public:
    enum VirionStates{Free_Virion, Dead, Contained};

public:
    /* Constructors */
    VirionAgent(repast::AgentId theId, double theLifespan, int theAge, double thePenetrationProb, double theClearanceProbability, double theClearanceProbScaler);
    VirionAgent(repast::AgentId theId, double theLifespan, int theAge, VirionStates theState, double thePenetrationProb, double theClearanceProbability, double theClearanceProbScaler);

    // Destructor
    ~VirionAgent();

    /* Required Getters */
    virtual repast::AgentId& getId(){                   return agentId;    }
    virtual const repast::AgentId& getId() const {      return agentId;    }

    /* Getters specific to the virion agents */
    int getVirionState(){                             return virionState;  }
    double getPenetrationProb(){                      return penetrationProbability; }
    double getClearanceProb(){                        return clearanceProbability; }
    double getClearanceProbScaler(){                  return clearanceProbScaler; }


    // Setter which sets all state variables and parameters of the agents. 
    // This setter is used only for updating agent copies at the buffer zone. It ensures that the non-local agents copies are always up-to-date with the original.
    void set(int currentRank, double newLifespan, double newAge, VirionStates newVirState, double newPenetrationProb, 
            double newClearanceProbability, double newClearanceProbScaler);

    // The function triggered on each step which makes the agent act.
    void doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    
private:    
    void move(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    void virionClearance( int countOfImmuneAgentsAtVirionLocation );
    void attemptToInfectCell(EpithelialCellAgent* theEpithelialCell, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);

private:
    VirionStates virionState;

    // The probability to penetrate an epithelial cell.
    double       penetrationProbability;

    // The probability of clearing the virus through unmodelled immune mechanisms
    double       clearanceProbability;

    // The scaler used to scale the base clearance probability. That scaler is applied for each immune cell at the vicinity of the virus
    // The more immune cells there are, the higher the chance there are also other unmodelled immune mechanisms at the location of the agent
    double       clearanceProbScaler; 
};

#endif // VIRION_AGENT