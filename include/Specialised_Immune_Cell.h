/*  Specialised_Immune_Cell.h  */
#ifndef SPECIALISED_IMMUNE_CELL
#define SPECIALISED_IMMUNE_CELL

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
*   The Speciallised Immune Cell Agent Class
**********************/
class SpecialisedImmuneCellAgent: public VirusCellInteractionAgents
{
public:

    // The immune cell states enum.
    enum SpecialisedImmuneCellStates{Healthy, Dead};

public:
    // Constructors
    SpecialisedImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, double theInfectedCellRecognitionProb, double theInfectedCellEliminationProb, double theSpecialisedImmuneCellRecruitRateOfSpecCell);
    
    SpecialisedImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, SpecialisedImmuneCellStates theState, double theInfectedCellRecognitionProb, double theInfectedCellEliminationProb, 
                                double theSpecialisedImmuneCellRecruitRateOfSpecCell, int theCountOfSpecCellsToRecruit, double theSpecCellsRecruitRemainder);

    // Destructor
    ~SpecialisedImmuneCellAgent();

    /* Required Getters */
    virtual repast::AgentId& getId(){                       return agentId;    }
    virtual const repast::AgentId& getId() const {          return agentId;    }

    /* Getters specific to the specialised immune cell agents */
    int getCellState(){                                     return specialisedImmuneCellState;  }
    double getInfCellRecognitionProb(){                     return infectedCellRecognitionProb; }
    double getInfCellEliminationProb(){                     return infectedCellEliminationProb; }
    double getSpecialisedImmuneCellRecruitRateOfSpecCell(){ return specialisedImmuneCellRecruitRateOfSpecCell;}
    int getCountOfSpecCellsToRecruit(){                     return countOfSpecCellsToRecruit;}
    double getSpecCellsRecruitRemainder(){                  return specCellsRecruitRemainder;}


    // Setter which sets all state variables and parameters of the agents. 
    // This setter is used only for updating agent copies at the buffer zone. It ensures that the non-local agents copies are always up-to-date with the original.
    void set(int currentRank, double newLifespan, double newAge, SpecialisedImmuneCellStates newSpecialisedImmuneCellState, double newInfectedCellRecognitionProb, 
        double newInfectedCellEliminationProb, double newSpecialisedImmuneCellRecruitRateOfSpecCell, double newCountOfSpecCellsToRecruit, double newSpecCellsRecruitRemainder);

    // The function triggered on each step which makes the agent act.
    void doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);

private:
    void move(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    void specialisedImmuneResponse(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    void specialisedCellRecruitingImmuneCells();

private:
    SpecialisedImmuneCellStates specialisedImmuneCellState;

    // Probability of detecting an infected epithelial cell.
    double infectedCellRecognitionProb;

    // Probability of eliminating an infected epithelial cell.
    double infectedCellEliminationProb;

    // The count of new specialised immune cell agents which are recruited by a specialised immune cell per discovered infection.
    double specialisedImmuneCellRecruitRateOfSpecCell;

    // The count of new specialised cells that this immune cell needs to recruit at this timestep.
    int countOfSpecCellsToRecruit;

    // The fractional remainder of specialised immune cells from the recruit rate (any decimal - if the recruit rate is 1.3, 
    // we'll keep 0.3 as remainder, so when it sums over time it will result in additional specialised immune cells to be recruited)
    double specCellsRecruitRemainder;
    
};

#endif // SPECIALISED_IMMUNE_CELL