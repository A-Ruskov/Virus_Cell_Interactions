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

#include "Virus_Cell_Agent.h"


/**********************
*   The Speciallised Immune Cell Agenrt Class
**********************/
class SpecialisedImmuneCellAgent: public VirusCellInteractionAgents
{
public:
    enum SpecialisedImmuneCellStates{healthy, dead};

public:
    // Constructors
    SpecialisedImmuneCellAgent(repast::AgentId theId);
    SpecialisedImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, double theInfectedCellRecognitionProb, double theinfectedCellEliminationProb);
    SpecialisedImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, SpecialisedImmuneCellStates theState, double theInfectedCellRecognitionProb, double theinfectedCellEliminationProb);

    // Destructor
    ~SpecialisedImmuneCellAgent();

    /* Required Getters */
    virtual repast::AgentId& getId(){                   return agentId;    }
    virtual const repast::AgentId& getId() const {      return agentId;    }

    /* Getters specific to this these specific agents */
    int getCellState(){                             return specialisedImmuneCellState;  }
    double getInfCellRecognitionProb(){             return infectedCellRecognitionProb; }
    double getInfCellEliminationProb(){             return infectedCellEliminationProb; }
    
    bool isToRecruitNewSpecImmuneCell(){            return toRecruitNewSpecImmuneCell;  }


    /* Setter */
    void set(int currentRank, double newLifespan, double newAge, SpecialisedImmuneCellStates newSpecialisedImmuneCellState, double newInfectedCellRecognitionProb, double newInfectedCellEliminationProb);

    void doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);

private:
    void move(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    void attemptToDetectInfectedCell(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);

private:
    SpecialisedImmuneCellStates specialisedImmuneCellState;

    double infectedCellRecognitionProb;
    double infectedCellEliminationProb;

    bool   toRecruitNewSpecImmuneCell;
    
};

#endif // SPECIALISED_IMMUNE_CELL