/*  Innate_Immune_Cell.h  */
#ifndef INNATE_IMMUNE_CELL
#define INNATE_IMMUNE_CELL

/**********************
*   Include files
**********************/
#include "repast_hpc/AgentId.h"
#include "repast_hpc/SharedContext.h"

// Spatial Projection includes
#include "repast_hpc/SharedDiscreteSpace.h"

#include "Virus_Cell_Agent.h"


/**********************
*   The Innate Immune Cell Agenrt Class
**********************/
class InnateImmuneCellAgent: public VirusCellInteractionAgents
{
public:
    enum InnateImmuneCellStates{healthy, dead};

public:
    // Constructors
    InnateImmuneCellAgent(repast::AgentId theId);
    InnateImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge,   double theInfectedCellRecognitionProb, double theInfectedCellEliminationProb, double theSpecialisedImmuneCellRecruitProb);
    InnateImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, InnateImmuneCellStates theState, double theInfectedCellRecognitionProb, double theInfectedCellEliminationProb, double theSpecialisedImmuneCellRecruitProb);

    // Destructor
    ~InnateImmuneCellAgent();

    /* Required Getters */
    virtual repast::AgentId& getId(){                   return agentId;    }
    virtual const repast::AgentId& getId() const {      return agentId;    }

    /* Getters specific to this these specific agents */
    int getCellState(){                             return innateImmuneCellState;  }
    double getInfCellRecognitionProb(){             return infectedCellRecognitionProb; }
    double getInfCellEliminationProb(){             return infectedCellEliminationProb; }
    double getSpecImmuneCellRecruitProb(){          return specialisedImmuneCellRecruitProb; }

    bool isToRecruitNewSpecImmuneCell(){            return toRecruitNewSpecImmuneCell;  }
    bool isToRecruitNewInnateImmunceCell(){         return toRecruitNewInnateImmuneCell;}


    /* Setter */
    void set(int currentRank, double newLifespan, double newAge, InnateImmuneCellStates newInnateImmuneCellState,  double newInfectedCellRecognitionProb, double newInfectedCellEliminationProb, double newSpecialisedImmuneCellRecruitProb);

    void doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    // void move(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);

private:
    void move(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    void attemptToDetectInfectedCell(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);

private:
    InnateImmuneCellStates innateImmuneCellState;

    double infectedCellRecognitionProb;
    double infectedCellEliminationProb;
    double specialisedImmuneCellRecruitProb;

    bool   toRecruitNewSpecImmuneCell;
    bool   toRecruitNewInnateImmuneCell;
    
};

#endif // INNATE_IMMUNE_CELL