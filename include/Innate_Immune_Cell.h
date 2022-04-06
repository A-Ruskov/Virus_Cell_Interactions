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

// Include the file containing the parent class.
#include "Virus_Cell_Agent.h"


/**********************
*   The Innate Immune Cell Agenrt Class
**********************/
class InnateImmuneCellAgent: public VirusCellInteractionAgents
{
public:
    enum InnateImmuneCellStates{Healthy, Dead};

public:
    // Constructors
    InnateImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, double theInfectedCellRecognitionProb, double theInfectedCellEliminationProb, 
                            double theSpecialisedImmuneCellRecruitProb, double theInnateImmuneCellRecruitRateOfInnateCell, double theSpecialisedImmuneCellRecruitRateOfInnateCell);

    InnateImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, InnateImmuneCellStates theState, double theInfectedCellRecognitionProb, 
                            double theInfectedCellEliminationProb, double theSpecialisedImmuneCellRecruitProb,  
                            double theInnateImmuneCellRecruitRateOfInnateCell, double theSpecialisedImmuneCellRecruitRateOfInnateCell,
                            double theCountOfInnateCellsToRecruit, double theInnateCellsRecruitRemainder,
                            double theCountOfSpecCellsToRecruit, double theSpecCellsRecruitRemainder);

    // Destructor
    ~InnateImmuneCellAgent();

    /* Required Getters */
    virtual repast::AgentId& getId(){                   return agentId;    }
    virtual const repast::AgentId& getId() const {      return agentId;    }

    /* Getters specific to the innate immune cell agents */
    int getCellState(){                             return innateImmuneCellState;  }
    double getInfCellRecognitionProb(){             return infectedCellRecognitionProb; }
    double getInfCellEliminationProb(){             return infectedCellEliminationProb; }
    double getSpecImmuneCellRecruitProb(){          return specialisedImmuneCellRecruitProb; }

    double getInnateCellRecruitRate(){              return innateImmuneCellRecruitRateOfInnateCell; }
    int getCountOfInnateCellsToRecruit(){           return countOfInnateCellsToRecruit; }
    double getInnateCellsRecruitRemainder(){        return innateCellsRecruitRemainder; }

    double getSpecialisedCellRecruitRate(){         return specialisedImmuneCellRecruitRateOfInnateCell; }
    int getCountOfSpecialisedCellsToRecruit(){      return countOfSpecCellsToRecruit; }
    double getSpecialisedCellsRecruitRemainder(){   return specCellsRecruitRemainder; }


    // Setter which sets all state variables and parameters of the agents. 
    // This setter is used only for updating agent copies at the buffer zone. It ensures that the non-local agents copies are always up-to-date with the original.
    void set(int currentRank, double newLifespan, double newAge, InnateImmuneCellStates newInnateImmuneCellState,  double newInfectedCellRecognitionProb, 
                double newInfectedCellEliminationProb, double newSpecialisedImmuneCellRecruitProb,
                double newInnateImmuneCellRecruitRateOfInnateCell, double newSpecialisedImmuneCellRecruitRateOfInnateCell,
                double newCountOfInnateCellsToRecruit, double newInnateCellsRecruitRemainder,
                double newCountOfSpecCellsToRecruit, double newSpecCellsRecruitRemainder);

    // The function triggered on each step which makes the agent act.
    void doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);

private:
    void move(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    void InnateImmuneResponse(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    void innateCellRecruitingImmuneCells( int specialisedImmCellsCountAtThisLocation );

private:
    InnateImmuneCellStates innateImmuneCellState;

    // Probability of detecting an infected epithelial cell.
    double infectedCellRecognitionProb;

    // Probability of eliminating an infected epithelial cell.
    double infectedCellEliminationProb;

    // The count of new innate immune cell agents which are recruited by an innate immune cell per discovered infection.
    double innateImmuneCellRecruitRateOfInnateCell;

    // The count of new innate cells that this immune cell needs to recruit at this timestep.
    int countOfInnateCellsToRecruit;

    // The fractional remainder of innate immune cells from the recruit rate (any decimal - if the recruit rate is 1.3, 
    // we'll keep 0.3 as remainder, so when it sums over time it will result in additional innate immune cells to be recruited)
    double innateCellsRecruitRemainder;


    /* Specialised cell recruitment variables */
    // Probability of an innate cell to recruit specialised immune cells when it detects an infection.
    double specialisedImmuneCellRecruitProb;

    // The count of new specialised immune cell agents which are recruited by an innate immune cell per discovered infection.
    double specialisedImmuneCellRecruitRateOfInnateCell;

    // The count of new specialised cells that this immune cell needs to recruit at this timestep.
    int countOfSpecCellsToRecruit;

    // The fractional remainder of specialised immune cells from the recruit rate (any decimal - if the recruit rate is 1.3, 
    // we'll keep 0.3 as remainder, so when it sums over time it will result in additional specialised immune cells to be recruited)
    double specCellsRecruitRemainder;
};


#endif // INNATE_IMMUNE_CELL