/*  Specialised_Immune_Cell.cpp  */
// Class for implementing the Specialised Immune Cell agents in the model.

/**********************
*   INCLUDE FILES
**********************/
#include "Specialised_Immune_Cell.h"
#include "Virus_Cell_Agent.h"
#include "Epithelial_Cell_Agent.h"

#include "repast_hpc/initialize_random.h"
#include "repast_hpc/Point.h"
#include "repast_hpc/Moore2DGridQuery.h"


/**********************
*   SpecialisedImmuneCellAgent::SpecialisedImmuneCellAgent - Constructor for the SpecialisedImmuneCellAgent class. 
**********************/
SpecialisedImmuneCellAgent::SpecialisedImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, double theInfectedCellRecognitionProb, double theInfectedCellEliminationProb, double theSpecialisedImmuneCellRecruitRateOfSpecCell):
VirusCellInteractionAgents(theId, theLifespan, theAge),
specialisedImmuneCellState(Healthy),
infectedCellRecognitionProb(theInfectedCellRecognitionProb),
infectedCellEliminationProb(theInfectedCellEliminationProb),
specialisedImmuneCellRecruitRateOfSpecCell(theSpecialisedImmuneCellRecruitRateOfSpecCell),
countOfSpecCellsToRecruit(0),
specCellsRecruitRemainder(0.0)
{
}


/**********************
*   SpecialisedImmuneCellAgent::SpecialisedImmuneCellAgent - Constructor for the SpecialisedImmuneCellAgent class.
**********************/
SpecialisedImmuneCellAgent::SpecialisedImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, SpecialisedImmuneCellStates theState, double theInfectedCellRecognitionProb, double theInfectedCellEliminationProb, 
                                double theSpecialisedImmuneCellRecruitRateOfSpecCell, int theCountOfSpecCellsToRecruit, double theSpecCellsRecruitRemainder):
VirusCellInteractionAgents(theId, theLifespan, theAge),
specialisedImmuneCellState(theState),
infectedCellRecognitionProb(theInfectedCellRecognitionProb),
infectedCellEliminationProb(theInfectedCellEliminationProb),
specialisedImmuneCellRecruitRateOfSpecCell(theSpecialisedImmuneCellRecruitRateOfSpecCell),
countOfSpecCellsToRecruit(theCountOfSpecCellsToRecruit),
specCellsRecruitRemainder(theSpecCellsRecruitRemainder)
{
}




/**********************
*   SpecialisedImmuneCellAgent::~SpecialisedImmuneCellAgent - Destructor for the SpecialisedImmuneCellAgent class.
**********************/
SpecialisedImmuneCellAgent::~SpecialisedImmuneCellAgent()
{ 

}



/**********************
*   SpecialisedImmuneCellAgent::set - Setter for the agent. Sets all state variables and parameters of the agents. 
*   This setter is used only for updating agent copies at the buffer zone. It ensures that the non-local agents copies are always up-to-date with the original.
**********************/
void SpecialisedImmuneCellAgent::set(int currentRank, double newLifespan, double newAge, SpecialisedImmuneCellStates newSpecialisedImmuneCellState, double newInfectedCellRecognitionProb, 
                                     double newInfectedCellEliminationProb, double newSpecialisedImmuneCellRecruitRateOfSpecCell, double newCountOfSpecCellsToRecruit, double newSpecCellsRecruitRemainder)
{
    agentId.currentRank(currentRank);
    agentLifespan = newLifespan;
    agentAge = newAge;
    specialisedImmuneCellState = newSpecialisedImmuneCellState;
    infectedCellRecognitionProb = newInfectedCellRecognitionProb;
    infectedCellEliminationProb = newInfectedCellEliminationProb;
    specialisedImmuneCellRecruitRateOfSpecCell = newSpecialisedImmuneCellRecruitRateOfSpecCell;
    countOfSpecCellsToRecruit = newCountOfSpecCellsToRecruit;
    specCellsRecruitRemainder = newSpecCellsRecruitRemainder;
}



/**********************
*   SpecialisedImmuneCellAgent::doStep - Function for an agent to do a step. Will be triggered on every step,
*   The agent will then do an action depending on the surrounding agents and its internal state.
**********************/
void SpecialisedImmuneCellAgent::doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    // Reset the count of specialised cells to recruit to the default 0.
    countOfSpecCellsToRecruit = 0;

    // Increment the agent's age and check if it has exceeded its lifespan. If that is the case, then set it to "Dead".
    ++agentAge;
    if(agentAge > agentLifespan && specialisedImmuneCellState != Dead)
    {
        specialisedImmuneCellState = Dead;
    }

    // If the agent is healthy, then it executes its SpecialisedImmuneResponse submodel.
    if( specialisedImmuneCellState == Healthy )
    {
        specialisedImmuneResponse(discreteGridSpace);
    }
}




/**********************
*   SpecialisedImmuneCellAgent::move - Function for moving an agent to a neighbouring grid cell.
**********************/
void SpecialisedImmuneCellAgent::move(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    if(discreteGridSpace == nullptr)
    {
        std::cout<<"An incorrect pointer to the grid projection was passed to InnateImmuneCellAgent::move! The move cannot be executed."<<std::endl;
        return;
    }

    // Get the current agent location from the grid projection.
    std::vector<int> currentLocation;
    discreteGridSpace->getLocation(agentId, currentLocation);
   
    repast::IntUniformGenerator moveGenerator = repast::Random::instance()->createUniIntGenerator(-1, 1);
    int moveX = 0;
    int moveY = 0;

    // Generate the agent's move stochasticly. The move is ensured to be to one of the 8 directly neighbouring cell.
    while(moveX == 0 && moveY == 0)
    {
        moveX = moveGenerator.next();
        moveY = moveGenerator.next();
    }

    // Create the new location coordinate.
    std::vector<int> newLocation;
    newLocation.push_back(currentLocation[0] + moveX);
    newLocation.push_back(currentLocation[1] + moveY);
    repast::Point<int> movePoint(newLocation);

    // Move the agent
    discreteGridSpace->moveTo(agentId, movePoint);
}



/**********************
*   SpecialisedImmuneCellAgent::specialisedImmuneResponse - Gets the epithelial cell at the current grid position, attempts to find if it is infected and eliminate it.
*   Executes the model's SpecialisedImmuneResponse submodel.
**********************/
void SpecialisedImmuneCellAgent::specialisedImmuneResponse(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    std::vector<int> immuneCellLoc;
    discreteGridSpace->getLocation(agentId, immuneCellLoc);
    repast::Point<int> locationPoint(immuneCellLoc);
    
    // Get the agents at this grid cell
    repast::Moore2DGridQuery<VirusCellInteractionAgents> theQuery(discreteGridSpace);
    std::vector<VirusCellInteractionAgents*> agentsAtThisCell;
    theQuery.query(locationPoint, 0, true, agentsAtThisCell);

    std::vector<VirusCellInteractionAgents*>::iterator agentsAtThisCellIter;
    for( agentsAtThisCellIter = agentsAtThisCell.begin(); agentsAtThisCellIter != agentsAtThisCell.end(); ++agentsAtThisCellIter )
    {
        // Find the epithelial cell which lives in this grid point
        if( (*agentsAtThisCellIter)->getId().agentType() == 0 )
        {
            EpithelialCellAgent* theEpithelialCell = static_cast<EpithelialCellAgent*>(*agentsAtThisCellIter);

            // If the cell seems to be helathy. (NOTE: It could be infected but is not expressing any proteins, so there is no sign of infection ).
            if( theEpithelialCell->getExternalState() == EpithelialCellAgent::ExternalState::DisplayingViralProtein )
            {
                repast::DoubleUniformGenerator probabilityGenerator = repast::Random::instance()->createUniDoubleGenerator( 0.0 , 1.0 );
                double toDetectInfectedCell = probabilityGenerator.next();

                // Check if "the immune cell has actually managed to detect the viral proteins" (probability of detecting the infected cell)
                if( toDetectInfectedCell > 1-infectedCellRecognitionProb )
                {
                    // If the agent has detected the infection, execute its SpecialisedCellRecruitingImmuneCells submodel for recruiting new immune cells.
                    specialisedCellRecruitingImmuneCells();

                    // Check if the immune cell manages to eliminate the virally-infected cell
                    double toEliminateCell = probabilityGenerator.next();
                    if( toEliminateCell > 1-infectedCellEliminationProb )
                    {
                        // Eliminate the infected cell.
                        theEpithelialCell->eliminate();
                    }
                }
                break;
            }
            break;
        }
    }

    // Move the specialised immune cell, then it could attepmt to detect another infected cell on the next time step.
    move(discreteGridSpace);
}



/**********************
*   SpecialisedImmuneCellAgent::specialisedCellRecruitingImmuneCells - Function for "recruiting" new specialised immune cells.
*   It calculates how many new cells need to be recruited. Then the Virus_Cell_Modell class will create the required count of new agents.
**********************/
void SpecialisedImmuneCellAgent::specialisedCellRecruitingImmuneCells()
{
    countOfSpecCellsToRecruit = int( specialisedImmuneCellRecruitRateOfSpecCell + specCellsRecruitRemainder );
    specCellsRecruitRemainder = ( specialisedImmuneCellRecruitRateOfSpecCell + specCellsRecruitRemainder ) - countOfSpecCellsToRecruit;
}