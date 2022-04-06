/*  Innate_Immune_Cell.cpp  */
// Class for implementing the Innate Immune Cell agents in the model.

/**********************
*   INCLUDE FILES
**********************/
#include "Innate_Immune_Cell.h"
#include "Virus_Cell_Agent.h"
#include "Epithelial_Cell_Agent.h"

#include "repast_hpc/initialize_random.h"
#include "repast_hpc/Point.h"
#include "repast_hpc/Moore2DGridQuery.h"


/**********************
*   InnateImmuneCellAgent::InnateImmuneCellAgent - Constructor for the InnateImmuneCellAgent class.
**********************/
InnateImmuneCellAgent::InnateImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, double theInfectedCellRecognitionProb, double theInfectedCellEliminationProb, 
                    double theSpecialisedImmuneCellRecruitProb, double theInnateImmuneCellRecruitRateOfInnateCell, double theSpecialisedImmuneCellRecruitRateOfInnateCell):
VirusCellInteractionAgents(theId, theLifespan, theAge),
innateImmuneCellState(Healthy),
infectedCellRecognitionProb(theInfectedCellRecognitionProb),
infectedCellEliminationProb(theInfectedCellEliminationProb),
specialisedImmuneCellRecruitProb(theSpecialisedImmuneCellRecruitProb),
innateImmuneCellRecruitRateOfInnateCell(theInnateImmuneCellRecruitRateOfInnateCell),
specialisedImmuneCellRecruitRateOfInnateCell(theSpecialisedImmuneCellRecruitRateOfInnateCell),
countOfInnateCellsToRecruit(0),
innateCellsRecruitRemainder(0.0),
countOfSpecCellsToRecruit(0),
specCellsRecruitRemainder(0.0)
{
}


/**********************
*   InnateImmuneCellAgent::InnateImmuneCellAgent - Constructor for the InnateImmuneCellAgent class.
**********************/
InnateImmuneCellAgent::InnateImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, InnateImmuneCellStates theState, double theInfectedCellRecognitionProb, 
                        double theInfectedCellEliminationProb, double theSpecialisedImmuneCellRecruitProb,  
                        double theInnateImmuneCellRecruitRateOfInnateCell, double theSpecialisedImmuneCellRecruitRateOfInnateCell,
                        double theCountOfInnateCellsToRecruit, double theInnateCellsRecruitRemainder,
                        double theCountOfSpecCellsToRecruit, double theSpecCellsRecruitRemainder):
VirusCellInteractionAgents(theId, theLifespan, theAge),
innateImmuneCellState(theState),
infectedCellRecognitionProb(theInfectedCellRecognitionProb),
infectedCellEliminationProb(theInfectedCellEliminationProb),
specialisedImmuneCellRecruitProb(theSpecialisedImmuneCellRecruitProb),
innateImmuneCellRecruitRateOfInnateCell(theInnateImmuneCellRecruitRateOfInnateCell),
specialisedImmuneCellRecruitRateOfInnateCell(theSpecialisedImmuneCellRecruitRateOfInnateCell),
countOfInnateCellsToRecruit(theCountOfInnateCellsToRecruit),
innateCellsRecruitRemainder(theInnateCellsRecruitRemainder),
countOfSpecCellsToRecruit(theCountOfSpecCellsToRecruit),
specCellsRecruitRemainder(theSpecCellsRecruitRemainder)
{
}




/**********************
*   InnateImmuneCellAgent::~InnateImmuneCellAgent - Destructor for the InnateImmuneCellAgent class.
**********************/
InnateImmuneCellAgent::~InnateImmuneCellAgent()
{ 

}



/**********************
*   InnateImmuneCellAgent::set - Setter for the agent. Sets all state variables and parameters of the agents. 
*   This setter is used only for updating agent copies at the buffer zone. It ensures that the non-local agents copies are always up-to-date with the original.
**********************/
void InnateImmuneCellAgent::set(int currentRank, double newLifespan, double newAge, InnateImmuneCellStates newInnateImmuneCellState,  double newInfectedCellRecognitionProb, 
                double newInfectedCellEliminationProb, double newSpecialisedImmuneCellRecruitProb,
                double newInnateImmuneCellRecruitRateOfInnateCell, double newSpecialisedImmuneCellRecruitRateOfInnateCell,
                double newCountOfInnateCellsToRecruit, double newInnateCellsRecruitRemainder,
                double newCountOfSpecCellsToRecruit, double newSpecCellsRecruitRemainder)
{
    agentId.currentRank(currentRank);
    agentLifespan = newLifespan;
    agentAge = newAge;
    innateImmuneCellState = newInnateImmuneCellState;
    infectedCellRecognitionProb = newInfectedCellRecognitionProb;
    infectedCellEliminationProb = newInfectedCellEliminationProb;
    specialisedImmuneCellRecruitProb = newSpecialisedImmuneCellRecruitProb;
    innateImmuneCellRecruitRateOfInnateCell = newInnateImmuneCellRecruitRateOfInnateCell;
    specialisedImmuneCellRecruitRateOfInnateCell = newSpecialisedImmuneCellRecruitRateOfInnateCell;
    countOfInnateCellsToRecruit = newCountOfInnateCellsToRecruit;
    innateCellsRecruitRemainder = newInnateCellsRecruitRemainder;
    countOfSpecCellsToRecruit = newCountOfSpecCellsToRecruit;
    specCellsRecruitRemainder = newSpecCellsRecruitRemainder;
}



/**********************
*   InnateImmuneCellAgent::doStep - Function for an agent to do a step. Will be triggered on every step,
*   The agent will then do an action depending on the surrounding agents and its internal state.
**********************/
void InnateImmuneCellAgent::doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    // Reset the counts of immune cells to be recruited in this step back to the default 0. We'll only recruit new cells if infection is detected.
    countOfInnateCellsToRecruit = 0;
    countOfSpecCellsToRecruit = 0;

    // Increment the agent's age and check if it has exceeded its lifespan. If that is the case, then set it to "Dead".
    ++agentAge;
    if(agentAge > agentLifespan && innateImmuneCellState != Dead)
    {
        innateImmuneCellState = Dead;
    }

    // If the agent is healthy, then it executes its InnateImmuneResponse submodel.
    if( innateImmuneCellState == Healthy )
    {
        InnateImmuneResponse(discreteGridSpace);
    }
}





/**********************
*   InnateImmuneCellAgent::move - Function for moving an agent to a neighbouring grid cell.
**********************/
void InnateImmuneCellAgent::move(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
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
*   InnateImmuneCellAgent::attemptToInfectCell - Gets the epithelial cell at the current grid position, attempts to find if it is infected and eliminate it.
**********************/
void InnateImmuneCellAgent::InnateImmuneResponse(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    std::vector<int> immuneCellLoc;
    discreteGridSpace->getLocation(agentId, immuneCellLoc);
    repast::Point<int> locationPoint(immuneCellLoc);
    
    // Get the agents at this grid cell
    repast::Moore2DGridQuery<VirusCellInteractionAgents> theQuery(discreteGridSpace);
    std::vector<VirusCellInteractionAgents*> agentsAtThisCell;
    theQuery.query(locationPoint, 0, true, agentsAtThisCell);

    std::vector<VirusCellInteractionAgents*>::iterator agentsAtThisCellIter;

    EpithelialCellAgent* theEpithelialCellToCheck = nullptr;
    int specialisedImmCellsCountHere = 0;


    for( agentsAtThisCellIter = agentsAtThisCell.begin(); agentsAtThisCellIter != agentsAtThisCell.end(); ++agentsAtThisCellIter )
    {
        // Find the epithelial cell which lives in this grid point
        if( (*agentsAtThisCellIter)->getId().agentType() == 0 )
        {
            theEpithelialCellToCheck = static_cast<EpithelialCellAgent*>(*agentsAtThisCellIter);
        }
        // Count how many specialised immune agents are situated at this grid location. This will then be used in the innateCellRecruitingImmuneCells (submodel)
        else if( (*agentsAtThisCellIter)->getId().agentType() == 3 )
        {
            specialisedImmCellsCountHere++;
        }

    }

    if( theEpithelialCellToCheck != nullptr )
    {
        // If the cell seems to be helathy. (NOTE: It could be infected but is not expressing any proteins, so there is no sign of infection ).
        if( theEpithelialCellToCheck->getExternalState() == EpithelialCellAgent::ExternalState::DisplayingViralProtein )
        {
            repast::DoubleUniformGenerator probabilityGenerator = repast::Random::instance()->createUniDoubleGenerator( 0.0 , 1.0 );
            double toDetectInfectedCell = probabilityGenerator.next();

            // Check if "the immune cell has actually managed to detect the viral proteins" (probability of detecting the infected cell)
            if( toDetectInfectedCell > 1-infectedCellRecognitionProb )
            {
                // If the agent has detected the infection, execute its SpecialisedCellRecruitingImmuneCells submodel for recruiting new immune cells.
                innateCellRecruitingImmuneCells( specialisedImmCellsCountHere );

                // Check if the immune cell manages to eliminate the virally-infected cell
                double toEliminateCell = probabilityGenerator.next();
                if( toEliminateCell > 1-infectedCellEliminationProb )
                {
                    // Eliminate the infected cell.
                    theEpithelialCellToCheck->eliminate();
                }
            }
        }
    }

    // Move the innate immune cell, then it could attepmt to detect another infected cell on the next time step.
    move(discreteGridSpace);
}



/**********************
*   InnateImmuneCellAgent::innateCellRecruitingImmuneCells - Function for "recruiting" new innate and specialised immune cells.
*   It calculates how many new cells of the two types need to be recruited. Then the Virus_Cell_Model class will create the required count of new agents.
**********************/
void InnateImmuneCellAgent::innateCellRecruitingImmuneCells( int specialisedImmCellsCountAtThisLocation )
{
    // If there are no specialised immune cells here, then recruit more innate cells. Innate cells act as first line of defense
    // They then go down in number when the specialised immune cells come in to kill the infection.
    if(specialisedImmCellsCountAtThisLocation == 0 )
    {
        // Calculate the count of innate immune cells which are to be recruited at this timestep and store the fractional remainder.
        countOfInnateCellsToRecruit = (int)(innateImmuneCellRecruitRateOfInnateCell + innateCellsRecruitRemainder);
        innateCellsRecruitRemainder = (innateImmuneCellRecruitRateOfInnateCell + innateCellsRecruitRemainder) - (double)countOfInnateCellsToRecruit;
    }


    repast::DoubleUniformGenerator newSpecImmuneCellRecruitProbGen = repast::Random::instance()->createUniDoubleGenerator(0.0, 1.0);
    double toRecruitSpecialisedCell = newSpecImmuneCellRecruitProbGen.next();

    // If the innate immune cell satisifies the probability for recruiting a specialised immune cell, then it can proceed to recruit them.
    // Only some innate cells are a part of the activation of the specialised (adaptive) immune response, so there is some probability of 
    // the cell being able to recruit specialised cells
    if( toRecruitSpecialisedCell > 1 - specialisedImmuneCellRecruitProb )
    {
        // Calculate the count of specialised immune cells which are to be recruited at this timestep and store the fractional remainder.
        countOfSpecCellsToRecruit = (int)(specialisedImmuneCellRecruitRateOfInnateCell + specCellsRecruitRemainder);
        specCellsRecruitRemainder = (specialisedImmuneCellRecruitRateOfInnateCell + specCellsRecruitRemainder) - countOfSpecCellsToRecruit;
    }
}