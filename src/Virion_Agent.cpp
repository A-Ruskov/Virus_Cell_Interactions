/*  Virion_Agent.cpp  */
// Class for implementing the Innate Immune Cell agents in the model.

/**********************
*   INCLUDE FILES
**********************/
#include "Virion_Agent.h"
#include "Virus_Cell_Agent.h"
#include "Epithelial_Cell_Agent.h"

#include "repast_hpc/initialize_random.h"
#include "repast_hpc/Point.h"
#include "repast_hpc/Moore2DGridQuery.h"


/**********************
*   VirionAgent::VirionAgent - Constructor for the VirionAgent class.
**********************/
VirionAgent::VirionAgent(repast::AgentId theId, double theLifespan, int theAge, double thePenetrationProb, double theClearanceProbability, double theClearanceProbScaler):
VirusCellInteractionAgents(theId, theLifespan, theAge),
virionState(Free_Virion),
penetrationProbability(thePenetrationProb),
clearanceProbability(theClearanceProbability),
clearanceProbScaler(theClearanceProbScaler)
{
}



/**********************
*   VirionAgent::VirionAgent - Constructor for the VirionAgent class.
**********************/
VirionAgent::VirionAgent(repast::AgentId theId, double theLifespan, int theAge, VirionStates theState, double thePenetrationProb, double theClearanceProbability, double theClearanceProbScaler):
VirusCellInteractionAgents(theId, theLifespan, theAge),
virionState(theState),
penetrationProbability(thePenetrationProb),
clearanceProbability(theClearanceProbability),
clearanceProbScaler(theClearanceProbScaler)
{
}



/**********************
*   VirionAgent::~VirionAgent - Destructor for the VirionAgent class.
**********************/
VirionAgent::~VirionAgent()
{ 
}



/**********************
*   VirionAgent::set - Setter for the agent. Sets all state variables and parameters of the agents. 
*   This setter is used only for updating agent copies at the buffer zone. It ensures that the non-local agents copies are always up-to-date with the original.
**********************/
void VirionAgent::set(int currentRank, double newLifespan, double newAge, VirionStates newVirState, double newPenetrationProb, double newClearanceProb, double newClearanceProbScaler)
{
    agentId.currentRank(currentRank);
    agentLifespan = newLifespan;
    agentAge = newAge;
    virionState = newVirState;
    penetrationProbability = newPenetrationProb;
    clearanceProbability = newClearanceProb;
    clearanceProbScaler = newClearanceProbScaler;
}



/**********************
*   VirionAgent::doStep - Function for an agent to do a step. Will be triggered on every step,
*   The agent will then do an action depending on the surrounding agents and its internal state.
**********************/
void VirionAgent::doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    // Increment the agent's age and check if it has exceeded its lifespan. If that is the case, then set it to "Dead".
    ++agentAge;
    if(agentAge > agentLifespan && virionState != Dead)
    {
        virionState = Dead;
        return;
    }
    
    if( virionState == Free_Virion )
    {
        std::vector<int> virionLocation;
        discreteGridSpace->getLocation(agentId, virionLocation);
        repast::Point<int> locationPoint(virionLocation);
    
        // Get the agents at this grid cell
        repast::Moore2DGridQuery<VirusCellInteractionAgents> theQuery(discreteGridSpace);
        std::vector<VirusCellInteractionAgents*> agentsAtThisCell;
        theQuery.query(locationPoint, 0, true, agentsAtThisCell);

        std::vector<VirusCellInteractionAgents*>::iterator agentsAtThisCellIter;

        EpithelialCellAgent* theEpithelialCellHere = nullptr;
        int immuneCellAgentsCount = 0;

        for( agentsAtThisCellIter = agentsAtThisCell.begin(); agentsAtThisCellIter != agentsAtThisCell.end(); ++agentsAtThisCellIter )
        {
            // Find the epithelial cell which lives in this grid point
            if( (*agentsAtThisCellIter)->getId().agentType() == 0 )
            {
               theEpithelialCellHere = static_cast<EpithelialCellAgent*>(*agentsAtThisCellIter);
            }
            // Count how many immune cell agents are in the same grid cell. This will then be used in the virionClearance function/submodel
            else if( (*agentsAtThisCellIter)->getId().agentType() == 2 || (*agentsAtThisCellIter)->getId().agentType() == 3 )
            {
                immuneCellAgentsCount++;
            }
        }

        // Execute the virion clearance submodel. /The virion might die due to unmodelled immune response mechanisms/
        virionClearance( immuneCellAgentsCount );

        // If the virus did not get cleared /it is not dead/, then attempt to infect a cell.
        if( virionState != Dead )
        {
            attemptToInfectCell( theEpithelialCellHere, discreteGridSpace );
        }
    }
}



/**********************
*   VirionAgent::virionClearance - Handles the possible death of a virion caused by unmodelled immune mechanisms. 
*   Based on a probability value the agent could get cleaned from the grid.
*   Implementation of the VirionClearance submodel.
**********************/
void VirionAgent::virionClearance( int countOfImmuneAgentsAtVirionLocation )
{
    // For each immune cell at the same grid position, the clearance probability should grow. The more modelled immune cells there are here,
    // The more the chance of other unmodelled immune mechanisms to happen there as well.
    repast::DoubleUniformGenerator clearanceGen = repast::Random::instance()->createUniDoubleGenerator(0.0, 1.0);
    double clearance =  clearanceGen.next();
    for( int i = 0; i < countOfImmuneAgentsAtVirionLocation; ++i )
    {
        clearance = clearance * clearanceProbScaler;
    }

    // If the virus get's cleared, mark it as Dead.
    if( clearance > 1 - clearanceProbability )
    {
        virionState = Dead;
    }
}



/**********************
*   VirionAgent::attemptToInfectCell - Gets the epithelial cell at the current grid position, and attempts to infect it. 
*   Implementation of the AttemptToInfectCell submodel.
**********************/
void VirionAgent::attemptToInfectCell(EpithelialCellAgent* theEpithelialCell, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    if(theEpithelialCell != nullptr)
    {
        // If the cell seems to be helathy. (NOTE: It could be infected but is not expressing any proteins, so there is no sign of infection ).
        if( theEpithelialCell->getExternalState() == EpithelialCellAgent::ExternalState::SeeminglyHealthy )
        {
            repast::DoubleUniformGenerator penetrationChanceGen = repast::Random::instance()->createUniDoubleGenerator( 0.0 , 1.0 );
            double toInfectTheCell = penetrationChanceGen.next();

            // Check if the virus manages to penetrate the cell based on the probability parameter.
            if( toInfectTheCell > 1-penetrationProbability )
            {
                // Infect the epithelial cell.
                theEpithelialCell->infect();

                // After infecting the cell we'll turn the virion agent to Contained, since it is no longer a free virion particle. 
                // Then the Contained virion particle will be removed from the simulation by the Virs_Cell_Model class, as it can no longer move and infect other cells.
                // This is to also ensure that the data tracking the count of free virions is correct.
                virionState = Contained;
            }
        }
    }
    else
    {
        std::cout<<"An incorrect pointer to the epithelial cell agent was passed to VirionAgent::attemptToInfectCell()! The penetration attempt cannot be executed."<<std::endl;
    }

    // If the virus is not contained (It could not infect the cell), then move to a neighbouring grid cell.
    if( virionState == Free_Virion)
    {
        move(discreteGridSpace);
    }
}



/**********************
*   VirionAgent::move - Function for moving an agent in the grid.
**********************/
void VirionAgent::move(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    if(discreteGridSpace == nullptr)
    {
        std::cout<<"An incorrect pointer to the grid projection was passed to VirionAgent::move! The move cannot be executed."<<std::endl;
        return;
    }

    // Get the current agent location from the grid projection.
    std::vector<int> currentLocation;
    discreteGridSpace->getLocation(agentId, currentLocation);

    repast::IntUniformGenerator moveGenerator = repast::Random::instance()->createUniIntGenerator(-1, 1);
    int moveX = 0;
    int moveY = 0;

    // Generate the agent's move stochasticly. The move is ensured to be to one of the 8 directly neighbouring cell. (Does not stay at the same spot)
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