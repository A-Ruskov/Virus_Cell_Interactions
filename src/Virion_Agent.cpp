/*  Virion_Agent.cpp  */
#include "Virion_Agent.h"
#include "Virus_Cell_Agent.h"
#include "Epithelial_Cell_Agent.h"

#include "repast_hpc/initialize_random.h"

#include "repast_hpc/Point.h"
#include "repast_hpc/Moore2DGridQuery.h"


/**********************
*   VirionAgent::VirionAgent - Basic Constructor for the VirionAgent class.
**********************/
VirionAgent::VirionAgent(repast::AgentId theId) : 
VirusCellInteractionAgents(theId),
virionState(free_virion),
penetrationProbability(0.4),
clearanceProbability(0.2),
clearanceProbScaler(1.2)
{ 
    agentLifespan = 15; 
    agentAge = 0;
}



/**********************
*   VirionAgent::VirionAgent - Constructor for the VirionAgent class.
**********************/
VirionAgent::VirionAgent(repast::AgentId theId, double theLifespan, int theAge, double thePenetrationProb, double theClearanceProbability, double theClearanceProbScaler):
VirusCellInteractionAgents(theId, theLifespan, theAge),
virionState(free_virion),
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
*   VirionAgent::set - Setter for the agent. Sets the current rank and updates the age and lifespan
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
    ++agentAge;
    if(agentAge > agentLifespan && virionState != dead)
    {
        virionState = dead;
    }
    else
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
            // Count how many immune cell agents are in the same grid cell
            else if( (*agentsAtThisCellIter)->getId().agentType() == 2 || (*agentsAtThisCellIter)->getId().agentType() == 3 )
            {
                immuneCellAgentsCount++;
            }
        }


        // For each immune cell at the same grid position, the clearance probability should grow. The more modelled immune cells there are here,
        // The more the chance of other unmodelled immune mechanisms to happen there as well.
        repast::DoubleUniformGenerator clearanceGen = repast::Random::instance()->createUniDoubleGenerator(0.0, 1.0);
        double clearance =  clearanceGen.next();
        for( int i = 0; i < immuneCellAgentsCount; ++i )
        {
            clearance = clearance * clearanceProbScaler;
        }

        // If the virus get's cleared, mark it as dead.
        if( clearance > 1 - clearanceProbability )
        {
            virionState = dead;
        }

        // If the virus is not "dead" (did not get cleared), then walk around and try to infect cells.
        if( virionState != dead )
        {
            attemptToInfectCell(theEpithelialCellHere);

            if(virionState != dead)
            {
                // Move the virus particle, if it could not infect the cell.
                move(discreteGridSpace);
            }
        }
    }


}



/**********************
*   VirionAgent::move - Function for moving an agent in the grid.
**********************/
void VirionAgent::move(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    std::vector<int> currentLocation;
    discreteGridSpace->getLocation(agentId, currentLocation);
    std::vector<int> newLocation;

    repast::IntUniformGenerator moveGenerator = repast::Random::instance()->createUniIntGenerator(-1, 1);
    newLocation.push_back(currentLocation[0] + moveGenerator.next());
    newLocation.push_back(currentLocation[1] + moveGenerator.next());

    repast::Point<int> movePoint(newLocation);

    discreteGridSpace->moveTo(agentId, movePoint);
}



/**********************
*   VirionAgent::attemptToInfectCell - Gets the epithelial cell at the current grid position, and attempts to infect it.
**********************/
void VirionAgent::attemptToInfectCell(EpithelialCellAgent* theEpithelialCell)
{
    if(theEpithelialCell != nullptr)
    {
        // If the cell seems to be helathy. (NOTE: It could be infected but is not expressing any proteins, so there is no sign of infection ).
        if( theEpithelialCell->getExternalState() == theEpithelialCell->seeminglyHealthy )
        {
            repast::DoubleUniformGenerator penetrationChanceGen = repast::Random::instance()->createUniDoubleGenerator( 0.0 , 1.0 );
            double toInfectTheCell = penetrationChanceGen.next();

            if( toInfectTheCell > 1-penetrationProbability )
            {
                theEpithelialCell->infect();

                // After infecting the cell we'll turn the virion agent to dead, since it is no longer a free virion particle. 
                // This is to ensure that the data tracking the count of free virions is correct.
                virionState = dead;
            }

        }
    }
}
