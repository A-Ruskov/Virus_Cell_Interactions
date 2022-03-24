/*  Specialised_Immune_Cell.cpp  */
#include "Specialised_Immune_Cell.h"
#include "Virus_Cell_Agent.h"
#include "Epithelial_Cell_Agent.h"

#include "repast_hpc/initialize_random.h"

#include "repast_hpc/Point.h"
#include "repast_hpc/Moore2DGridQuery.h"


/**********************
*   SpecialisedImmuneCellAgent::SpecialisedImmuneCellAgent - Basic Constructor for the SpecialisedImmuneCellAgent class.
**********************/
SpecialisedImmuneCellAgent::SpecialisedImmuneCellAgent(repast::AgentId theId) : 
VirusCellInteractionAgents(theId),
specialisedImmuneCellState(healthy),
toRecruitNewSpecImmuneCell(false)
{ 
    agentLifespan = 25; 
    agentAge = 0;
    infectedCellRecognitionProb = 0.7;
    infectedCellEliminationProb = 0.8;
}



/**********************
*   SpecialisedImmuneCellAgent::SpecialisedImmuneCellAgent - Constructor for the SpecialisedImmuneCellAgent class.
**********************/
SpecialisedImmuneCellAgent::SpecialisedImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, double theInfectedCellRecognitionProb, double theInfectedCellEliminationProb):
VirusCellInteractionAgents(theId, theLifespan, theAge),
specialisedImmuneCellState(healthy),
infectedCellRecognitionProb(theInfectedCellRecognitionProb),
infectedCellEliminationProb(theInfectedCellEliminationProb),
toRecruitNewSpecImmuneCell(false)
{
}


/**********************
*   SpecialisedImmuneCellAgent::SpecialisedImmuneCellAgent - Constructor for the SpecialisedImmuneCellAgent class.
**********************/
SpecialisedImmuneCellAgent::SpecialisedImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, SpecialisedImmuneCellStates theState, double theInfectedCellRecognitionProb, double theInfectedCellEliminationProb):
VirusCellInteractionAgents(theId, theLifespan, theAge),
specialisedImmuneCellState(theState),
infectedCellRecognitionProb(theInfectedCellRecognitionProb),
infectedCellEliminationProb(theInfectedCellEliminationProb),
toRecruitNewSpecImmuneCell(false)
{
}




/**********************
*   SpecialisedImmuneCellAgent::~SpecialisedImmuneCellAgent - Destructor for the SpecialisedImmuneCellAgent class.
**********************/
SpecialisedImmuneCellAgent::~SpecialisedImmuneCellAgent()
{ 

}



/**********************
*   SpecialisedImmuneCellAgent::set - Setter for the agent. Sets the current rank and updates the age, lifespan and cell state
**********************/
void SpecialisedImmuneCellAgent::set(int currentRank, double newLifespan, double newAge, SpecialisedImmuneCellStates newCellState, double newInfectedCellRecognitionProb, double newInfectedCellEliminationProb)
{
    agentId.currentRank(currentRank);
    agentLifespan = newLifespan;
    agentAge = newAge;
    specialisedImmuneCellState = newCellState;
    infectedCellRecognitionProb = newInfectedCellRecognitionProb;
    infectedCellEliminationProb = newInfectedCellEliminationProb;
}



/**********************
*   SpecialisedImmuneCellAgent::doStep - Function for an agent to do a step. Will be triggered on every step,
*   The agent will then do an action depending on the surrounding agents and its internal state.
**********************/
void SpecialisedImmuneCellAgent::doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    ++agentAge;
    if(agentAge > agentLifespan && specialisedImmuneCellState != dead)
    {
        specialisedImmuneCellState = dead;
    }

    if( specialisedImmuneCellState != dead )
    {
        attemptToDetectInfectedCell(discreteGridSpace);

        // Move the specialised immune cell, then it could attepmt to detect another infected cell on the next time step.
        move(discreteGridSpace);
    }
}




/**********************
*   SpecialisedImmuneCellAgent::move - Function for moving an agent in the grid.
**********************/
void SpecialisedImmuneCellAgent::move(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
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
*   SpecialisedImmuneCellAgent::attemptToInfectCell - Gets the epithelial cell at the current grid position, attempts to find if it is infected and eliminate it.
**********************/
void SpecialisedImmuneCellAgent::attemptToDetectInfectedCell(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    toRecruitNewSpecImmuneCell = false;

    std::vector<int> immuneCellLoc;
    discreteGridSpace->getLocation(agentId, immuneCellLoc);
    repast::Point<int> locationPoint(immuneCellLoc);
    
    // Get the at this grid cell
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
            if( theEpithelialCell->getExternalState() == theEpithelialCell->displayingViralProtein )
            {
                repast::DoubleUniformGenerator detectionProbGen = repast::Random::instance()->createUniDoubleGenerator( 0.0 , 1.0 );
                double toDetectInfectedCell = detectionProbGen.next();

                if( toDetectInfectedCell > 1-infectedCellRecognitionProb )
                {
                    toRecruitNewSpecImmuneCell = true;

                    repast::DoubleUniformGenerator infectedCellElimProbGen = repast::Random::instance()->createUniDoubleGenerator( 0.0 , 1.0 );
                    double toEliminateCell = infectedCellElimProbGen.next();

                    if( toEliminateCell > 1-infectedCellEliminationProb )
                    {
                        theEpithelialCell->eliminate();
                    }

                }

                return;
            }
            
            break;
        }
    }
}