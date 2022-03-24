/*  Innate_Immune_Cell.cpp  */
#include "Innate_Immune_Cell.h"
#include "Virus_Cell_Agent.h"
#include "Epithelial_Cell_Agent.h"

#include "repast_hpc/initialize_random.h"

#include "repast_hpc/Point.h"
#include "repast_hpc/Moore2DGridQuery.h"


/**********************
*   InnateImmuneCellAgent::InnateImmuneCellAgent - Basic Constructor for the InnateImmuneCellAgent class.
**********************/
InnateImmuneCellAgent::InnateImmuneCellAgent(repast::AgentId theId) : 
VirusCellInteractionAgents(theId),
innateImmuneCellState(healthy),
toRecruitNewSpecImmuneCell(false),
toRecruitNewInnateImmuneCell(false),
specialisedImmuneCellRecruitProb(0.5)
{ 
    agentLifespan = 20; 
    agentAge = 0;
    infectedCellRecognitionProb = 0.5;
    infectedCellEliminationProb = 0.3;
}



/**********************
*   InnateImmuneCellAgent::InnateImmuneCellAgent - Constructor for the InnateImmuneCellAgent class.
**********************/
InnateImmuneCellAgent::InnateImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, double theInfectedCellRecognitionProb, double theInfectedCellEliminationProb, double theSpecialisedImmuneCellRecruitProb ):
VirusCellInteractionAgents(theId, theLifespan, theAge),
innateImmuneCellState(healthy),
infectedCellRecognitionProb(theInfectedCellRecognitionProb),
infectedCellEliminationProb(theInfectedCellEliminationProb),
toRecruitNewSpecImmuneCell(false),
toRecruitNewInnateImmuneCell(false),
specialisedImmuneCellRecruitProb(theSpecialisedImmuneCellRecruitProb)
{
}


/**********************
*   InnateImmuneCellAgent::InnateImmuneCellAgent - Constructor for the InnateImmuneCellAgent class.
**********************/
InnateImmuneCellAgent::InnateImmuneCellAgent(repast::AgentId theId, double theLifespan, int theAge, InnateImmuneCellStates theState, double theInfectedCellRecognitionProb, double theInfectedCellEliminationProb, double theSpecialisedImmuneCellRecruitProb):
VirusCellInteractionAgents(theId, theLifespan, theAge),
innateImmuneCellState(theState),
infectedCellRecognitionProb(theInfectedCellRecognitionProb),
infectedCellEliminationProb(theInfectedCellEliminationProb),
toRecruitNewSpecImmuneCell(false),
toRecruitNewInnateImmuneCell(false),
specialisedImmuneCellRecruitProb(theSpecialisedImmuneCellRecruitProb)
{
}




/**********************
*   InnateImmuneCellAgent::~InnateImmuneCellAgent - Destructor for the InnateImmuneCellAgent class.
**********************/
InnateImmuneCellAgent::~InnateImmuneCellAgent()
{ 

}



/**********************
*   InnateImmuneCellAgent::set - Setter for the agent. Sets the current rank and updates the age, lifespan and cell state
**********************/
void InnateImmuneCellAgent::set(int currentRank, double newLifespan, double newAge, InnateImmuneCellStates newCellState, double newInfectedCellRecognitionProb, double newInfectedCellEliminationProb, double newSpecialisedImmuneCellRecruitProb)
{
    agentId.currentRank(currentRank);
    agentLifespan = newLifespan;
    agentAge = newAge;
    innateImmuneCellState = newCellState;
    infectedCellRecognitionProb = newInfectedCellRecognitionProb;
    infectedCellEliminationProb = newInfectedCellEliminationProb;
    specialisedImmuneCellRecruitProb = newSpecialisedImmuneCellRecruitProb;
}



/**********************
*   InnateImmuneCellAgent::doStep - Function for an agent to do a step. Will be triggered on every step,
*   The agent will then do an action depending on the surrounding agents and its internal state.
**********************/
void InnateImmuneCellAgent::doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    ++agentAge;
    if(agentAge > agentLifespan && innateImmuneCellState != dead)
    {
        innateImmuneCellState = dead;
    }

    if( innateImmuneCellState != dead )
    {
        attemptToDetectInfectedCell(discreteGridSpace);

        // Move the innate immune cell, then it could attepmt to detect another infected cell on the next time step.
        move(discreteGridSpace);
    }
}





/**********************
*   InnateImmuneCellAgent::move - Function for moving an agent in the grid.
**********************/
void InnateImmuneCellAgent::move(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
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
*   InnateImmuneCellAgent::attemptToInfectCell - Gets the epithelial cell at the current grid position, attempts to find if it is infected and eliminate it.
**********************/
void InnateImmuneCellAgent::attemptToDetectInfectedCell(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    toRecruitNewSpecImmuneCell = false;
    toRecruitNewInnateImmuneCell = false;

    std::vector<int> immuneCellLoc;
    discreteGridSpace->getLocation(agentId, immuneCellLoc);
    repast::Point<int> locationPoint(immuneCellLoc);
    
    // Get the at this grid cell
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
        else if( (*agentsAtThisCellIter)->getId().agentType() == 3 )
        {
            specialisedImmCellsCountHere++;
        }

    }


    if( theEpithelialCellToCheck != nullptr )
    {
        // If the cell seems to be helathy. (NOTE: It could be infected but is not expressing any proteins, so there is no sign of infection ).
        if( theEpithelialCellToCheck->getExternalState() == EpithelialCellAgent::ExternalState::displayingViralProtein )
        {
            repast::DoubleUniformGenerator detectionProbGen = repast::Random::instance()->createUniDoubleGenerator( 0.0 , 1.0 );
            double toDetectInfectedCell = detectionProbGen.next();

            if( toDetectInfectedCell > 1-infectedCellRecognitionProb )
            {
                // If there are no specialised immune cells here, then recruit more innate cells. Innate cells act as first line of defense
                // They then go down in number when the specialised immune cells come in to kill the infection.
                if(specialisedImmCellsCountHere == 0 )
                {
                    toRecruitNewInnateImmuneCell = true;
                }


                repast::DoubleUniformGenerator newSpecImmuneCellRecruitProbGen = repast::Random::instance()->createUniDoubleGenerator(0.0, 1.0);
                double toRecruitSpecialisedCell = newSpecImmuneCellRecruitProbGen.next();

                if( toRecruitSpecialisedCell > 1 - specialisedImmuneCellRecruitProb )
                {
                    toRecruitNewSpecImmuneCell = true;
                }

                repast::DoubleUniformGenerator infectedCellElimProbGen = repast::Random::instance()->createUniDoubleGenerator( 0.0 , 1.0 );
                double toEliminateCell = infectedCellElimProbGen.next();

                if( toEliminateCell > 1-infectedCellEliminationProb )
                {
                    theEpithelialCellToCheck->eliminate();
                }

            }
        }
    }
}


