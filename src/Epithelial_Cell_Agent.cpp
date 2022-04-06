/* Epithelial_Cell_Agent.cpp */
// Implements the Epithelial Cell Agents

/**********************
*   INCLUDE FILES
**********************/
#include "Epithelial_Cell_Agent.h"
#include "Virus_Cell_Agent.h"
#include "repast_hpc/Moore2DGridQuery.h"
#include "repast_hpc/Point.h"


/**********************
*   EpithelialCellAgent::EpithelialCellAgent - Constructor for the EpithelialCellAgent class.
**********************/
EpithelialCellAgent::EpithelialCellAgent(repast::AgentId theId, double theLifespan, int theAge, double theInfectedLifespan, double theDivisionRate, int theTimeSinceLastDivision,
                                         double theReleaseDelay, double theDisplayVirProtDelay, double theExtracellularReleaseProb, double theCellToCellTransmissionProb,
                                         double theVirionReleaseRate):
VirusCellInteractionAgents(theId, theLifespan, theAge) ,
internalState(Healthy),
externalState(SeeminglyHealthy),
timeInfected(0),
infectedLifespan(theInfectedLifespan),
divisionRate(theDivisionRate),
timeSinceLastDivision(theTimeSinceLastDivision),
modificationToNeighbCell(NoModification),
neighbouringCellToModify(-1, -1, -1, -1),
idForNoNeighbourModification(-1, -1, -1, -1),
releaseDelay(theReleaseDelay),
displayVirProteinsDelay(theDisplayVirProtDelay),
extracellularReleaseProb(theExtracellularReleaseProb),
cellToCellTransmissionProb(theCellToCellTransmissionProb),
virionReleaseRate(theVirionReleaseRate),
countOfVirionsToRelease(0),
virionReleaseRemainder(0.0)
{
}



/**********************
*   EpithelialCellAgent::EpithelialCellAgent - Constructor for the EpithelialCellAgent class.
**********************/
EpithelialCellAgent::EpithelialCellAgent(repast::AgentId theId, double theLifespan, int theAge, InternalState theInternalState, 
                        ExternalState theExternalState, double theInfectedLifespan, int theInfectedTime, double theDivisionRate, int theTimeSinceLastDivision, double theReleaseDelay, double theDisplayVirProtDelay,
                        NeighbouringCellModificationType theModificationToNeighbCell, repast::AgentId theNeighbouringCellToModify,
                        double theExtracellularReleaseProb, double theCellToCellTransmissionProb,
                        double theVirionReleaseRate, int theCountOfVirionsToRelease, double theVirionReleaseRemainder):
VirusCellInteractionAgents(theId, theLifespan, theAge),
internalState(theInternalState),
externalState(theExternalState),
infectedLifespan(theInfectedLifespan),
timeInfected(theInfectedTime),
divisionRate(theDivisionRate),
timeSinceLastDivision(theTimeSinceLastDivision),
modificationToNeighbCell(theModificationToNeighbCell),
neighbouringCellToModify(theNeighbouringCellToModify),
idForNoNeighbourModification(-1, -1, -1, -1),
releaseDelay(theReleaseDelay),
displayVirProteinsDelay(theDisplayVirProtDelay),
extracellularReleaseProb(theExtracellularReleaseProb),
cellToCellTransmissionProb(theCellToCellTransmissionProb),
virionReleaseRate(theVirionReleaseRate),
countOfVirionsToRelease(theCountOfVirionsToRelease),
virionReleaseRemainder(theVirionReleaseRemainder)
{

}



/**********************
*   EpithelialCellAgent::~EpithelialCellAgent - Destructor for the EpithelialCellAgent class.
**********************/
EpithelialCellAgent::~EpithelialCellAgent()
{ 

}



/**********************
*   EpithelialCellAgent::set - Setter for the agent. Sets all related variable values. This is needed in order to update the agents in the buffer zone.
**********************/
void EpithelialCellAgent::set(int currentRank, double newLifespan, double newAge, InternalState newInternalState, ExternalState newExtState, 
                              double newInfectedLifespan ,int newInfectedTime, double newDivisionRate, int newTimeSinceLastDivision, double newReleaseDelay, 
                              double newDisplayVirProteinsDelay, NeighbouringCellModificationType newModificationToNeighbCell, repast::AgentId newNeighbouringCellToModify, 
                              double newExtracellularReleaseProb, double newCellToCellTransmissionProb,
                              double newVirionReleaseRate, int newCountOfVirionsToRelease, double newVirionReleaseRemainder)
{
    agentId.currentRank(currentRank);
    agentLifespan = newLifespan;
    agentAge = newAge;
    internalState = newInternalState;
    externalState = newExtState;
    infectedLifespan = newInfectedLifespan;
    timeInfected = newInfectedTime;
    divisionRate = newDivisionRate;
    timeSinceLastDivision = newTimeSinceLastDivision;
    releaseDelay = newReleaseDelay;
    displayVirProteinsDelay = newDisplayVirProteinsDelay;
    modificationToNeighbCell = newModificationToNeighbCell;
    neighbouringCellToModify = newNeighbouringCellToModify;
    extracellularReleaseProb = newExtracellularReleaseProb;
    cellToCellTransmissionProb = newCellToCellTransmissionProb;
    virionReleaseRate = newVirionReleaseRate;
    countOfVirionsToRelease = newCountOfVirionsToRelease;
    virionReleaseRemainder = newVirionReleaseRemainder;
}



/**********************
*   EpithelialCellAgent::doStep - Function for an agent to do a step. Will be triggered on every step,
*   The agent will then do an action depending on the surrounding agents and its internal state.
**********************/
void EpithelialCellAgent::doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    // Reset the identifiers for releasing virions and modifying neighbouring cells to the default values. These will be set to specific values if needed.
    countOfVirionsToRelease = 0;
    modificationToNeighbCell = NoModification;
    neighbouringCellToModify = idForNoNeighbourModification;

    // Increment the epithelial cell agent age.
    ++agentAge;
    // If the age exceeds the cell's lifespan, the cell dies. Change both internal and external cell states to dead.
    if(agentAge > agentLifespan && internalState != Dead)
    {
        internalState = Dead;
        externalState = DeadCell;
    }

    // If the cell is not dead then it can act
    if( internalState != Dead )
    {
        // If the cell is Infected
        if( internalState == Infected )
        {
            actInfected(discreteGridSpace);
        }
        // If the cell is healthy, then act healthy.
        else
        {
            actHealthy(discreteGridSpace);
        } 
    }
}



/**********************
*   EpithelialCellAgent::actHealthy - If the epithelial cell agent is heallthy, then act normally.
*   Track its division rate and if it is ready to divide, choose a neighbouring cell to divide into.
**********************/
void EpithelialCellAgent::actHealthy(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    ++timeSinceLastDivision;
    // Check if the cell is ready to divide. If it is find which of its 8 neighbouring epithelial cells is dead.
    // If there is at least one, then the cell will divide into that position.
    // The division of a cell is modelled as reviving a dead cell, however the new cell at that place will use the same EpithelialCellAgent object.
    // Will have different parameters to the previously dead cell, in order to actually represent a new epithelial cell.
    // The reset of the cell will be executed by the Virus_Cell_Model class. Here we just need to identify the cell which is to be revived.
    if(timeSinceLastDivision > divisionRate)
    {
        std::vector<int> myLocation;
        discreteGridSpace->getLocation(agentId, myLocation);
        repast::Point<int> locationPoint(myLocation);
        
        // Get the agents in the surrounding area
        repast::Moore2DGridQuery<VirusCellInteractionAgents> theQuery(discreteGridSpace);
        std::vector<VirusCellInteractionAgents*> surroundingAgents;
        theQuery.query(locationPoint, 1, false, surroundingAgents);

        std::vector<VirusCellInteractionAgents*>::iterator surroundingAgentsIter;

        // We only need to consider epithelial cell agents which are dead, as we consider those as empty space, where the cell can divide into
        // Find those cells and store them in the vector.
        std::vector<VirusCellInteractionAgents*> neighbouringEpithelialCellsToConsider;
        for( surroundingAgentsIter = surroundingAgents.begin(); surroundingAgentsIter != surroundingAgents.end(); ++surroundingAgentsIter )
        {
            if( (*surroundingAgentsIter)->getId().agentType() == 0 )
            {
                EpithelialCellAgent* neighbouringEpithelialCell = static_cast<EpithelialCellAgent*>(*surroundingAgentsIter);
                if( neighbouringEpithelialCell->getExternalState() == DeadCell )
                {
                    neighbouringEpithelialCellsToConsider.push_back(*surroundingAgentsIter);
                }  
            }
        }

        // Arbitrarily choose one of the possible cells to divide into
        if( neighbouringEpithelialCellsToConsider.size() > 0 )
        {
            repast::IntUniformGenerator neighbouringCellChoiceGen = repast::Random::instance()->createUniIntGenerator(0, neighbouringEpithelialCellsToConsider.size()-1);
            int cellToDivideIntoIndex = neighbouringCellChoiceGen.next();
            VirusCellInteractionAgents* theAgentToRevive = neighbouringEpithelialCellsToConsider[cellToDivideIntoIndex];

            // Store the id of the agent whic is to be "revived" (turned into a brand new cell)
            // Set the modification to neighbouring cell to be "ToDivideInto". That means that the cell wants to modify a neighbouring agent 
            //(that could be local or non-local agent). Then set what this modification should be. In this case we want to divide into that cell, 
            // which means revive the cell with new parameters. The model class will get the stored id and then use it in order to create the new cell into the existing dead cell object.
            modificationToNeighbCell = ToDivideInto;
            neighbouringCellToModify = theAgentToRevive->getId();
        }

        timeSinceLastDivision = 0;
    } 
}



/**********************
*   EpithelialCellAgent::actInfected - If the epithelial cell agent is infected, then act accordingly.
**********************/
void EpithelialCellAgent::actInfected(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    ++timeInfected;

    // If the cell has been infected for more than its infected lifespan, then the cell's states should be turned to dead.
    // Then return since a dead cell cannot perform any further actions.
    if( timeInfected > infectedLifespan )
    {
        internalState = Dead;
        externalState = DeadCell;
        return;
    }

    if( timeInfected > displayVirProteinsDelay &&  externalState != DisplayingViralProtein)
    {
        externalState = DisplayingViralProtein;
    }

    // If the cell is ready to release new virus particles, then mark itself as ready to release.
    if( timeInfected > releaseDelay )
    {
        externalState = DisplayingViralProtein;

        // Perform new virion release and/or cell to cell infection.
        releaseProgenyVirus();
        cellToCellInfection(discreteGridSpace);        
    }
}



/**********************
*   EpithelialCellAgent::releaseProgenyVirus - Calculates the count of new virus particles which need to be released in the extracellular space.
*   The Virus_Cell_Model class will then add the required number of agents to the simulation.
*   Implements the ReleaseProgenyVirus submodel.
**********************/
void EpithelialCellAgent::releaseProgenyVirus()
{
    repast::DoubleUniformGenerator releaseProbabilitiesGen = repast::Random::instance()->createUniDoubleGenerator(0.0, 1.0);

    // Attempt releasing a new virus particle in the extracellular space. 
    // The cell may/may not release new virus particles depending on the specifics of the Virus-Cell Interaction
    if( releaseProbabilitiesGen.next() > 1 - extracellularReleaseProb )
    {
        countOfVirionsToRelease = (int)(virionReleaseRate + virionReleaseRemainder);
        virionReleaseRemainder = (virionReleaseRate + virionReleaseRemainder) - (double)countOfVirionsToRelease;
    }
}



/**********************
*   EpithelialCellAgent::cellToCellInfection - Function for infecting a directly neighbouring cell which seems to be healthy.
*   Implements the CellToCellInfection submodel.
**********************/
void EpithelialCellAgent::cellToCellInfection(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    if(discreteGridSpace == nullptr)
    {
        std::cout<<"An incorrect pointer to the grid projection was passed to EpithelialCellAgent::cellToCellInfection! The cell to cell infection attempt cannot be executed."<<std::endl;
        return;
    }

    repast::DoubleUniformGenerator releaseProbabilitiesGen = repast::Random::instance()->createUniDoubleGenerator(0.0, 1.0);
    // Attempt a cell to cell infection of a neighbouring cell. Thaat will depend on a probability value.
    if( releaseProbabilitiesGen.next() > 1 - cellToCellTransmissionProb )
    {
        std::vector<int> myLocation;
        discreteGridSpace->getLocation(agentId, myLocation);
        repast::Point<int> locationPoint(myLocation);
        
        // Get the agents in the surrounding area
        repast::Moore2DGridQuery<VirusCellInteractionAgents> theQuery(discreteGridSpace);
        std::vector<VirusCellInteractionAgents*> surroundingAgents;
        theQuery.query(locationPoint, 1, false, surroundingAgents);

        std::vector<VirusCellInteractionAgents*>::iterator surroundingAgentsIter;

        // We only need to consider epithelial cell agents which appear to be healthy as we want to infect them
        // Find those cells and store them in the vector.
        std::vector<VirusCellInteractionAgents*> neighbouringEpithelialCellsToConsider;
        for( surroundingAgentsIter = surroundingAgents.begin(); surroundingAgentsIter != surroundingAgents.end(); ++surroundingAgentsIter )
        {
            if( (*surroundingAgentsIter)->getId().agentType() == 0 )
            {
                EpithelialCellAgent* neighbouringEpithelialCell = static_cast<EpithelialCellAgent*>(*surroundingAgentsIter);
                if( neighbouringEpithelialCell->getExternalState() == SeeminglyHealthy )
                {
                    neighbouringEpithelialCellsToConsider.push_back(*surroundingAgentsIter);
                }  
            }
        }

        // Arbitrarily choose one of the possible cells to infect
        if( neighbouringEpithelialCellsToConsider.size() > 0 )
        {
            repast::IntUniformGenerator neighbouringCellChoiceGen = repast::Random::instance()->createUniIntGenerator(0, neighbouringEpithelialCellsToConsider.size()-1);
            int cellToInfectIndex = neighbouringCellChoiceGen.next();
            VirusCellInteractionAgents* theAgentToInfect = neighbouringEpithelialCellsToConsider[cellToInfectIndex];

            // Set the modification to neighbouring cell to be "ToInfect". That means that the cell wants to modify a neighbouring agent 
            //(that could be local or non-local agent). Then set what this modification should be. In this case we want to infect that cell, 
            // The model class will get the stored id and then use it in order to infect the agent.
            modificationToNeighbCell = ToInfect;
            neighbouringCellToModify = theAgentToInfect->getId();;
        }
    }
}