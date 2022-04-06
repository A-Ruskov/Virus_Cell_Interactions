/* DataCollection.cpp */
// Implements all DataSource classes.

/***********************************
********** DATA COLLECTION *********
***********************************/
#include "Data_Collection.h"

/******************************************
* Data Source class for tracking the count of alive epithelial cells in the model
******************************************/


/**********************
*   DataSource_EpithelialCellsCount::DataSource_EpithelialCellsCount - Constructor
**********************/
DataSource_EpithelialCellsCount::DataSource_EpithelialCellsCount(repast::SharedContext<VirusCellInteractionAgents>* theContext):
context(theContext)
{

}



/**********************
*   DataSource_EpithelialCellsCount::getData - Gets the agents of type Epithelial cell on this process and counts how many are alive
**********************/
int DataSource_EpithelialCellsCount::getData()
{
    int aliveEpithCellsCount = 0;

    std::vector<VirusCellInteractionAgents*> theEpithelialCells;
    
    // Get all local epithelial cell agents
    context->selectAgents(repast::SharedContext<VirusCellInteractionAgents>::LOCAL, theEpithelialCells, 0, false);
    std::vector<VirusCellInteractionAgents*>::iterator iter;

    // Count only the alive epithelial cells
    for( iter = theEpithelialCells.begin(); iter != theEpithelialCells.end(); ++iter )
    {

        EpithelialCellAgent* currentEpithCell = static_cast<EpithelialCellAgent*>(*iter);
        if(currentEpithCell->getInternalState() != EpithelialCellAgent::InternalState::Dead)
        {
            ++aliveEpithCellsCount;
        }

    } 
    
    return aliveEpithCellsCount;
}




/******************************************
* Data Source class for tracking the count of virions in the model
******************************************/

/**********************
*   DataSource_VirionsCount::DataSource_VirionsCount - Constructor
**********************/
DataSource_VirionsCount::DataSource_VirionsCount(repast::SharedContext<VirusCellInteractionAgents>* theContext):
context(theContext)
{

}


/**********************
*   DataSource_VirionsCount::getData - Gets the count of Virion Agents on this process
**********************/
int DataSource_VirionsCount::getData()
{
    std::vector<VirusCellInteractionAgents*> theVirions;
    
    context->selectAgents(repast::SharedContext<VirusCellInteractionAgents>::LOCAL, theVirions, 1, false, -1);
    int virionCount = theVirions.size();
    
    return virionCount;
}



/******************************************
* Data Source class for tracking the count of innate immune cells in the model
******************************************/

/**********************
*   DataSource_InnateImmuneCellsCount::DataSource_InnateImmuneCellsCount - Constructor
**********************/
DataSource_InnateImmuneCellsCount::DataSource_InnateImmuneCellsCount(repast::SharedContext<VirusCellInteractionAgents>* theContext):
context(theContext)
{

}


/**********************
*   DataSource_InnateImmuneCellsCount::getData - Gets the count of Innate Immune Cell Agents on this process
**********************/
int DataSource_InnateImmuneCellsCount::getData()
{
    std::vector<VirusCellInteractionAgents*> theInnateImmuneCells;
    
    context->selectAgents(repast::SharedContext<VirusCellInteractionAgents>::LOCAL, theInnateImmuneCells, 2, false, -1);
    int innateImmuneCellsCount = theInnateImmuneCells.size();
    
    return innateImmuneCellsCount;
}



/******************************************
* Data Source class for tracking the count of specialised immune cells in the model
******************************************/

/**********************
*   DataSource_SpecialisedImmuneCellsCount::DataSource_SpecialisedImmuneCellsCount - Constructor
**********************/
DataSource_SpecialisedImmuneCellsCount::DataSource_SpecialisedImmuneCellsCount(repast::SharedContext<VirusCellInteractionAgents>* theContext):
context(theContext)
{

}


/**********************
*   DataSource_SpecialisedImmuneCellsCount::getData - Gets the count of Specialised Immune Cell Agents on this process
**********************/
int DataSource_SpecialisedImmuneCellsCount::getData()
{
    std::vector<VirusCellInteractionAgents*> theSpecialisedImmuneCells;
    
    context->selectAgents(repast::SharedContext<VirusCellInteractionAgents>::LOCAL, theSpecialisedImmuneCells, 3, false, -1);
    int specialisedImmuneCellsCount = theSpecialisedImmuneCells.size();
    
    return specialisedImmuneCellsCount;
}



/******************************************
* Data Source class for tracking the count of infected epithelial cells in the model
******************************************/

/**********************
*   DataSource_InfectedEpithelialCellsCount::DataSource_InfectedEpithelialCellsCount - Constructor
**********************/
DataSource_InfectedEpithelialCellsCount::DataSource_InfectedEpithelialCellsCount(repast::SharedContext<VirusCellInteractionAgents>* theContext):
context(theContext)
{

}


/**********************
*   DataSource_InfectedEpithelialCellsCount::getData - Gets the count of Infected Epithelial Cell Agents on this process
**********************/
int DataSource_InfectedEpithelialCellsCount::getData()
{
    int infectedEpithCells = 0;

    std::vector<VirusCellInteractionAgents*> theEpithelialCells;
    
    // Get all local epithelial cell agents
    context->selectAgents(repast::SharedContext<VirusCellInteractionAgents>::LOCAL, theEpithelialCells, 0, false);
    std::vector<VirusCellInteractionAgents*>::iterator iter;

    // Count only the infected epithelial cells
    for( iter = theEpithelialCells.begin(); iter != theEpithelialCells.end(); ++iter )
    {

        EpithelialCellAgent* currentEpithCell = static_cast<EpithelialCellAgent*>(*iter);
        if(currentEpithCell->getInternalState() == EpithelialCellAgent::InternalState::Infected)
        {
            ++infectedEpithCells;
        }

    } 
    return infectedEpithCells;
}




/******************************************
* Data Source class for tracking the count of dead epithelial cells in the model
******************************************/

/**********************
*   DataSource_DeadEpithelialCellsCount::DataSource_DeadEpithelialCellsCount - Constructor
**********************/
DataSource_DeadEpithelialCellsCount::DataSource_DeadEpithelialCellsCount(repast::SharedContext<VirusCellInteractionAgents>* theContext):
context(theContext)
{

}


/**********************
*   DataSource_DeadEpithelialCellsCount::getData - Gets the count of Dead Epithelial Cell Agents on this process
**********************/
int DataSource_DeadEpithelialCellsCount::getData()
{
    int deadEpithCells = 0;

    std::vector<VirusCellInteractionAgents*> theEpithelialCells;
    
    // Get all local epithelial cell agents
    context->selectAgents(repast::SharedContext<VirusCellInteractionAgents>::LOCAL, theEpithelialCells, 0, false);
    std::vector<VirusCellInteractionAgents*>::iterator iter;

    // Count only the dead epithelial cells
    for( iter = theEpithelialCells.begin(); iter != theEpithelialCells.end(); ++iter )
    {

        EpithelialCellAgent* currentEpithCell = static_cast<EpithelialCellAgent*>(*iter);
        if(currentEpithCell->getInternalState() == EpithelialCellAgent::InternalState::Dead)
        {
            ++deadEpithCells;
        }

    } 
    return deadEpithCells;
}




/******************************************
* Data Source class for tracking the total count of agents in the model
******************************************/

/**********************
*   DataSource_TotalAgentsCount::DataSource_TotalAgentsCount - Constructor
**********************/
DataSource_TotalAgentsCount::DataSource_TotalAgentsCount(repast::SharedContext<VirusCellInteractionAgents>* theContext):
context(theContext)
{

}


/**********************
*   DataSource_TotalAgentsCount::getData - Gets the total count of agents on this process
**********************/
int DataSource_TotalAgentsCount::getData()
{
    std::vector<VirusCellInteractionAgents*> theAgents;
    
    // Get all local agents
    context->selectAgents(repast::SharedContext<VirusCellInteractionAgents>::LOCAL, theAgents, false);
    
    return theAgents.size();
}