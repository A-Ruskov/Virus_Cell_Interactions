/* Virus_Cell_Agent.cpp */

#include "Virus_Cell_Agent.h"



/**********************
*   VirusCellInteractionAgents::VirusCellInteractionAgents - Basic Constructor for the VirusCellInteractionAgents class.
**********************/
VirusCellInteractionAgents::VirusCellInteractionAgents(repast::AgentId theId): 
agentId(theId), 
agentLifespan(50), 
agentAge(0)
{ }



/**********************
*   VirusCellInteractionAgents::VirusCellInteractionAgents - Constructor for the VirusCellInteractionAgents class.
**********************/
VirusCellInteractionAgents::VirusCellInteractionAgents(repast::AgentId theId, double theLifespan, int theAge): 
agentId(theId), 
agentLifespan(theLifespan), 
agentAge(theAge)
{ }



/**********************
*   VirusCellInteractionAgents::~VirusCellInteractionAgents - Destructor for the VirusCellInteractionAgents class.
**********************/
VirusCellInteractionAgents::~VirusCellInteractionAgents()
{ 

}



/**********************
*   VirusCellInteractionAgents::set - Setter for the agent. Sets the current rank and updates the age and lifespan
**********************/
void VirusCellInteractionAgents::set(int currentRank, double newLifespan, double newAge)
{
    agentId.currentRank(currentRank);
    agentLifespan = newLifespan;
    agentAge = newAge;
}



/**********************
*   VirusCellInteractionAgents::doStep - Function for an agent to do a step. Will be triggered on every step,
*   The agent will then do an action depending on the surrounding agents and its internal state.
**********************/
void VirusCellInteractionAgents::doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
    agentAge++;
    std::cout << "I " << agentId << " aged " << agentAge << " am doing a step" << std::endl;
}




/**********************
* Serializable Agent Package Data
**********************/
VirusCellInteractionAgentPackage::VirusCellInteractionAgentPackage(){ }

VirusCellInteractionAgentPackage::VirusCellInteractionAgentPackage(int _id, int _rank, int _type, int _currentRank, int _lifespan, int _age, int _privateState, int _publicState, int _infectedLifespan, 
                int _timeInfected, int _divisionRate, int _timeSinceLastDivision, double _releaseDelay,  double _displayVirProteinsDelay, 
                double _extracellularReleaseProb, double _cellToCellTransmissionProb,
                double _penetrationProbability, double _clearanceProbability, double _clearanceProbScaler,
                double _infectedCellRecognitionProb, double _infectedCellEliminationProb,
                double _specialisedImmuneCellRecruitProb,
                int _typeOfChangeToMakeToAgent, int _agentToEditID, int _agentToEditStartRank, int _agentToEditType, int _agentToEditCurrRank):
id(_id), 
rank(_rank), 
type(_type),
currentRank(_currentRank), 
lifespan(_lifespan), 
age(_age),
privateState(_privateState),

// Specific to epithelial cell agents
publicState(_publicState),
infectedLifespan(_infectedLifespan),
timeInfected(_timeInfected),
divisionRate(_divisionRate),
timeSinceLastDivision(_timeSinceLastDivision),
releaseDelay(_releaseDelay),
displayVirProteinsDelay(_displayVirProteinsDelay),
extracellularReleaseProb(_extracellularReleaseProb),
cellToCellTransmissionProb(_cellToCellTransmissionProb),

// Specific to virion agents
penetrationProbability(_penetrationProbability),
clearanceProbability(_clearanceProbability),
clearanceProbScaler(_clearanceProbScaler),

// Specific to the two immune cell types
infectedCellRecognitionProb(_infectedCellRecognitionProb),
infectedCellEliminationProb(_infectedCellEliminationProb),

// Specific to the innate immune cell agents
specialisedImmuneCellRecruitProb(_specialisedImmuneCellRecruitProb),

// Specific for handling buffer zone agent changes propagating back to original agents
typeOfChangeToMakeToAgent(_typeOfChangeToMakeToAgent),
agentToEditID(_agentToEditID),
agentToEditStartRank(_agentToEditStartRank),
agentToEditType(_agentToEditType),
agentToEditCurrRank(_agentToEditCurrRank)
{ 

}