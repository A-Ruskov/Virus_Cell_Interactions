/* Virus_Cell_Agent.cpp */
// Implements the parent class for all agents found in the model.
// Implements the Agent Serializable Package Struct used for moving agents between processes and updating agent copies in the bufferzones.

#include "Virus_Cell_Agent.h"

/****************************
* VirusCellInteractionAgents CLASS Implementation
****************************/

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
*   VirusCellInteractionAgents::doStep - Function for an agent to do a step. Will be triggered on every step,
*   The agent will then do an action depending on the surrounding agents and its internal state. 
*   Left empty as each specific agent type will have its own implementation.
**********************/
void VirusCellInteractionAgents::doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace)
{
}