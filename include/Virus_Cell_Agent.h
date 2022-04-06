/* Virus_Cell_Agent.h */
#ifndef VIRUS_CELL_AGENT
#define VIRUS_CELL_AGENT

/**********************
*   Include files
**********************/
#include "repast_hpc/AgentId.h"
#include "repast_hpc/SharedContext.h"

// Spatial Projection includes
#include "repast_hpc/SharedDiscreteSpace.h"
#include "repast_hpc/GridComponents.h"

/**********************
* Base Agent Class 
**********************/
class VirusCellInteractionAgents
{

protected:
    repast::AgentId     agentId;
    double              agentLifespan;
    int                 agentAge;

public:
    // Constructors
    VirusCellInteractionAgents(repast::AgentId theId, double theLifespan, int theAge);

    // Destructor
    virtual ~VirusCellInteractionAgents();

    /* Required Getters */
    virtual repast::AgentId& getId(){                   return agentId;    }
    virtual const repast::AgentId& getId() const {      return agentId;    }

    /* Getters for varibales which are shared by all agent types: the age and lifespan. */
    double getLifespan(){                                      return agentLifespan;      }
    double getAge(){                                  return agentAge;  }

    virtual void doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);

protected:
    virtual void move(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace) {}
};
#endif // VIRUS_CELL_AGENT