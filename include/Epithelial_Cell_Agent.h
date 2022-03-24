/* Epithelial_Cell_agent.h */
#ifndef EPITHELIAL_CELL_AGENT
#define EPITHELIAL_CELL_AGENT

/**********************
*   Include files
**********************/
#include "repast_hpc/AgentId.h"
#include "repast_hpc/SharedContext.h"

// Spatial Projection includes
#include "repast_hpc/SharedDiscreteSpace.h"

#include "Virus_Cell_Agent.h"

class EpithelialCellAgent: public VirusCellInteractionAgents
{
public:
    enum ExternalState{ displayingViralProtein, seeminglyHealthy, deadCell };
    enum InternalState{ dead, healthy, infected };
    enum NeighbouringCellModificationType{ noModification, toDivideInto, toInfect };


public:
    // Constructors
    EpithelialCellAgent(repast::AgentId theId);
    EpithelialCellAgent(repast::AgentId theId, double theLifespan, int theAge, int theInfectedLifespan, int theDivisionRate, int theTimeSinceLastDivision, 
                        double theReleaseDelay, double theDisplayVirProtDelay, double theExtracellularReleaseProb, double theCellToCellTransmissionProb);
    EpithelialCellAgent(repast::AgentId theId, double theLifespan, int theAge, InternalState theInternalState, 
                        ExternalState theExternalState, int theInfectedLifespan, int theInfectedTime, int theDivisionRate, int theTimeSinceLastDivision, double theReleaseDelay,  double theDisplayVirProtDelay,
                        NeighbouringCellModificationType theModificationToNeighbCell, repast::AgentId theNeighbouringCellToModify,  double theExtracellularReleaseProb, double theCellToCellTransmissionProb);

    // Destructor
    ~EpithelialCellAgent();

    /* Required Getters */
    virtual repast::AgentId& getId(){                   return agentId;    }
    virtual const repast::AgentId& getId() const {      return agentId;    }

    /* Getters specific to this these specific agents */
    int getInternalState(){                             return internalState;           }
    int getExternalState(){                             return externalState;           }
    int getInfectedLifespan(){                          return infectedLifespan;        }
    int getTimeInfected(){                              return timeInfected;            }
    int getTimeSinceLastDivision(){                     return timeSinceLastDivision;   }
    int getDivisionRate(){                              return divisionRate;            }

    int getTypeOfModifToNeighbCell(){                   return modificationToNeighbCell;}
    repast::AgentId getNeighbouringCellToModify(){      return neighbouringCellToModify;}

    double getReleaseDelay(){                           return releaseDelay;            }
    bool isCellToReleaseVirion(){                       return toReleaseVirion;         }

    double getDisplayVirProteinsDelay(){                return displayVirProteinsDelay; }

    double getExtracellularReleaseProb(){               return extracellularReleaseProb;}
    double getCellToCellTransmissionProb(){             return cellToCellTransmissionProb;}

    /* Setter */
    void set(int currentRank, double newLifespan, double newAge, InternalState newInternalState, ExternalState newExtState, int newInfectedLifespan ,int newInfectedTime, int newDivisionRate, int newTimeSinceLastDivision, double newReleaseDelay, 
            double newDisplayVirProteinsDelay, NeighbouringCellModificationType newModificationToNeighbCell, repast::AgentId newNeighbouringCellToModify,
            double newExtracellularReleaseProb, double newCellToCellTransmissionProb);

    void doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    void actHealthy(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    void actInfected(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    void infect(){                  internalState=infected;     }
    void eliminate(){               internalState=dead;     externalState = deadCell;}

private:
    InternalState internalState;
    ExternalState externalState;

    int infectedLifespan;
    int timeInfected;
  
    int timeSinceLastDivision;
    int divisionRate;

    NeighbouringCellModificationType modificationToNeighbCell;
    repast::AgentId neighbouringCellToModify;
    repast::AgentId idForNoNeighbourModification;

    double releaseDelay;
    bool toReleaseVirion;

    double displayVirProteinsDelay;

    // Probability of releasing a new virus particle in the extracellular space, when the cell starts producing the progeny virus
    double extracellularReleaseProb;

    // Probability of directly infecting a neighbouring cell as a form of new virus particle release
    double cellToCellTransmissionProb;
        
    
};

#endif