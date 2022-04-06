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

// Include the file containing the parent class
#include "Virus_Cell_Agent.h"

class EpithelialCellAgent: public VirusCellInteractionAgents
{
public:
    // The enum with the set of external states of the cell. The states which can be sensed by other agetns
    enum ExternalState{ DisplayingViralProtein, SeeminglyHealthy, DeadCell };

    // The enum with the set of internal states of the cell.
    enum InternalState{ Dead, Healthy, Infected };

    // The enum holding all potential types of modification that an epithelial cell can do to a neighbouring cell.
    enum NeighbouringCellModificationType{ NoModification, ToDivideInto, ToInfect };


public:
    /* Constructors */
    EpithelialCellAgent(repast::AgentId theId, double theLifespan, int theAge, double theInfectedLifespan, double theDivisionRate, int theTimeSinceLastDivision, 
                        double theReleaseDelay, double theDisplayVirProtDelay, double theExtracellularReleaseProb, double theCellToCellTransmissionProb, 
                        double theVirionReleaseRate);
    EpithelialCellAgent(repast::AgentId theId, double theLifespan, int theAge, InternalState theInternalState, 
                        ExternalState theExternalState, double theInfectedLifespan, int theInfectedTime, double theDivisionRate, int theTimeSinceLastDivision, double theReleaseDelay,  double theDisplayVirProtDelay,
                        NeighbouringCellModificationType theModificationToNeighbCell, repast::AgentId theNeighbouringCellToModify,  double theExtracellularReleaseProb, double theCellToCellTransmissionProb,
                        double theVirionReleaseRate, int theCountOfVirionsToRelease, double theVirionReleaseRemainder);

    // Destructor
    ~EpithelialCellAgent();

    /* Required Getters */
    virtual repast::AgentId& getId(){                   return agentId;    }
    virtual const repast::AgentId& getId() const {      return agentId;    }

    /* Getters specific to this these specific agents */
    int getInternalState(){                             return internalState;           }
    int getExternalState(){                             return externalState;           }

    double getInfectedLifespan(){                          return infectedLifespan;        }
    int getTimeInfected(){                              return timeInfected;            }

    int getTimeSinceLastDivision(){                     return timeSinceLastDivision;   }
    double getDivisionRate(){                              return divisionRate;            }

    int getTypeOfModifToNeighbCell(){                   return modificationToNeighbCell;}
    repast::AgentId getNeighbouringCellToModify(){      return neighbouringCellToModify;}

    double getReleaseDelay(){                           return releaseDelay;            }
    double getDisplayVirProteinsDelay(){                return displayVirProteinsDelay; }

    double getExtracellularReleaseProb(){               return extracellularReleaseProb;}
    double getCellToCellTransmissionProb(){             return cellToCellTransmissionProb;}

    double getVirionReleaseRate(){                      return virionReleaseRate; }
    int getVirionCountToRelease(){                      return countOfVirionsToRelease; }
    double getVirionReleaseRemainder(){                 return virionReleaseRemainder; }

    // Setter which sets all state variables and parameters of the agents. 
    // This setter is used only for updating agent copies at the buffer zone. It ensures that the non-local agents copies are always up-to-date with the original.
    void set(int currentRank, double newLifespan, double newAge, InternalState newInternalState, ExternalState newExtState,
            double newInfectedLifespan ,int newInfectedTime, double newDivisionRate, int newTimeSinceLastDivision, double newReleaseDelay, 
            double newDisplayVirProteinsDelay, NeighbouringCellModificationType newModificationToNeighbCell, repast::AgentId newNeighbouringCellToModify,
            double newExtracellularReleaseProb, double newCellToCellTransmissionProb, 
            double newVirionReleaseRate, int newCountOfVirionsToRelease, double newVirionReleaseRemainder);

    void doStep(repast::SharedContext<VirusCellInteractionAgents>* context, repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    
    // Function which the virion agents can use to infect an epithelial cell agent.
    void infect(){                  internalState = Infected;     }

    // Function which the two immune cell agent types can use to eliminate the epithelial cell agent when it is infected.
    void eliminate(){               internalState = Dead;     externalState = DeadCell;}

private:
    void actHealthy(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    void actInfected(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);
    void releaseProgenyVirus();
    void cellToCellInfection(repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace);

private:
    // The internal state of the cell. It dictates the way the agent acts.
    InternalState internalState;

    // The external state of the cell. It is what the other agent types can sense, in order to make their decisions on whether to do a particular action.
    // This is what the state of the cell appears to be from the outside.
    ExternalState externalState;

    // The amount of time the cell can live for when it is infected.
    double infectedLifespan;

    // The amount of time the cell has been infected for.
    int timeInfected;
  
    // The time which has passed since the last division attempt of the cell.
    int timeSinceLastDivision;
    
    // The rate at which the cell divides - the required time between two division attempts of a cell.
    double divisionRate;

    // The type of modification that the cell wants to do to a neighbouring cell. (It could divide into it, infect it, or not modify it)
    NeighbouringCellModificationType modificationToNeighbCell;

    // The id of the  neighbouring epithelial cell which this cell wants to modify.
    repast::AgentId neighbouringCellToModify;

    // Holds a default agent id value, which will be the value of neighbouringCellToModify, when no modification is required.
    repast::AgentId idForNoNeighbourModification;

    // The amount of time which needs to pass after the cell gets infected before it starts releasing new viruses/infecting neighbouring cells.
    double releaseDelay;

    // The amount of time which needs to pass after the cell gets infected before it starts displaying virus proteins on its surface 
    // (before it changes its externalState to DisplayingViralPeptide)
    double displayVirProteinsDelay;

    // Probability of releasing a new virus particle in the extracellular space, when the cell starts producing the progeny virus
    double extracellularReleaseProb;

    // Probability of directly infecting a neighbouring cell as a form of new virus particle release
    double cellToCellTransmissionProb;

    // The virus count rate which a virus producing cell releases each hour.
    double virionReleaseRate;

    // The count of viruses that a virus producing cell will release at the current step.
    int countOfVirionsToRelease;

    // The remainder of virus particles from the release rate (any decimal - if the release rate is 1.3, we'll keep 0.3 as remainder, so when it sums over time
    // it will result in additional virus particles to be released)
    double virionReleaseRemainder;
};
#endif