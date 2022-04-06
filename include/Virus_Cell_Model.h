/* Virus_Cell_Model.h */
#ifndef VIRUS_CELL_MODEL
#define VIRUS_CELL_MODEL

/**********************
*   Include files
**********************/
#include <boost/mpi.hpp>
#include "repast_hpc/Schedule.h"
#include "repast_hpc/Properties.h"
#include "repast_hpc/SharedContext.h"
#include "repast_hpc/AgentRequest.h"
#include "repast_hpc/TDataSource.h"
#include "repast_hpc/SVDataSet.h"

// Spatial Projection includes
#include "repast_hpc/SharedDiscreteSpace.h"
#include "repast_hpc/GridComponents.h"

// Include Data Collection File (Containing the datasources)
#include "Data_Collection.h"

// Include agent related files
#include "Virus_Cell_Agent.h"
#include "Epithelial_Cell_Agent.h"
#include "Virion_Agent.h"
#include "Innate_Immune_Cell.h"
#include "Specialised_Immune_Cell.h"

// Include the file which contains all agent package syncrhonisation implementation
#include "Agent_Synchronisation_Package_Pattern.h"



/**********************
*   Forward class declarations
**********************/
class VirusCellInteractionAgentsPackageProvider;
class VirusCellInteractionAgentsPackageReceiver;


/**********************
*   The model class
**********************/
class VirusCellModel{ 
    int countVirsWhichManagedToInfectACell = 0;
	
    int stopAt;
    int countOfVirionAgents;
    int countOfInnateImmuneCellAgents;
    int countOfSpecialisedImmuneCellAgents;
    int gridDimensionSize;

    // Epithelial cell agents parameters
    int epithCellAvgLifespan;
    int epithCellLifespanStdev;
    int epithCellInfectedLifespanAvg;
    int epithCellInfectedLifespanStdev;
    int epithCellDivisionRateAvg;
    int epithCellDivisionRateStdev;
    double epithCellVirionReleaseDelayAvg;
    double epithCellVirionReleaseDelayStdev;
    double epithCellDispViralPeptidesDelayAvg;
    double epithCellDispViralPeptidesDelayStdev;
    double extracellularVirusReleaseProb;
    double cellToCellTransmissionProb;
    double epithCellVirionReleaseRateAvg;
    double epithCellVirionReleaseRateStdev;

    // Virion (Virus Particle) agents parameters.
    double virionAvgLifespan;
    double virionLifespanStdev;
    double virionPenetrationProbability;
    double virionClearanceProbability;
    double virionClearanceProbabilityScaler;

    // Innate immune cell agents parameters.
    double innateImmuneCellAvgLifespan;
    double innateImmuneCellLifespanStdev;
    double innateImmuneCellInfectedCellRecognitionProb;
    double innateImmuneCellInfectedCellEliminationProb;
    double innateImmuneCellRecruitSpecImmuneCellProb;
    double innateImmuneCellRecruitRateOfInnateCell;
    double specialisedImmuneCellRecruitRateOfInnateCell;
    
    // Specialised immune cell agents parameters.
    double specialisedImmuneCellAvgLifespan;
    double specialisedImmuneCellLifespanStdev;
    double specialisedImmuneCellInfectedCellRecognitionProb;
    double specialisedImmuneCellInfectedCellEliminationProb;
    double specialisedImmuneCellRecruitRateOfSpecCell;

    // Trackers of the last used index for each of the agent types. Used in order to ensi
    int currVirionAgentId;
    int currInnateImmuneCellAgendId;
    int currSpecialisedImmuneCellAgentId;
    

	repast::Properties* props;
    repast::SharedContext<VirusCellInteractionAgents> context;

    // Agent package provider and receiver
    VirusCellInteractionAgentsPackageProvider* agentProvider;
	VirusCellInteractionAgentsPackageReceiver* agentReceiver;

    // Data collection variables
    repast::SVDataSet* agentsData;

    // The grid shared by the processes.
    repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace;
public:
	VirusCellModel(std::string propsFile, int argc, char** argv, boost::mpi::communicator* comm);
	~VirusCellModel();
	void init();
	void initSchedule(repast::ScheduleRunner& runner);

private:
    void printEndOfTimestep();
	void executeTimestep();
    void recordResults();

    void initialiseEpithelialCellAgent( int epithelialCellIndex, int xCoor, int yCoor, bool isExistingAgentObject, EpithelialCellAgent* theExistingCellObject );
    void initialiseVirionAgent(int virionIndex, bool isAReleasedVirus, int xCoor, int yCoor);
    void initialiseInnateImmuneCellAgent( int immuneCellId, bool isFreshCell );
    void initialiseSpecialisedImmuneCellAgent( int immuneCellId, bool isFreshCell );

    void checkForCellDivision(VirusCellInteractionAgents* theEpithelialCellAgent);
    void checkForCellToCellInfection(VirusCellInteractionAgents* theEpithelialCellAgent);
    void checkForCellVirionRelease(VirusCellInteractionAgents* theEpithelialCellAgent);
    void checkForSpecialisedImmuneCellRecruitement(VirusCellInteractionAgents* theRecruitingImmuneCell);
    void checkForInnateImmuneCellRecruitment(VirusCellInteractionAgents* theRecruitingImmuneCell);
    void removeLocalAgentIfDead(VirusCellInteractionAgents* theAgent);
};

#endif // #ifndef VIRUS_CELL_MODEL