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

// Include agent related files
#include "Virus_Cell_Agent.h"
#include "Epithelial_Cell_Agent.h"
#include "Virion_Agent.h"
#include "Innate_Immune_Cell.h"
#include "Specialised_Immune_Cell.h"



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
    int countOfEpithelialCellAgents;
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

    // Virion (Virus Particle) agents parameters.
    double virionAvgLifespan;
    double virionLifespanStdev;
    double virionCellPenetrationProbability;
    double virionClearanceProbability;
    double virionClearanceProbabilityScaler;

    // Innate immune cell agents parameters.
    double innateImmuneCellAvgLifespan;
    double innateImmuneCellLifespanStdev;
    double innateImmuneCellInfectedCellRecognitionProb;
    double innateImmuneCellInfectedCellEliminationProb;
    double innateImmuneCellRecruitSpecImmuneCellProb;
    
    // Specialised immune cell agents parameters.
    double specialisedImmuneCellAvgLifespan;
    double specialisedImmuneCellLifespanStdev;
    double specialisedImmuneCellInfectedCellRecognitionProb;
    double specialisedImmuneCellInfectedCellEliminationProb;

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

    // Shared space
    repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents> >* discreteGridSpace;

    std::vector<EpithelialCellAgent*> epithelialCellsToRevive;
public:
	VirusCellModel(std::string propsFile, int argc, char** argv, boost::mpi::communicator* comm);
	~VirusCellModel();
	void init();
    void requestAgents();
	void doSomething();
	void initSchedule(repast::ScheduleRunner& runner);
	void recordResults();

private:
    void initialiseEpithelialCellAgent( int epithelialCellIndex, int xCoor, int yCoor, bool isExistingAgentObject, EpithelialCellAgent* theExistingCellObject );
    void initialiseVirionCellAgent(int virionIndex, bool isAReleasedVirus, int xCoor, int yCoor);
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





/**********************
*   Agent Package Provider class
**********************/
class VirusCellInteractionAgentsPackageProvider {
	
private:
    repast::SharedContext<VirusCellInteractionAgents>* agentsContext;
	
public:
	
    VirusCellInteractionAgentsPackageProvider(repast::SharedContext<VirusCellInteractionAgents>* contextPtr);
	
    void providePackage(VirusCellInteractionAgents * agent, std::vector<VirusCellInteractionAgentPackage>& out);
	
    void provideContent(repast::AgentRequest req, std::vector<VirusCellInteractionAgentPackage>& out);
	
};





/**********************
* Agent Package Receiver class
**********************/
class VirusCellInteractionAgentsPackageReceiver {
	
private:
    repast::SharedContext<VirusCellInteractionAgents>* agentsContext;
	
public:
	
    VirusCellInteractionAgentsPackageReceiver(repast::SharedContext<VirusCellInteractionAgents>* agentPtr);
	
    VirusCellInteractionAgents * createAgent(VirusCellInteractionAgentPackage package);
	
    void updateAgent(VirusCellInteractionAgentPackage package);
	
};



/***********************************
********** DATA COLLECTION *********
***********************************/


/**********************
* Data Source class for tracking the count of alive epithelial cells in the model
**********************/
class DataSource_EpithelialCellsCount : public repast::TDataSource<int>{
private:
	repast::SharedContext<VirusCellInteractionAgents>* context;
    
public:
	DataSource_EpithelialCellsCount(repast::SharedContext<VirusCellInteractionAgents>* theContext);
	int getData();
};


/**********************
* Data Source class for tracking the count of alive epithelial cells in the model
**********************/
class DataSource_VirionsCount : public repast::TDataSource<int>{
private:
	repast::SharedContext<VirusCellInteractionAgents>* context;
    
public:
	DataSource_VirionsCount(repast::SharedContext<VirusCellInteractionAgents>* theContext);
	int getData();
};


/**********************
* Data Source class for tracking the count of innate immune cells in the model
**********************/
class DataSource_InnateImmuneCellsCount : public repast::TDataSource<int>{
private:
	repast::SharedContext<VirusCellInteractionAgents>* context;
    
public:
	DataSource_InnateImmuneCellsCount(repast::SharedContext<VirusCellInteractionAgents>* theContext);
	int getData();
};



/**********************
* Data Source class for tracking the count of specialised immune cells in the model
**********************/
class DataSource_SpecialisedImmuneCellsCount : public repast::TDataSource<int>{
private:
	repast::SharedContext<VirusCellInteractionAgents>* context;
    
public:
	DataSource_SpecialisedImmuneCellsCount(repast::SharedContext<VirusCellInteractionAgents>* theContext);
	int getData();
};



/**********************
* Data Source class for tracking the count of infected epithelial cells in the model
**********************/
class DataSource_InfectedEpithelialCellsCount : public repast::TDataSource<int>
{
private:
	repast::SharedContext<VirusCellInteractionAgents>* context;
    
public:
	DataSource_InfectedEpithelialCellsCount(repast::SharedContext<VirusCellInteractionAgents>* theContext);
	int getData();
};