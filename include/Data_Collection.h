/***********************************
***** DATA COLLECTION CLASSES *****
***********************************/
#include "repast_hpc/TDataSource.h"
#include "repast_hpc/SVDataSet.h"
#include "repast_hpc/SharedContext.h"

// Include agent related files
#include "Virus_Cell_Agent.h"
#include "Epithelial_Cell_Agent.h"
#include "Virion_Agent.h"
#include "Innate_Immune_Cell.h"
#include "Specialised_Immune_Cell.h"


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
* Data Source class for tracking the count of free virion agents in the model
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



/**********************
* Data Source class for tracking the count of dead epithelial cells in the model
**********************/
class DataSource_DeadEpithelialCellsCount : public repast::TDataSource<int>
{
private:
	repast::SharedContext<VirusCellInteractionAgents>* context;
    
public:
	DataSource_DeadEpithelialCellsCount(repast::SharedContext<VirusCellInteractionAgents>* theContext);
	int getData();
};



/**********************
* Data Source class for tracking the total count of agents in the model
**********************/
class DataSource_TotalAgentsCount : public repast::TDataSource<int>
{
private:
	repast::SharedContext<VirusCellInteractionAgents>* context;
    
public:
	DataSource_TotalAgentsCount(repast::SharedContext<VirusCellInteractionAgents>* theContext);
	int getData();
};