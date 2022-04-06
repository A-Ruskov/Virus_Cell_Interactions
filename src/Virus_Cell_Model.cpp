/* VirusCellModel.cpp */

#include <stdio.h>
#include <vector>
#include <boost/mpi.hpp>
#include "repast_hpc/AgentId.h"
#include "repast_hpc/RepastProcess.h"
#include "repast_hpc/Utilities.h"
#include "repast_hpc/Properties.h"
#include "repast_hpc/initialize_random.h"
#include "repast_hpc/Point.h"
#include "repast_hpc/SVDataSetBuilder.h"

#include "Virus_Cell_Model.h"

/**********************
*   VirusCellModel::VirusCellModel - Constructor for the VirusCellModel class.
**********************/
VirusCellModel::VirusCellModel(std::string propsFile, int argc, char** argv, boost::mpi::communicator* comm):
context(comm)
{
    // Read in the passed properties of the simmulation. Properties are the parameter values passed.
    props = new repast::Properties(propsFile, argc, argv, comm);

    // Get the index of the final timestep of the simulation.
    stopAt = repast::strToInt(props->getProperty("stop.at"));    

    // Get the counts of agents which are to be initially created for each process.
    countOfVirionAgents = repast::strToInt(props->getProperty("count.of.virions"));
    countOfInnateImmuneCellAgents = repast::strToInt(props->getProperty("count.of.innate.immune.cells"));
    countOfSpecialisedImmuneCellAgents = repast::strToInt(props->getProperty("count.of.specialised.immune.cells"));

    // Since we will be creating new virion agents, innate and specialised immune cell agents we need to track the last id index which was used.
    // That is to be incremented for each new agent of any of those types, to ensure that any new agents will have a unique id.
    currVirionAgentId = 0;
    currInnateImmuneCellAgendId = 0;
    currSpecialisedImmuneCellAgentId = 0;

    // Epithelial cell agents parameters read.
    epithCellAvgLifespan = repast::strToInt(props->getProperty("epithelial.cell.average.lifespan"));
    epithCellLifespanStdev = repast::strToInt(props->getProperty("epithelial.cell.lifespan.standard.dev"));
    epithCellInfectedLifespanAvg = repast::strToInt(props->getProperty("epithelial.cell.infected.lifespan.average"));
    epithCellInfectedLifespanStdev = repast::strToInt(props->getProperty("epithelial.cell.infected.lifespan.standard.dev"));
    epithCellDivisionRateAvg = repast::strToInt(props->getProperty("epithelial.cell.division.rate.average"));
    epithCellDivisionRateStdev = repast::strToInt(props->getProperty("epithelial.cell.division.rate.standard.dev"));
    epithCellVirionReleaseDelayAvg = repast::strToDouble(props->getProperty("epithelial.cell.virion.release.delay.average"));
    epithCellVirionReleaseDelayStdev = repast::strToDouble(props->getProperty("epithelial.cell.virion.release.delay.standard.dev"));
    epithCellDispViralPeptidesDelayAvg = repast::strToDouble(props->getProperty("epithelial.cell.display.viral.peptides.delay.average"));
    epithCellDispViralPeptidesDelayStdev = repast::strToDouble(props->getProperty("epithelial.cell.display.viral.peptides.delay.standard.dev"));
    extracellularVirusReleaseProb = repast::strToDouble(props->getProperty("release.virus.in.extracellular.space.probability"));
    cellToCellTransmissionProb = repast::strToDouble(props->getProperty("cell.to.cell.transmission.probability"));
    epithCellVirionReleaseRateAvg = repast::strToDouble(props->getProperty("epithelial.cell.infected.virion.release.rate.average"));
    epithCellVirionReleaseRateStdev = repast::strToDouble(props->getProperty("epithelial.cell.infected.virion.release.rate.standard.dev"));

    // Virion (Virus Particle) agents parameters read.
    virionAvgLifespan = repast::strToDouble(props->getProperty("virion.average.lifespan"));
    virionLifespanStdev = repast::strToDouble(props->getProperty("virion.lifespan.standard.dev"));
    virionPenetrationProbability = repast::strToDouble(props->getProperty("virion.cell.penetration.probability"));
    virionClearanceProbability = repast::strToDouble(props->getProperty("virion.clearance.probability"));
    virionClearanceProbabilityScaler = repast::strToDouble(props->getProperty("virion.clearance.scaler"));

    // Innate immune cell agents parameters read.
    innateImmuneCellAvgLifespan = repast::strToDouble(props->getProperty("innate.immune.cell.average.lifespan"));
    innateImmuneCellLifespanStdev = repast::strToDouble(props->getProperty("innate.immune.cell.lifespan.stdev"));
    innateImmuneCellInfectedCellRecognitionProb = repast::strToDouble(props->getProperty("innate.immune.cell.infected.cell.recognition.probability"));
    innateImmuneCellInfectedCellEliminationProb = repast::strToDouble(props->getProperty("innate.immune.cell.infected.cell.elimination.probability"));
    innateImmuneCellRecruitSpecImmuneCellProb = repast::strToDouble(props->getProperty("innate.immune.cell.recruit.specialised.immune.cell.probability"));
    innateImmuneCellRecruitRateOfInnateCell = repast::strToDouble(props->getProperty("innate.immune.cell.recruit.rate.of.innate.cell"));
    // This is the rate of specialised immune cells which an innate immune cell recruits per detection of infected epithelial cell agent.
    specialisedImmuneCellRecruitRateOfInnateCell = repast::strToDouble(props->getProperty("specialised.immune.cell.recruit.rate.of.innate.cell"));

    // Specialised immune cell agents parameters read.
    specialisedImmuneCellAvgLifespan = repast::strToDouble(props->getProperty("specialised.immune.cell.average.lifespan"));
    specialisedImmuneCellLifespanStdev = repast::strToDouble(props->getProperty("specialised.immune.cell.lifespan.stdev"));
    specialisedImmuneCellInfectedCellRecognitionProb = repast::strToDouble(props->getProperty("specialised.immune.cell.infected.cell.recognition.probability"));
    specialisedImmuneCellInfectedCellEliminationProb = repast::strToDouble(props->getProperty("specialised.immune.cell.infected.cell.elimination.probability"));
    specialisedImmuneCellRecruitRateOfSpecCell = repast::strToDouble(props->getProperty("specialised.immune.cell.recruit.rate.of.specialised.cell"));

    // Initialize the random singleton with the distributions and random seed provided in the properties.
    initializeRandom(*props, comm);

    if(repast::RepastProcess::instance()->rank() == 1)
    {
        props->writeToSVFile("./output/simulation_parameters_record.csv");
    }    


    // Take the grid dimension sizes, number of dimensions and create the grid projection which will be inhabited by the agents.
    gridDimensionSize = repast::strToInt(props->getProperty("grid.dimension"));
    int originCoordinate = -gridDimensionSize / 2;
    repast::Point<double> origin(originCoordinate, originCoordinate);
    repast::Point<double> extent(gridDimensionSize, gridDimensionSize);

    repast::GridDimensions gridDimensions(origin, extent);

    // Get the parameters about how we want to split the grid across the parallel processes
    int processesCountXAxis = repast::strToInt(props->getProperty("count.of.processes.X.axis"));
    int processesCountYAxis = repast::strToInt(props->getProperty("count.of.processes.Y.axis"));
    std::vector<int> processDimensions;
    processDimensions.push_back(processesCountXAxis);
    processDimensions.push_back(processesCountYAxis);

    // The grid projection will contain agents of type VirusCellInteractionAgents (the parent Agent class), so that it can facilitate all agents types
    // Then we can use the agent type identifier in each agent ID, to cast them to the correct type of agent.
    discreteGridSpace = new repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents>>("AgentsDeiscreteSpace", gridDimensions, processDimensions, 1, comm);

    std::cout << "RANK " << repast::RepastProcess::instance()->rank() << " BOUNDS: " << discreteGridSpace->dimensions().origin() << " " << discreteGridSpace->dimensions().extents() << std::endl;
    
    // Add the grid to the shared context.
   	context.addProjection(discreteGridSpace);


    // Create the agents' package providers and receivers which will be used for agent synchronisation across processes.
    agentProvider = new VirusCellInteractionAgentsPackageProvider(&context);
	agentReceiver = new VirusCellInteractionAgentsPackageReceiver(&context);


    // Initialise Data collection
	// Create the data set builder
	std::string fileOutputName("./output/agents_data.csv");
	repast::SVDataSetBuilder dataBuilder(fileOutputName.c_str(), ",", repast::RepastProcess::instance()->getScheduleRunner().schedule());
	
	// Create the individual data sets to be added to the builder
	DataSource_EpithelialCellsCount* aliveEpithCellsCount_DataSource = new DataSource_EpithelialCellsCount(&context);
	dataBuilder.addDataSource(createSVDataSource("# Alive Epithelial Cells", aliveEpithCellsCount_DataSource, std::plus<int>()));

    DataSource_InfectedEpithelialCellsCount* infectedEpithelialCellsCount_DataSource = new DataSource_InfectedEpithelialCellsCount(&context);
    dataBuilder.addDataSource(createSVDataSource("# Infected Epithelial Cells", infectedEpithelialCellsCount_DataSource, std::plus<int>()));
    
    DataSource_DeadEpithelialCellsCount* deadEpithelialCellsCount_DataSource = new DataSource_DeadEpithelialCellsCount(&context);
    dataBuilder.addDataSource(createSVDataSource("# Dead Epithelial Cells", deadEpithelialCellsCount_DataSource, std::plus<int>()));

	DataSource_VirionsCount* virionsCount_DataSource = new DataSource_VirionsCount(&context);
	dataBuilder.addDataSource(createSVDataSource("# Free Virions", virionsCount_DataSource, std::plus<int>()));

    DataSource_InnateImmuneCellsCount* innateImmuneCellsCount_DataSource = new DataSource_InnateImmuneCellsCount(&context);
    dataBuilder.addDataSource(createSVDataSource("# Innate Immune Cells", innateImmuneCellsCount_DataSource, std::plus<int>()));

    DataSource_SpecialisedImmuneCellsCount* specialisedImmuneCellsCount_DataSource = new DataSource_SpecialisedImmuneCellsCount(&context);
    dataBuilder.addDataSource(createSVDataSource("# Specialised Immune Cells", specialisedImmuneCellsCount_DataSource, std::plus<int>()));

    DataSource_TotalAgentsCount* totalAgentsCount_DataSource = new DataSource_TotalAgentsCount(&context);
    dataBuilder.addDataSource(createSVDataSource("# Agents In Total", totalAgentsCount_DataSource, std::plus<int>()));

	// Use the builder to create the data set
	agentsData = dataBuilder.createDataSet();   
}



/**********************
*   VirusCellModel::~VirusCellModel - Destructor for the VirusCellModel class.
**********************/
VirusCellModel::~VirusCellModel(){
		delete props;
        delete agentProvider;
        delete agentReceiver;

        // Deleting the dataset, will also automatically delete all individual datasets/datasources
        delete agentsData;
}



/**********************
*   VirusCellModel::initSchedule - Function which schedules all events that need to be run during the simulation.
**********************/
void VirusCellModel::initSchedule(repast::ScheduleRunner& runner)
{
	runner.scheduleStop(stopAt);

    // Schedule agents to act on each tick.
    runner.scheduleEvent(1, 1, repast::Schedule::FunctorPtr(new repast::MethodFunctor<VirusCellModel> (this, &VirusCellModel::executeTimestep)));

    // Schedule Data collection. Record data on each tick and write all recorded records on every 3 ticks.
	runner.scheduleEvent(0.1, 1, repast::Schedule::FunctorPtr(new repast::MethodFunctor<repast::DataSet>(agentsData, &repast::DataSet::record)));
	runner.scheduleEvent(0.2, 3, repast::Schedule::FunctorPtr(new repast::MethodFunctor<repast::DataSet>(agentsData, &repast::DataSet::write)));

    runner.scheduleEvent(1.3, 1, repast::Schedule::FunctorPtr(new repast::MethodFunctor<VirusCellModel> (this, &VirusCellModel::printEndOfTimestep)));

	runner.scheduleEndEvent(repast::Schedule::FunctorPtr(new repast::MethodFunctor<repast::DataSet>(agentsData, &repast::DataSet::write)));
}



/**********************
*   VirusCellModel::Init - Initialiser function for the model.
**********************/
void VirusCellModel::init()
{
    int rank = repast::RepastProcess::instance()->rank();

    // Create the epithelial cell agents in the model
    int epithelialCellIndex = 0;
    for( int x = 0; x < discreteGridSpace->dimensions().extents().getX(); ++x)
    {
        for( int y = 0; y < discreteGridSpace->dimensions().extents().getY(); ++y)
        {
            initialiseEpithelialCellAgent(epithelialCellIndex, x, y, false, nullptr);
            ++epithelialCellIndex;
        }
    }

    // Create the initial virion agents in the model
    for( int i = 0; i < countOfVirionAgents; ++i )
    {
        initialiseVirionAgent(i, false, -1, -1);
        currVirionAgentId++;
    }

    // Create the initial innate immune cell agents in the model
    for( int i = 0; i < countOfInnateImmuneCellAgents; ++i )
    {
        initialiseInnateImmuneCellAgent(i, false);
        currInnateImmuneCellAgendId++;
    }

    // Create the initial specialised immune cell agents in the model
    for( int i = 0; i < countOfSpecialisedImmuneCellAgents; ++i )
    {
        initialiseSpecialisedImmuneCellAgent(i, false);   
        currSpecialisedImmuneCellAgentId++; 
    }
}



/**********************
*   VirusCellModel::initialiseEpithelialCellAgent - Initialises an epithelial cell agent and places it on the grid.
**********************/
void VirusCellModel::initialiseEpithelialCellAgent( int epithelialCellIndex, int xCoor, int yCoor, bool isExistingAgentObject, EpithelialCellAgent* theExistingCellObject )
{
    int rank = repast::RepastProcess::instance()->rank();

    // Allocate an arbitrary lifespan to the cell
    repast::NormalGenerator lifespanGenerator = repast::Random::instance()->createNormalGenerator(epithCellAvgLifespan, epithCellLifespanStdev);
    int cellLifespan = 0;
    while(cellLifespan < 1)
    {
        cellLifespan = lifespanGenerator.next();
    }         


    int cellAge = 0;
    // Alloacte an arbitrary age to the cell if there is not an already existing object. If there is no existing object, that means we are just creating
    // the objects of all agents at the start of the simulation. Then we have arbitrary ages for each cell, rather than 0.
    // The age will remain 0, if there is an existing agent object, which means we are handling division. That will be a completely fresh cell with 0 age.
    if(!isExistingAgentObject)
    {
        repast::IntUniformGenerator ageGenerator = repast::Random::instance()->createUniIntGenerator(0, cellLifespan);
        cellAge = ageGenerator.next();
    }
   

    // Alloacte an arbitrary infected cell lifespan. This is how long can the cell live when it is infected with a virion.
    repast::NormalGenerator infectedCellLifespanGen = repast::Random::instance()->createNormalGenerator(epithCellInfectedLifespanAvg, epithCellInfectedLifespanStdev);
    int infectedCellLifespan = 0;
    while( infectedCellLifespan < 1 )
    {
        infectedCellLifespan =  infectedCellLifespanGen.next();
    }

    // The new cell is not infected so reset the infectedTime to 0
    int infectedTime = 0;

    // Allocate an arbitrary division rate
    repast::NormalGenerator divisionRateGen = repast::Random::instance()->createNormalGenerator(epithCellDivisionRateAvg, epithCellDivisionRateStdev);
    int divisionRate = 0;
    while(divisionRate < 1)
    {
        divisionRate = divisionRateGen.next();
    }

   
    // Allocate an arbitrary time since the last division if there is no already existing object. If there is no existing object, that means we are just creating
    // the objects of all agents at the start of the simulation. Then we have times since their last division rather than 0.
    // The time since last division will remain 0, if there is an existing agent object, which means we are handling a newly divided cell. 
    // That will be a completely fresh cell with 0 time since last division.
    int timeSinceLastDivision = 0;
    if(!isExistingAgentObject)
    {
        repast::IntUniformGenerator timeSinceLastDivisionGen = repast::Random::instance()->createUniIntGenerator(0, divisionRate);
        timeSinceLastDivision = timeSinceLastDivisionGen.next();
    }

    // Allocate ana arbitrary delay time, until the cell starts displaying viral proteins (Starts looking infected)
    repast::NormalGenerator displayVirProtDelayGen = repast::Random::instance()->createNormalGenerator(epithCellDispViralPeptidesDelayAvg, epithCellDispViralPeptidesDelayStdev);
    double displayVirProtDelay = displayVirProtDelayGen.next();
    while( displayVirProtDelay < 1 )
    {
        displayVirProtDelay = displayVirProtDelayGen.next();
    }

    // Allocate an arbitrary release delay. The release delay should be larger than the "displayVirProtDelay".
    repast::NormalGenerator releaseDelayGen = repast::Random::instance()->createNormalGenerator(epithCellVirionReleaseDelayAvg, epithCellVirionReleaseDelayStdev);
    double releaseDelay = releaseDelayGen.next();
    while(releaseDelay < 1 && releaseDelay <= displayVirProtDelay)
    {
        releaseDelay = releaseDelayGen.next();
    }


    // Allocate an arbitrary release rate. That is how many new viruses it would produce on a timestep, taken it releases them in the grid.
    repast::NormalGenerator virionReleaseRateGen = repast::Random::instance()->createNormalGenerator(epithCellVirionReleaseRateAvg, epithCellVirionReleaseRateStdev);
    double releaseRate = virionReleaseRateGen.next();
    while(releaseRate < 0.1 )
    {
        releaseRate = virionReleaseRateGen.next();
    }

    // Create the new Epithelial cell agent
    if( isExistingAgentObject == false )
    {
        repast::AgentId newAgentId(epithelialCellIndex, rank, 0);
        newAgentId.currentRank(rank);
        EpithelialCellAgent* newEpithelialCell = new EpithelialCellAgent(newAgentId, cellLifespan, cellAge, infectedCellLifespan, divisionRate, timeSinceLastDivision, releaseDelay, displayVirProtDelay, extracellularVirusReleaseProb, cellToCellTransmissionProb, releaseRate);
        context.addAgent(newEpithelialCell);

        // Set the new cell location in the section of the grid handled by this rank (process). The exact coordinates are passed as an attribute to the function.
        repast::Point<int> agentLocation(xCoor + discreteGridSpace->dimensions().origin().getX() , yCoor + discreteGridSpace->dimensions().origin().getY() );
        discreteGridSpace->moveTo(newAgentId, agentLocation);
    }
    // If the epithelial cell agent object already exists, then we are handling division. In such case, that means, that we just need to reinitialise the
    // existing agent object with the newly generated parameters. Then the existing object will be representing a new cell.
    else if( theExistingCellObject != nullptr )
    {
        // Set the current count of virions to release to 0 and the remainder to 0, as this is a new healthy cell, which has no virions to release
        int countOfVirionsToRelease = 0;
        double virionReleaseRemainder = 0.0;

        theExistingCellObject->set(theExistingCellObject->getId().currentRank(), cellLifespan, cellAge, 
                            EpithelialCellAgent::InternalState::Healthy, EpithelialCellAgent::ExternalState::SeeminglyHealthy, 
                            infectedCellLifespan, infectedTime, divisionRate,  timeSinceLastDivision, releaseDelay, displayVirProtDelay, 
                            EpithelialCellAgent::NeighbouringCellModificationType::NoModification, repast::AgentId(-1, -1, -1, -1), 
                            extracellularVirusReleaseProb, cellToCellTransmissionProb,
                            releaseRate, countOfVirionsToRelease, virionReleaseRemainder);
    }
}



/**********************
*   VirusCellModel::initialiseVirionAgent - Creates a virion (virus particle) agent, sets all its parameters and places it on the grid.
**********************/
void VirusCellModel::initialiseVirionAgent(int virionIndex, bool isAReleasedVirus, int xCoor, int yCoor)
{  
    int rank = repast::RepastProcess::instance()->rank();

    repast::AgentId newVirionId(virionIndex, rank, 1);
    newVirionId.currentRank(rank);

    // Assign an arbitrary lifespan to the virion
    repast::NormalGenerator lifespanGenerator = repast::Random::instance()->createNormalGenerator(virionAvgLifespan, virionLifespanStdev);
    int virionLifespan = lifespanGenerator.next();
    while(virionLifespan < 1)
    {
        virionLifespan = lifespanGenerator.next();
    }

    // Assign an arbitrary age to the virion if it is from the starting population (simulation initialisation). If it is a released virus, then the age should be 0.
    int virionAge = 0;
    if( !isAReleasedVirus )
    {
        repast::IntUniformGenerator ageGenerator = repast::Random::instance()->createUniIntGenerator(0, virionLifespan);
        virionAge = ageGenerator.next();
    }

    // Assign the target cell penetration probability to the virion
    double virionPenetrationProb = virionPenetrationProbability;

    // Assign the virion clearance probability (that is probability that a free virus particle is cleared off through unmodelled immune mechanisms)
    double clearanceProb = virionClearanceProbability;

    // The coefficient of increasing the clearance probability, for each immune cell which is at the current location.
    // The more immune cells are present at that locaiton, the higher the probability that there are also 
    // more unmodelled immune cell mechanisms carried out there, hence higher clearance probability.
    double clearanceProbScaler = virionClearanceProbabilityScaler;

    VirionAgent* newVirion = new VirionAgent(newVirionId, virionLifespan, virionAge, virionPenetrationProb, clearanceProb, clearanceProbScaler);
    context.addAgent(newVirion);

    // Place the agent in the grid spatial projection. If it is a released virus then use the provided coordinates, to place it in the grid.
    // That is since, the released virions are initialised at the position of the cell which has released them
    if( isAReleasedVirus )
    {
        repast::Point<int> virionLocation(xCoor, yCoor);
        discreteGridSpace->moveTo(newVirionId, virionLocation);
    }
    else
    {
        // If it is a virion agent from the starting population, then set its position stochastically.
        // This will place the agent somewhere in the bounds of the part of the grid handled by this process/rank.
        repast::IntUniformGenerator gridXCoorGenerator = repast::Random::instance()->createUniIntGenerator(0, discreteGridSpace->dimensions().extents().getX() - 1);
        int virionXCoor = discreteGridSpace->dimensions().origin().getX() + gridXCoorGenerator.next();

        repast::IntUniformGenerator gridYCoorGenerator = repast::Random::instance()->createUniIntGenerator(0, discreteGridSpace->dimensions().extents().getY() - 1);
        int virionYCoor = discreteGridSpace->dimensions().origin().getY() + gridYCoorGenerator.next();

        repast::Point<int> virionLocation(virionXCoor, virionYCoor);
        discreteGridSpace->moveTo(newVirionId, virionLocation);
    }

}



/**********************
*   VirusCellModel::initialiseInnateImmuneCellAgent - Creates an innate immune cell agent, sets all its parameters and places it on the grid.
**********************/
void VirusCellModel::initialiseInnateImmuneCellAgent( int immuneCellId, bool isRecruitedCell )
{
    int rank = repast::RepastProcess::instance()->rank();

    repast::AgentId newInnateImmuneCellId(immuneCellId, rank, 2);
    newInnateImmuneCellId.currentRank(rank);

    // Stochastically assign the lifepsan of the innate immune cell based on the provided lifespan parameters.
    repast::NormalGenerator lifespanGenerator = repast::Random::instance()->createNormalGenerator(innateImmuneCellAvgLifespan, innateImmuneCellLifespanStdev);
    int innateImmuneCellLifespan = lifespanGenerator.next();
    while( innateImmuneCellLifespan < 1 )
    {
        innateImmuneCellLifespan = lifespanGenerator.next();
    }
  
    // Assign an arbitrary age (between 0 and the lifespan) to the innate immune cell if it is from the starting population (simulation initialisation).
    // If it is a newly recruited cell, then the age should be 0.
    int cellAge = 0;
    if( !isRecruitedCell )
    {
        repast::IntUniformGenerator ageGenerator = repast::Random::instance()->createUniIntGenerator(0, innateImmuneCellLifespan);
        cellAge = ageGenerator.next();
    }

    // Set the probabilities for recognising and eliminating infected cells.
    double infectedCellRecognitionProb = innateImmuneCellInfectedCellRecognitionProb;
    double infectedCellEliminationProb = innateImmuneCellInfectedCellEliminationProb;

    // Set the probability of an innate cell recruiting a specialised immune cell when infection is detected.
    double specialisedImmuneCellRecruitProb = innateImmuneCellRecruitSpecImmuneCellProb;

    // Create the agent.
    InnateImmuneCellAgent* newInnateImmuneCell = new InnateImmuneCellAgent(newInnateImmuneCellId, innateImmuneCellLifespan, cellAge ,infectedCellRecognitionProb, 
    infectedCellEliminationProb, specialisedImmuneCellRecruitProb, innateImmuneCellRecruitRateOfInnateCell, specialisedImmuneCellRecruitRateOfInnateCell);
    context.addAgent(newInnateImmuneCell);

    // Place the agent in the grid spatial projection stochastically.
    // This will place the agent somewhere in the bounds of the part of the grid handled by this process/rank.
    repast::IntUniformGenerator gridXCoorGenerator = repast::Random::instance()->createUniIntGenerator(0, discreteGridSpace->dimensions().extents().getX() - 1);
    int innateImmuneCellX = discreteGridSpace->dimensions().origin().getX() + gridXCoorGenerator.next();

    repast::IntUniformGenerator gridYCoorGenerator = repast::Random::instance()->createUniIntGenerator(0, discreteGridSpace->dimensions().extents().getY() - 1);
    int innateImmuneCellY = discreteGridSpace->dimensions().origin().getY() + gridYCoorGenerator.next();

    repast::Point<int> immuneCellLoc(innateImmuneCellX, innateImmuneCellY);
    discreteGridSpace->moveTo(newInnateImmuneCellId, immuneCellLoc);
}



/**********************
*   VirusCellModel::initialiseSpecialisedImmuneCellAgent - Creates a specialised immune cell agent, sets all its parameters and places it on the grid.
**********************/
void VirusCellModel::initialiseSpecialisedImmuneCellAgent( int immuneCellId, bool isRecruitedCell )
{
    int rank = repast::RepastProcess::instance()->rank();

    repast::AgentId newSpecialisedImmuneCellId(immuneCellId, rank, 3);
    newSpecialisedImmuneCellId.currentRank(rank);

    // Stochastically assign the lifepsan of the specialised immune cell based on the provided lifespan parameters.
    repast::NormalGenerator lifespanGenerator = repast::Random::instance()->createNormalGenerator(specialisedImmuneCellAvgLifespan, specialisedImmuneCellLifespanStdev);
    int specialisedImmuneCellLifespan = lifespanGenerator.next();
    while( specialisedImmuneCellLifespan < 1 )
    {
        specialisedImmuneCellLifespan = lifespanGenerator.next();
    }

    // Assign an arbitrary age (between 0 and the lifespan) to the specialised immune cell if it is from the starting population (simulation initialisation).
    // If it is a newly recruited cell, then the age should be 0.
    int cellAge = 0;
    if( !isRecruitedCell )
    {
        repast::IntUniformGenerator ageGenerator = repast::Random::instance()->createUniIntGenerator(0, specialisedImmuneCellLifespan);
        ageGenerator.next();
    }

    // Set the probabilities for recognising and eliminating infected cells.
    double infectedCellRecognitionProb = specialisedImmuneCellInfectedCellRecognitionProb;
    double infectedCellEliminationProb = specialisedImmuneCellInfectedCellEliminationProb;

    // Create the agent object.
    SpecialisedImmuneCellAgent* newSpecialisedImmuneCell = new SpecialisedImmuneCellAgent(newSpecialisedImmuneCellId, specialisedImmuneCellLifespan, cellAge, infectedCellRecognitionProb, infectedCellEliminationProb, specialisedImmuneCellRecruitRateOfSpecCell);
    context.addAgent(newSpecialisedImmuneCell);

    // Place the agent in the grid spatial projection stochastically. 
    // This will place the agent somewhere in the bounds of the part of the grid handled by this process/rank.
    repast::IntUniformGenerator gridXCoorGenerator = repast::Random::instance()->createUniIntGenerator(0, discreteGridSpace->dimensions().extents().getX() - 1);
    int specialisedImmuneCellX = discreteGridSpace->dimensions().origin().getX() + gridXCoorGenerator.next();

    repast::IntUniformGenerator gridYCoorGenerator = repast::Random::instance()->createUniIntGenerator(0, discreteGridSpace->dimensions().extents().getX() - 1);
    int specialisedImmuneCellY = discreteGridSpace->dimensions().origin().getY() + gridYCoorGenerator.next();

    repast::Point<int> immuneCellLoc(specialisedImmuneCellX, specialisedImmuneCellY);
    discreteGridSpace->moveTo(newSpecialisedImmuneCellId, immuneCellLoc);
}



/**********************
*   VirusCellModel::printEndOfTimestep - Prints a statement that a timestep has finished.
**********************/
void VirusCellModel::printEndOfTimestep(){
    // Get the rank
	int rank = repast::RepastProcess::instance()->rank();
    // Print it only for one rank, as otherwise we'll get duplication of print statements.
    if(rank == 0)
    {
        std::cout<<"Timestep: "<<(int)repast::RepastProcess::instance()->getScheduleRunner().currentTick()<<std::endl;
    }
}



/**********************
*   VirusCellModel::executeTimestep - Function which will execute the timestep. It will trigger all agents to act on each timestep.
**********************/
void VirusCellModel::executeTimestep()
{
    std::vector<VirusCellInteractionAgents*> theBufferZoneAgents;
    context.selectAgents(repast::SharedContext<VirusCellInteractionAgents>::NON_LOCAL, theBufferZoneAgents, 0, false);

    std::vector<VirusCellInteractionAgents*>::iterator iter;

    // First we need to check if the  agents in the bufferzones have requested division/infection to agents which are local to this rank.
    // This is the only way to propagate changes to the original agents. 
    // Otherwise, these modification would have been done to the copies of the agents rather than the originals and the changes would not reach the originals.
    // In this way we ensure that the agents will act in an environmen where all agents are at their most up-to date state.
    for( iter = theBufferZoneAgents.begin(); iter != theBufferZoneAgents.end(); ++iter )
    {
        if( (*iter)->getId().agentType() == 0 )
        {
            checkForCellDivision(*iter);
            checkForCellToCellInfection(*iter);
        }
    }


    std::vector<VirusCellInteractionAgents*> theLocalAgents;
    // Get the local agents and make them perform a step. This will include all agents: Epithelial cells, Virions, Specialised and Non-Specialised immune cells.
    // They will also be in random order which will provide the stochasticity we need, as there is no way to make them act synchronously.
    context.selectAgents(repast::SharedContext<VirusCellInteractionAgents>::LOCAL, theLocalAgents);

    for( iter = theLocalAgents.begin(); iter != theLocalAgents.end(); ++iter )
    {
        (*iter)->doStep(&context, discreteGridSpace);

        // For each specific type of agent we need if they have requested any change to the environment, which is only handled by the Virus_Cell_Model clas.
        if( (*iter)->getId().agentType() == 0 )
        {
            // Check if the cell has requested division, viral release or infection of a neighbouring cell.
            checkForCellDivision(*iter);
            checkForCellVirionRelease(*iter);
            checkForCellToCellInfection(*iter);
        }
        else if ((*iter)->getId().agentType() == 2 )
        {
            // Check if the innate immune cell has requested recruitment of other immune cells.
            checkForInnateImmuneCellRecruitment(*iter);
            checkForSpecialisedImmuneCellRecruitement(*iter);
        }
        else if( (*iter)->getId().agentType() == 3 )
        {
            // Check if the specialised immune cell has requested recruitment of other specialised immune cells.
            checkForSpecialisedImmuneCellRecruitement(*iter);
        }
        
        
        // Check if the agent has died during the timestep and remove it from the simulation.
        if( (*iter)->getId().agentType() != 0 )
        {
            removeLocalAgentIfDead(*iter);
        }
    } 

    // Balancing the grid will identify the agents which have crossed the boundaries of their rank and need to be moved. 
    discreteGridSpace->balance();

    // Synchronising the agent status will move the agents to the correct process.
    repast::RepastProcess::instance()->synchronizeAgentStatus<VirusCellInteractionAgents, VirusCellInteractionAgentPackage, VirusCellInteractionAgentsPackageProvider, 
        VirusCellInteractionAgentsPackageReceiver>(context, *agentProvider, *agentReceiver, *agentReceiver);

    // Synchronise the data about the section of the whole grid handled by each process. 
    repast::RepastProcess::instance()->synchronizeProjectionInfo<VirusCellInteractionAgents, VirusCellInteractionAgentPackage, VirusCellInteractionAgentsPackageProvider, 
        VirusCellInteractionAgentsPackageReceiver>(context, *agentProvider, *agentReceiver, *agentReceiver);

    // Synchronise all agents which are non-local to this process (The copies of non-local agents which this process owns).
    // Ensures the buffer zone agents are most up-to-date copies of their original agents.
    repast::RepastProcess::instance()->synchronizeAgentStates<VirusCellInteractionAgentPackage, VirusCellInteractionAgentsPackageProvider, 
        VirusCellInteractionAgentsPackageReceiver>(*agentProvider, *agentReceiver);
}



/**********************
*   VirusCellModel::checkForCellDivision - Function which checks if an Epithelial Cell agent is ready to divide
*   If it is, then it check where it will divide, and sets the new epithelial cell with its parameters
*   This is handled by the VirusCellModel opposed to the EpithelialCellAgent class itself in order to allow 
*   epithelial cells to divide into cells which are native to another process. Otherwise division would be applied to a copy of an agent
*   found in the buffer zone and the division would not propagate back to the original agent.
**********************/
void VirusCellModel::checkForCellDivision(VirusCellInteractionAgents* theEpithelialCellAgent)
{
    EpithelialCellAgent* epithelialCell = static_cast<EpithelialCellAgent*>(theEpithelialCellAgent);
    if( epithelialCell == nullptr )
    {
        std::cout<<"Incorrect agent object was provided to VirusCellModel::checkForCellDivision!"<<std::endl;
        return;
    }

    // If the epithelial cell wants to modify a neighbouring cell and the modification it wants to do is to divide into it, 
    // then there needs to be a division into the stored cell
    if( epithelialCell->getTypeOfModifToNeighbCell() == EpithelialCellAgent::NeighbouringCellModificationType::ToDivideInto )
    {
        int rank = repast::RepastProcess::instance()->rank();

        // Get the id of the neighbouring cell where the division will happen.
        repast::AgentId cellToDivideIntoId = epithelialCell->getNeighbouringCellToModify();
        if(cellToDivideIntoId.id() != -1 && cellToDivideIntoId.agentType() == 0 && rank == cellToDivideIntoId.currentRank())
        {
            VirusCellInteractionAgents* cellToBeRevivedBaseClass = context.getAgent(cellToDivideIntoId);  
            EpithelialCellAgent* cellToBeRevived = static_cast<EpithelialCellAgent*>(cellToBeRevivedBaseClass);

            // Dividing, will reset the agent's computational object with new parameters (revives it with new parameters).
            // Ensure we revive only dead cells. Since cross-process agent synchronisation happens at the end of the step,
            // There could potentially be some inconsistencies (local cell trying to revive a cell and a buffer zone cell trying to revive the same cell)
            // These inconsistencies are prevented by checking the cell state and ensuring it is dead.
            if(cellToBeRevived->getInternalState() == EpithelialCellAgent::InternalState::Dead)
            {
                // Reinitialise the epithelial cell agent (Carry the division out).
                initialiseEpithelialCellAgent(-1, -1, -1, true, cellToBeRevived);
            }
        }
    }
}



/**********************
*   VirusCellModel::checkForCellToCellInfection - Function which checks if an Epithelial Cell agent is to infect a neighbouring cell through 
*   Cell-to-Cell virus transmission release. If it is, then it infects the corresponding neighbouring cell.
*   This is handled by the VirusCellModel opposed to the EpithelialCellAgent class itself in order to allow 
*   epithelial cells to infect neighbouring cells which are native to another process. 
*   Otherwise the cell-cell infection would be applied to a copy of an agent
*   found in the buffer zone and the infection would not propagate back to the original agent.
**********************/
void VirusCellModel::checkForCellToCellInfection(VirusCellInteractionAgents* theEpithelialCellAgent)
{
    EpithelialCellAgent* epithelialCell = static_cast<EpithelialCellAgent*>(theEpithelialCellAgent);
    if( epithelialCell == nullptr )
    {
        std::cout<<"Incorrect agent object was provided!"<<std::endl;
        return;
    }

    // If the epithelial cell wants to modify a neighbouring cell and the modification it wants to do is to infect it, 
    if( epithelialCell->getTypeOfModifToNeighbCell() == EpithelialCellAgent::NeighbouringCellModificationType::ToInfect && epithelialCell->getInternalState() == EpithelialCellAgent::InternalState::Infected )
    {
        int rank = repast::RepastProcess::instance()->rank();

        // Get the id of the neighbouring cell which is to be infected.
        repast::AgentId cellToInfectId = epithelialCell->getNeighbouringCellToModify();
        if(cellToInfectId.id() != -1 && cellToInfectId.agentType() == 0 && rank == cellToInfectId.currentRank())
        {
            VirusCellInteractionAgents* cellToInfectBaseClass = context.getAgent(cellToInfectId);  
            EpithelialCellAgent* cellToBeInfected = static_cast<EpithelialCellAgent*>(cellToInfectBaseClass);

            // Ensure we infect only healthy cells. Since cross-process agent synchronisation happens at the end of the step,
            // There could potentially be some inconsistencies (local cell trying to revive a cell and a buffer zone cell trying to infect the same cell)
            // These inconsistencies are prevented by checking the cell state and ensuring it is healthy.
            if(cellToBeInfected->getInternalState() == EpithelialCellAgent::InternalState::Healthy)
            {
                cellToBeInfected->infect();
            }
        }
    }
}



/**********************
*   VirusCellModel::checkForCellVirionRelease - Function which checks if an Epithelial Cell agent has requested the release of new virions in the extracellular space.
*   This is handled by the VirusCellModel opposed to the EpithelialCellAgent class as only the Model is able to create and initialise new agents.
**********************/
void VirusCellModel::checkForCellVirionRelease(VirusCellInteractionAgents* theEpithelialCellAgent)
{
    EpithelialCellAgent* epithelialCell = static_cast<EpithelialCellAgent*>(theEpithelialCellAgent);
    if( epithelialCell == nullptr )
    {
        std::cout<<"The theEpithelialCellAgent passed to VirusCellModel::checkForCellVirionRelease was of incorrect type! Could not cast to EpithelialCellAgent!"<<std::endl;
        return;
    }

    // Get the count of new virion agents to be released.
    int numVirionsToRelease = epithelialCell->getVirionCountToRelease();
    // If there are any new virus particles to be released, the required count of new virion agents will be created. 
    for( int i = 0; i < numVirionsToRelease; ++i)
    {
        // Place the agent in the grid spatial projection at the location of the epithelial cell which releases it.
        std::vector<int> epithelialCellLoc;
        discreteGridSpace->getLocation(epithelialCell->getId(), epithelialCellLoc);

        initialiseVirionAgent(++currVirionAgentId, true, epithelialCellLoc[0], epithelialCellLoc[1]);
    }
}



/**********************
*   VirusCellModel::checkForSpecialisedImmuneCellRecruitement - Function which checks if an immune cell agent (innate/specialised) is recruiting new specialised immune cells.
*   This is handled by the VirusCellModel opposed to the EpithelialCellAgent class as only the Model is able to create and initialise new agents.
**********************/
void VirusCellModel::checkForSpecialisedImmuneCellRecruitement(VirusCellInteractionAgents* theRecruitingImmuneCell)
{
    repast::AgentId theRecruitingCellId = theRecruitingImmuneCell->getId();
    int countOfSpecCellsToRecruit = 0;

    if(theRecruitingCellId.agentType() == 2)
    {
        InnateImmuneCellAgent* theRecruitingInnateImmCell = static_cast<InnateImmuneCellAgent*>(theRecruitingImmuneCell);
        if( theRecruitingInnateImmCell == nullptr )
        {
            std::cout<<"The theRecruitingImmuneCell passed to VirusCellModel::checkForSpecialisedImmuneCellRecruitement was of incorrect type! Could not cast to InnateImmuneCellAgent!"<<std::endl;
            return;
        }

        countOfSpecCellsToRecruit = theRecruitingInnateImmCell->getCountOfSpecialisedCellsToRecruit();
    }
    else if (theRecruitingCellId.agentType() == 3)
    {
        SpecialisedImmuneCellAgent* theRecruitingSpecialisedImmCell = static_cast<SpecialisedImmuneCellAgent*>(theRecruitingImmuneCell);
        if( theRecruitingSpecialisedImmCell == nullptr )
        {
            std::cout<<"The theRecruitingImmuneCell passed to VirusCellModel::checkForSpecialisedImmuneCellRecruitement was of incorrect type! Could not cast to SpecialisedImmuneCellAgent!"<<std::endl;
            return;
        }

        countOfSpecCellsToRecruit = theRecruitingSpecialisedImmCell->getCountOfSpecCellsToRecruit();  
    }

    // If there are any new virus particles to be released, the required count of new virion agents will be created and placed in the grid stochastically.
    for( int i = 0; i < countOfSpecCellsToRecruit; ++i )
    {
        initialiseSpecialisedImmuneCellAgent(currSpecialisedImmuneCellAgentId++, true);
    }            
    
}



/**********************
*   VirusCellModel::checkForInnateImmuneCellRecruitment - Function which checks if an innate immune cell agent is recruiting more innate immune cells.
*   This is handled by the VirusCellModel opposed to the EpithelialCellAgent class as only the Model is able to create and initialise new agents.
**********************/ 
void VirusCellModel::checkForInnateImmuneCellRecruitment(VirusCellInteractionAgents* theRecruitingImmuneCell)
{
    InnateImmuneCellAgent* theRecruitingInnateImmCell = static_cast<InnateImmuneCellAgent*>(theRecruitingImmuneCell);
    if( theRecruitingImmuneCell == nullptr )
    {
        std::cout<<"Incorrect pointer to innate immune cell agent was passed to checkForInnateImmuneCellRecruitment()!"<<std::endl;
        return;
    }

    // If there are any new virus particles to be released, the required count of new virion agents will be created and placed in the grid stochastically.
    int countOfInnateCellsToRecruit = theRecruitingInnateImmCell->getCountOfInnateCellsToRecruit();
    for( int i = 0; i < countOfInnateCellsToRecruit; ++i )
    {
        initialiseInnateImmuneCellAgent(currInnateImmuneCellAgendId++, true);
    }
}



/**********************
*   VirusCellModel::removeLocalAgentIfDead - Function which checks if a Virion/Innate Immune Cell/Specailise Immune Cell agent 
*   has the dead state and removes it from the simulation.
*   We will ignore the Epithelial cell agents, since they are not removed from the simulation if they die.
*   That is since, epithelial cells can divide, and the division of a cell will basically "revive" a dead cell.
*   This will remove the need to create a new epithelial cell agent in the simulation.
**********************/
void VirusCellModel::removeLocalAgentIfDead(VirusCellInteractionAgents* theAgent)
{
    repast::AgentId theAgentId = theAgent->getId();
    int theAgentType = theAgentId.agentType();

    bool removeAgent = false;

    // Check if the agent state is dead and if it is remove it from the process
    if( theAgentType == 1 )
    {
        VirionAgent* theVirion = static_cast<VirionAgent*>(theAgent);
        // For virion agents we will also remove contained viruses (ones which have managed to infect a cell), 
        // as they cannot continue moving and attempting to infect new cells.
        removeAgent = (theVirion->getVirionState() == VirionAgent::VirionStates::Dead || theVirion->getVirionState() == VirionAgent::VirionStates::Contained);
    }
    else if( theAgentType == 2 )
    {
        InnateImmuneCellAgent* theInnateImmuneCell = static_cast<InnateImmuneCellAgent*>(theAgent);
        removeAgent = (theInnateImmuneCell->getCellState() == InnateImmuneCellAgent::InnateImmuneCellStates::Dead);

        // Innate immune cells need to be kepts at a constant minimum level, as they are always present in the organism.
        // If their count drops under the initial count of innate immune cells at the start of the simulation, 
        // then we need to create a new innate immune cell in the place of the dead one.
        if(removeAgent)
        {
            std::vector<VirusCellInteractionAgents*> theLocalInnateImmuneCells;
    
            // Get all local innate immune cell agents
            context.selectAgents(repast::SharedContext<VirusCellInteractionAgents>::LOCAL, theLocalInnateImmuneCells, 2, false);

            // Until we reach the minimum thrshold count of innate immune cells, create new ones,
            while( theLocalInnateImmuneCells.size() <= countOfInnateImmuneCellAgents )
            {
                context.selectAgents(repast::SharedContext<VirusCellInteractionAgents>::LOCAL, theLocalInnateImmuneCells, 2, false);
                initialiseInnateImmuneCellAgent(currInnateImmuneCellAgendId++, true);
            }
        }
    }
    else if( theAgentType == 3 )
    {
        SpecialisedImmuneCellAgent* theSpecialisedImmuneCell = static_cast<SpecialisedImmuneCellAgent*>(theAgent);
        removeAgent = (theSpecialisedImmuneCell->getCellState() == SpecialisedImmuneCellAgent::SpecialisedImmuneCellStates::Dead);
    }

    // Remove the agent from the simulation.
    if( removeAgent )
    {
        repast::RepastProcess::instance()->agentRemoved( theAgentId );
        context.removeAgent( theAgentId );
    }

}