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
    props = new repast::Properties(propsFile, argc, argv, comm);
    stopAt = repast::strToInt(props->getProperty("stop.at"));
    countOfEpithelialCellAgents = repast::strToInt(props->getProperty("count.of.epithelial.cells"));
    countOfVirionAgents = repast::strToInt(props->getProperty("count.of.virions"));
    countOfInnateImmuneCellAgents = repast::strToInt(props->getProperty("count.of.innate.immune.cells"));
    countOfSpecialisedImmuneCellAgents = repast::strToInt(props->getProperty("count.of.specialised.immune.cells"));

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
    // Set the probabilities of the 2 virus release mechanisms.
    // An infected cell which has passed the release delay, will follow these probabilities in order to "produce" the new virus particles
    // It could release it in the extracellular space (the grid), or could directly infect a neighbouring cell (known as cell-to-cell transmission release)
    extracellularVirusReleaseProb = repast::strToDouble(props->getProperty("release.virus.in.extracellular.space.probability"));
    cellToCellTransmissionProb = repast::strToDouble(props->getProperty("cell.to.cell.transmission.probability"));

    // Virion (Virus Particle) agents parameters read.
    virionAvgLifespan = repast::strToDouble(props->getProperty("virion.average.lifespan"));
    virionLifespanStdev = repast::strToDouble(props->getProperty("virion.lifespan.standard.dev"));
    virionCellPenetrationProbability = repast::strToDouble(props->getProperty("virion.cell.penetration.probability"));
    virionClearanceProbability = repast::strToDouble(props->getProperty("virion.clearance.probability"));
    virionClearanceProbabilityScaler = repast::strToDouble(props->getProperty("virion.clearance.scaler"));

    // Innate immune cell agents parameters read.
    innateImmuneCellAvgLifespan = repast::strToDouble(props->getProperty("innate.immune.cell.average.lifespan"));
    innateImmuneCellLifespanStdev = repast::strToDouble(props->getProperty("innate.immune.cell.lifespan.stdev"));
    innateImmuneCellInfectedCellRecognitionProb = repast::strToDouble(props->getProperty("innate.immune.cell.infected.cell.recognition.probability"));
    innateImmuneCellInfectedCellEliminationProb = repast::strToDouble(props->getProperty("innate.immune.cell.infected.cell.elimination.probability"));
    innateImmuneCellRecruitSpecImmuneCellProb = repast::strToDouble(props->getProperty("innate.immune.cell.recruit.specialised.immune.cell.probability"));

    // Specialised immune cell agents parameters read.
    specialisedImmuneCellAvgLifespan = repast::strToDouble(props->getProperty("specialised.immune.cell.average.lifespan"));
    specialisedImmuneCellLifespanStdev = repast::strToDouble(props->getProperty("specialised.immune.cell.lifespan.stdev"));
    specialisedImmuneCellInfectedCellRecognitionProb = repast::strToDouble(props->getProperty("specialised.immune.cell.infected.cell.recognition.probability"));
    specialisedImmuneCellInfectedCellEliminationProb = repast::strToDouble(props->getProperty("specialised.immune.cell.infected.cell.elimination.probability"));

    // Initialize the random singleton with the distributions and random seed provided in the properties.
    initializeRandom(*props, comm);

    if(repast::RepastProcess::instance()->rank() == 1)
    {
        props->writeToSVFile("./output/record.csv");
    }

    currVirionAgentId = 0;
    currInnateImmuneCellAgendId = 0;
    currSpecialisedImmuneCellAgentId = 0;
    

    // Create the agents' package providers and receivers.
    agentProvider = new VirusCellInteractionAgentsPackageProvider(&context);
	agentReceiver = new VirusCellInteractionAgentsPackageReceiver(&context);

    // Take the grid dimension sizes, number of dimensions and create the grid projection which will be inhabited by the agents.
    gridDimensionSize = repast::strToInt(props->getProperty("grid.dimension"));
    int originCoordinate = -gridDimensionSize / 2;
    repast::Point<double> origin(originCoordinate, originCoordinate);
    repast::Point<double> extent(gridDimensionSize, gridDimensionSize);

    repast::GridDimensions gridDimensions(origin, extent);

    std::vector<int> processDimensions;
    processDimensions.push_back(8);
    processDimensions.push_back(8);

    // The grid projection will contain agents of type VirusCellInteractionAgents, so that it can facilitate all agents types
    // Then we can use the agent type identifier in each agent ID, to cast them to the correct type of agent.
    discreteGridSpace = new repast::SharedDiscreteSpace<VirusCellInteractionAgents, repast::WrapAroundBorders, repast::SimpleAdder<VirusCellInteractionAgents>>("AgentsDeiscreteSpace", gridDimensions, processDimensions, 1, comm);

    std::cout << "RANK " << repast::RepastProcess::instance()->rank() << " BOUNDS: " << discreteGridSpace->dimensions().origin() << " " << discreteGridSpace->dimensions().extents() << std::endl;
    
   	context.addProjection(discreteGridSpace);


    // Data collection
	// Create the data set builder
	std::string fileOutputName("./output/agents_data.csv");
	repast::SVDataSetBuilder dataBuilder(fileOutputName.c_str(), ",", repast::RepastProcess::instance()->getScheduleRunner().schedule());
	
	// Create the individual data sets to be added to the builder
	DataSource_EpithelialCellsCount* aliveEpithCellsCount_DataSource = new DataSource_EpithelialCellsCount(&context);
	dataBuilder.addDataSource(createSVDataSource("# Alive Epithelial Cells", aliveEpithCellsCount_DataSource, std::plus<int>()));

    DataSource_InfectedEpithelialCellsCount* infectedEpithelialCellsCount_DataSource = new DataSource_InfectedEpithelialCellsCount(&context);
    dataBuilder.addDataSource(createSVDataSource("# Infected Epithelial Cells", infectedEpithelialCellsCount_DataSource, std::plus<int>()));
    
	DataSource_VirionsCount* virionsCount_DataSource = new DataSource_VirionsCount(&context);
	dataBuilder.addDataSource(createSVDataSource("# Virions", virionsCount_DataSource, std::plus<int>()));

    DataSource_InnateImmuneCellsCount* innateImmuneCellsCount_DataSource = new DataSource_InnateImmuneCellsCount(&context);
    dataBuilder.addDataSource(createSVDataSource("# Innate Immune Cells", innateImmuneCellsCount_DataSource, std::plus<int>()));

    DataSource_SpecialisedImmuneCellsCount* specialisedImmuneCellsCount_DataSource = new DataSource_SpecialisedImmuneCellsCount(&context);
    dataBuilder.addDataSource(createSVDataSource("# Specialised Immune Cells", specialisedImmuneCellsCount_DataSource, std::plus<int>()));

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

        delete agentsData;
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


    // Create the virion agents in the model
    for( int i = 0; i < countOfVirionAgents; ++i )
    {
        initialiseVirionCellAgent(i, false, -1, -1);
        currVirionAgentId++;
    }


    // Create the innate immune cell agents in the model
    for( int i = 0; i < countOfInnateImmuneCellAgents; ++i )
    {
        initialiseInnateImmuneCellAgent(i, false);
        currInnateImmuneCellAgendId++;
    }

    // Create the specialised immune cell agents in the model
    for( int i = 0; i < countOfSpecialisedImmuneCellAgents; ++i )
    {
        initialiseSpecialisedImmuneCellAgent(i, false);   
        currSpecialisedImmuneCellAgentId++; 
    }
}



/**********************
*   VirusCellModel::initialiseEpithelialCellAgent - Creates an epithelial cell agent and places it on the grid.
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
    // The time since last division will remain 0, if there is an existing agent object, which means we are handling a newly divided cell. That will be a completely fresh cell with 0 time since last division.
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

    // Allocate an arbitrary release delay
    repast::NormalGenerator releaseDelayGen = repast::Random::instance()->createNormalGenerator(epithCellVirionReleaseDelayAvg, epithCellVirionReleaseDelayStdev);
    double releaseDelay = releaseDelayGen.next();
    while(releaseDelay < 1 && releaseDelay < displayVirProtDelay)
    {
        releaseDelay = releaseDelayGen.next();
    }

    // Create the new Epithelial cell agent
    if( isExistingAgentObject == false )
    {
        repast::AgentId newAgentId(epithelialCellIndex, rank, 0);
        newAgentId.currentRank(rank);
        EpithelialCellAgent* newEpithelialCell = new EpithelialCellAgent(newAgentId, cellLifespan, cellAge, infectedCellLifespan, divisionRate, timeSinceLastDivision, releaseDelay, displayVirProtDelay, extracellularVirusReleaseProb, cellToCellTransmissionProb);
        context.addAgent(newEpithelialCell);

        // Set the new cell location
        repast::Point<int> agentLocation(xCoor + discreteGridSpace->dimensions().origin().getX() , yCoor + discreteGridSpace->dimensions().origin().getY() );
        discreteGridSpace->moveTo(newAgentId, agentLocation);
    }
    // If the epithelial cell agent object already exists, then we are handling division. In such case, that means, that we just need to reinitialise the
    // existing agent object with the newly generated parameters. Then the existing object will be representing a new cell.
    else if( theExistingCellObject != nullptr )
    {
        theExistingCellObject->set(theExistingCellObject->getId().currentRank(), cellLifespan, cellAge, theExistingCellObject->healthy, theExistingCellObject->seeminglyHealthy, 
        infectedCellLifespan, infectedTime, divisionRate, timeSinceLastDivision, releaseDelay, displayVirProtDelay, 
        theExistingCellObject->noModification, repast::AgentId(-1, -1, -1, -1), 
        extracellularVirusReleaseProb, cellToCellTransmissionProb);
    }

}



/**********************
*   VirusCellModel::initialiseVirionCellAgent - Creates a virion (virus particle) agent, sets all its parameters and places it on the grid.
**********************/
void VirusCellModel::initialiseVirionCellAgent(int virionIndex, bool isAReleasedVirus, int xCoor, int yCoor)
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

    // Assign an arbitrary age to the virion if it is from the starting population. If it is a released virus, then the age should be 0.
    int virionAge = 0;
    if( !isAReleasedVirus )
    {
        repast::IntUniformGenerator ageGenerator = repast::Random::instance()->createUniIntGenerator(0, virionLifespan);
        virionAge = ageGenerator.next();
    }

    // Assign the target cell penetration probability to the virion
    double virionPenetrationProb = virionCellPenetrationProbability;

    // Assign the virion clearance probability (that is probability that a free virus particle is cleared off through unmodelled immune mechanisms)
    double clearanceProb = virionClearanceProbability;

    // The coefficient of increasing the clearance probability, for each immune cell which is at the current location.
    // The more immune cells are present at that locaiton, the higher the probability that there are also 
    // more unmodelled immune cell mechanisms carried out there, hence higher clearance probability.
    double clearanceProbScaler = virionClearanceProbabilityScaler;

    VirionAgent* newVirion = new VirionAgent(newVirionId, virionLifespan, virionAge, virionPenetrationProb, clearanceProb, clearanceProbScaler);
    context.addAgent(newVirion);

    // Place the agent in the grid spatial projection. If it is a released virus then use the provided coordinates, to place it in the grid.
    // That is since, the released virions are initialised at the position of the cell which has released the,
    if( isAReleasedVirus )
    {
        repast::Point<int> virionLocation(xCoor, yCoor);
        discreteGridSpace->moveTo(newVirionId, virionLocation);
    }
    else
    {
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
void VirusCellModel::initialiseInnateImmuneCellAgent( int immuneCellId, bool isFreshCell )
{
    int rank = repast::RepastProcess::instance()->rank();

    repast::AgentId newInnateImmuneCellId(immuneCellId, rank, 2);
    newInnateImmuneCellId.currentRank(rank);

    repast::NormalGenerator lifespanGenerator = repast::Random::instance()->createNormalGenerator(innateImmuneCellAvgLifespan, innateImmuneCellLifespanStdev);
    int innateImmuneCellLifespan = lifespanGenerator.next();
    while( innateImmuneCellLifespan < 1 )
    {
        innateImmuneCellLifespan = lifespanGenerator.next();
    }
  

    int cellAge = 0;
    if( !isFreshCell )
    {
        repast::IntUniformGenerator ageGenerator = repast::Random::instance()->createUniIntGenerator(0, innateImmuneCellLifespan);
        cellAge = ageGenerator.next();
    }

    double infectedCellRecognitionProb = innateImmuneCellInfectedCellRecognitionProb;
    double infectedCellEliminationProb = innateImmuneCellInfectedCellEliminationProb;

    double specialisedImmuneCellRecruitProb = innateImmuneCellRecruitSpecImmuneCellProb;

    InnateImmuneCellAgent* newInnateImmuneCell = new InnateImmuneCellAgent(newInnateImmuneCellId, innateImmuneCellLifespan, cellAge ,infectedCellRecognitionProb, infectedCellEliminationProb, specialisedImmuneCellRecruitProb);
    context.addAgent(newInnateImmuneCell);

    // Place the agent in the grid spatial projection.
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
void VirusCellModel::initialiseSpecialisedImmuneCellAgent( int immuneCellId, bool isFreshCell )
{
    int rank = repast::RepastProcess::instance()->rank();

    repast::AgentId newSpecialisedImmuneCellId(immuneCellId, rank, 3);
    newSpecialisedImmuneCellId.currentRank(rank);

    repast::NormalGenerator lifespanGenerator = repast::Random::instance()->createNormalGenerator(specialisedImmuneCellAvgLifespan, specialisedImmuneCellLifespanStdev);
    int specialisedImmuneCellLifespan = lifespanGenerator.next();
    while( specialisedImmuneCellLifespan < 1 )
    {
        specialisedImmuneCellLifespan = lifespanGenerator.next();
    }

    int cellAge = 0;
    if( !isFreshCell )
    {
        repast::IntUniformGenerator ageGenerator = repast::Random::instance()->createUniIntGenerator(0, specialisedImmuneCellLifespan);
        ageGenerator.next();
    }

    double infectedCellRecognitionProb = specialisedImmuneCellInfectedCellRecognitionProb;
    double infectedCellEliminationProb = specialisedImmuneCellInfectedCellEliminationProb;

    SpecialisedImmuneCellAgent* newSpecialisedImmuneCell = new SpecialisedImmuneCellAgent(newSpecialisedImmuneCellId, specialisedImmuneCellLifespan, cellAge, infectedCellRecognitionProb, infectedCellEliminationProb);
    context.addAgent(newSpecialisedImmuneCell);

    // Place the agent in the grid spatial projection.
    repast::IntUniformGenerator gridXCoorGenerator = repast::Random::instance()->createUniIntGenerator(0, discreteGridSpace->dimensions().extents().getX() - 1);
    int specialisedImmuneCellX = discreteGridSpace->dimensions().origin().getX() + gridXCoorGenerator.next();

    repast::IntUniformGenerator gridYCoorGenerator = repast::Random::instance()->createUniIntGenerator(0, discreteGridSpace->dimensions().extents().getX() - 1);
    int specialisedImmuneCellY = discreteGridSpace->dimensions().origin().getY() + gridYCoorGenerator.next();

    repast::Point<int> immuneCellLoc(specialisedImmuneCellX, specialisedImmuneCellY);
    discreteGridSpace->moveTo(newSpecialisedImmuneCellId, immuneCellLoc);
}



/**********************
*   VirusCellModel::requestAgents - Function for requesting agents across processes
**********************/
void VirusCellModel::requestAgents(){
    // Get the rank
	int rank = repast::RepastProcess::instance()->rank();
    if(rank == 0)
    {

        std::vector<VirusCellInteractionAgents*> theAgents;
        
        // Get the local agents and make them perform a step. This will include all agents: Epithelial cells, Virions, Specialised and Non-Specialised immune cells.
        // They will also be in random order which will provide the stochasticity we need.
        context.selectAgents(repast::SharedContext<VirusCellInteractionAgents>::NON_LOCAL, theAgents, 0, false, -1);

        std::cout<<"Non Local Agentss: "<<theAgents.size()<<std::endl;
        std::vector<VirusCellInteractionAgents*>::iterator iter;

        // for( iter = theAgents.begin(); iter != theAgents.end(); ++iter )
        // {
        //     std::cout<<"Non local agent: "<<(*iter)->getId()<<std::endl;
        // }
    }



	// int worldSize= repast::RepastProcess::instance()->worldSize();
    // repast::AgentRequest req(rank);

	// repast::AgentRequest theRequest(rank);
	// for(int i = 0; i < worldSize; i++){                              // For each process
	// 	if(i != rank){                                           // ... except this one
	// 		std::vector<VirusCellInteractionAgents*> agents;        
	// 		context.selectAgents(3, agents, 0, false);                 // Choose 5 local agents randomly

	// 		for(size_t j = 0; j < agents.size(); j++){
    //             EpithelialCellAgent* theEpithelialCell = dynamic_cast<EpithelialCellAgent*>(agents[j]);
	// 			repast::AgentId local = theEpithelialCell->getId();          // Transform each local agent's id into a matching non-local one
	// 			repast::AgentId other(local.id(), i, 0);
	// 			other.currentRank(i);
	// 			req.addRequest(other);                      // Add it to the agent request
	// 		}
	// 	}
	// }
	// repast::RepastProcess::instance()->requestAgents<VirusCellInteractionAgents, VirusCellInteractionAgentPackage, VirusCellInteractionAgentsPackageProvider, VirusCellInteractionAgentsPackageReceiver>(context, req, *agentProvider, *agentReceiver, *agentReceiver);
}



/**********************
*   VirusCellModel::doSomething - Function which will do something on each step.
**********************/
void VirusCellModel::doSomething()
{
    std::vector<VirusCellInteractionAgents*> theBufferZoneAgents;
    context.selectAgents(repast::SharedContext<VirusCellInteractionAgents>::NON_LOCAL, theBufferZoneAgents);

    std::vector<VirusCellInteractionAgents*>::iterator iter;

    for( iter = theBufferZoneAgents.begin(); iter != theBufferZoneAgents.end(); ++iter )
    {
        if( (*iter)->getId().agentType() == 0 )
        {
            checkForCellDivision(*iter);
        }
    }


    std::vector<VirusCellInteractionAgents*> theAgents;
    
    // Get the local agents and make them perform a step. This will include all agents: Epithelial cells, Virions, Specialised and Non-Specialised immune cells.
    // They will also be in random order which will provide the stochasticity we need.
    context.selectAgents(repast::SharedContext<VirusCellInteractionAgents>::LOCAL, theAgents);

    for( iter = theAgents.begin(); iter != theAgents.end(); ++iter )
    {
        (*iter)->doStep(&context, discreteGridSpace);

        if( (*iter)->getId().agentType() == 0 )
        {
            checkForCellDivision(*iter);
            checkForCellVirionRelease(*iter);
            checkForCellToCellInfection(*iter);
        }
        else if ((*iter)->getId().agentType() == 2 )
        {
            checkForInnateImmuneCellRecruitment(*iter);
            checkForSpecialisedImmuneCellRecruitement(*iter);
        }
        else if( (*iter)->getId().agentType() == 3 )
        {
            checkForSpecialisedImmuneCellRecruitement(*iter);
        }
        
        
        if( (*iter)->getId().agentType() != 0 )
        {
            removeLocalAgentIfDead(*iter);
        }
    } 

    discreteGridSpace->balance();

    repast::RepastProcess::instance()->synchronizeAgentStatus<VirusCellInteractionAgents, VirusCellInteractionAgentPackage, VirusCellInteractionAgentsPackageProvider, 
        VirusCellInteractionAgentsPackageReceiver>(context, *agentProvider, *agentReceiver, *agentReceiver);


    repast::RepastProcess::instance()->synchronizeProjectionInfo<VirusCellInteractionAgents, VirusCellInteractionAgentPackage, VirusCellInteractionAgentsPackageProvider, 
        VirusCellInteractionAgentsPackageReceiver>(context, *agentProvider, *agentReceiver, *agentReceiver);

    // Synchronise all agents which are non-local to this process (The copies of non-local agents which this process owns)
    repast::RepastProcess::instance()->synchronizeAgentStates<VirusCellInteractionAgentPackage, VirusCellInteractionAgentsPackageProvider, 
        VirusCellInteractionAgentsPackageReceiver>(*agentProvider, *agentReceiver);
}



/**********************
*   VirusCellModel::checkForCellDivision - Function which checks if an Epithelial Cell agent is ready to divide
*       If it is, then it check where it will divide, and sets the new epithelial cell with its parameters
**********************/
void VirusCellModel::checkForCellDivision(VirusCellInteractionAgents* theEpithelialCellAgent)
{
    EpithelialCellAgent* epithelialCell = static_cast<EpithelialCellAgent*>(theEpithelialCellAgent);
    if( epithelialCell == nullptr )
    {
        std::cout<<"Incorrect agent object was provided!"<<std::endl;
        return;
    }

    // If the epithelial cell wants to modify a neighbouring cell and the modification it wants to do is to divide into it, 
    // then there needs to be a division into the stored cell
    if( epithelialCell->getTypeOfModifToNeighbCell() == epithelialCell->toDivideInto )
    {
        int rank = repast::RepastProcess::instance()->rank();

        // Get the id of the neighbouring cell where the division will happen.
        repast::AgentId cellToDivideIntoId = epithelialCell->getNeighbouringCellToModify();
        if(cellToDivideIntoId.id() != -1 && cellToDivideIntoId.agentType() == 0 && rank == cellToDivideIntoId.currentRank())
        {
            VirusCellInteractionAgents* cellToBeRevivedBaseClass = context.getAgent(cellToDivideIntoId);  
            EpithelialCellAgent* cellToBeRevived = static_cast<EpithelialCellAgent*>(cellToBeRevivedBaseClass);

            // Ensure we revive only dead cells. Since cross-process agent synchronisation happens at the end of the step,
            // There could potentially be some inconsistencies (local cell trying to revive a cell and a buffer zone cell trying to revive the same cell)
            if(cellToBeRevived->getInternalState() == cellToBeRevived->dead)
            {
                initialiseEpithelialCellAgent(-1, -1, -1, true, cellToBeRevived);
            }
        }
    }
}



/**********************
*   VirusCellModel::checkForCellToCellInfection - Function which checks if an Epithelial Cell agent is to infect a neighbouring cell through Cell-to-Cell virus transmission release
*       If it is, then it infects the corresponding neighbouring cell.
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
    if( epithelialCell->getTypeOfModifToNeighbCell() == EpithelialCellAgent::NeighbouringCellModificationType::toInfect && epithelialCell->getInternalState() == EpithelialCellAgent::InternalState::infected )
    {
        int rank = repast::RepastProcess::instance()->rank();

        // Get the id of the neighbouring cell which is to be infected.
        repast::AgentId cellToInfectId = epithelialCell->getNeighbouringCellToModify();
        if(cellToInfectId.id() != -1 && cellToInfectId.agentType() == 0 && rank == cellToInfectId.currentRank())
        {
            VirusCellInteractionAgents* cellToInfectBaseClass = context.getAgent(cellToInfectId);  
            EpithelialCellAgent* cellToBeInfected = static_cast<EpithelialCellAgent*>(cellToInfectBaseClass);

            // Ensure we revive only dead cells. Since cross-process agent synchronisation happens at the end of the step,
            // There could potentially be some inconsistencies (local cell trying to revive a cell and a buffer zone cell trying to revive the same cell)
            if(cellToBeInfected->getInternalState() == EpithelialCellAgent::InternalState::healthy)
            {
                cellToBeInfected->infect();
            }
        }
    }
}



/**********************
*   VirusCellModel::checkForCellVirionRelease - Function which checks if an Epithelial Cell agent is ready to release a new virion in the extracellular space.
**********************/
void VirusCellModel::checkForCellVirionRelease(VirusCellInteractionAgents* theEpithelialCellAgent)
{
    EpithelialCellAgent* epithelialCell = static_cast<EpithelialCellAgent*>(theEpithelialCellAgent);

    if( epithelialCell->isCellToReleaseVirion() == true )
    {
        repast::IntUniformGenerator releaseCountGen = repast::Random::instance()->createUniIntGenerator(1, 1);
        int numVirionsToRelease = 1;

        for( int i = 0; i < numVirionsToRelease; ++i)
        {
            ++currVirionAgentId;

            // Place the agent in the grid spatial projection.
            std::vector<int> epithelialCellLoc;
            discreteGridSpace->getLocation(epithelialCell->getId(), epithelialCellLoc);

            initialiseVirionCellAgent(currVirionAgentId, true, epithelialCellLoc[0], epithelialCellLoc[1]);
        }
    }

}



/**********************
*   VirusCellModel::checkForSpecialisedImmuneCellRecruitement - Function which checks if an immune cell agent is recruiting a new specialised immune cell.
**********************/
void VirusCellModel::checkForSpecialisedImmuneCellRecruitement(VirusCellInteractionAgents* theRecruitingImmuneCell)
{
    repast::AgentId theRecruitingCellId = theRecruitingImmuneCell->getId();

    if(theRecruitingCellId.agentType() == 2)
    {
        InnateImmuneCellAgent* theRecruitingInnateImmCell = static_cast<InnateImmuneCellAgent*>(theRecruitingImmuneCell);
        if(theRecruitingInnateImmCell->isToRecruitNewSpecImmuneCell())
        {
            // std::cout<<"Recruited a new Spec imm cell!"<<std::endl;

            initialiseSpecialisedImmuneCellAgent(currSpecialisedImmuneCellAgentId++, true);
        }
    }
    else if (theRecruitingCellId.agentType() == 3)
    {
        SpecialisedImmuneCellAgent* theRecruitingSpecialisedImmCell = static_cast<SpecialisedImmuneCellAgent*>(theRecruitingImmuneCell);
        if(theRecruitingSpecialisedImmCell->isToRecruitNewSpecImmuneCell())
        {
            initialiseSpecialisedImmuneCellAgent(currSpecialisedImmuneCellAgentId++, true);
        }
    }
    
}


 
void VirusCellModel::checkForInnateImmuneCellRecruitment(VirusCellInteractionAgents* theRecruitingImmuneCell)
{
    InnateImmuneCellAgent* theRecruitingInnateImmCell = static_cast<InnateImmuneCellAgent*>(theRecruitingImmuneCell);
    if( theRecruitingImmuneCell == nullptr )
    {
        std::cout<<"Incorrect pointer to innate immune cell agent was passed to checkForInnateImmuneCellRecruitment()!"<<std::endl;
        return;
    }

    if(theRecruitingInnateImmCell->isToRecruitNewInnateImmunceCell())
    {
        // std::cout<<"Recruited a new Innate imm cell!"<<std::endl;

        initialiseInnateImmuneCellAgent(currInnateImmuneCellAgendId++, true);
    }
}



/**********************
*   VirusCellModel::removeLocalAgentIfDead - Function which checks if a Virion/Innate Immune Cell/Specailise Immune Cell agent 
*       has the dead state and removes it from the simulation
**********************/
void VirusCellModel::removeLocalAgentIfDead(VirusCellInteractionAgents* theAgent)
{
    // We will ignore the Epithelial cell agents, since they are not removed from the simulation if they die.
    // That is since, epithelial cells can divide, and the division of a cell will basically "revive" a dead cell.
    // This will remove the need to create a new epithelial cell agent in the simulation.

    repast::AgentId theAgentId = theAgent->getId();
    int theAgentType = theAgentId.agentType();

    bool removeAgent = false;


    // Check if the agent state is dead and if it is remove it from the process
    if( theAgentType == 1 )
    {
        VirionAgent* theVirion = static_cast<VirionAgent*>(theAgent);
        removeAgent = (theVirion->getVirionState() == theVirion->dead);
    }
    else if( theAgentType == 2 )
    {
        InnateImmuneCellAgent* theInnateImmuneCell = static_cast<InnateImmuneCellAgent*>(theAgent);
        removeAgent = (theInnateImmuneCell->getCellState() == theInnateImmuneCell->dead);

        if(removeAgent)
        {
            std::vector<VirusCellInteractionAgents*> theLocalInnateImmuneCells;
    
            // Get all local innate immune cell agents
            context.selectAgents(repast::SharedContext<VirusCellInteractionAgents>::LOCAL, theLocalInnateImmuneCells, 2, false);

            while( theLocalInnateImmuneCells.size() <= countOfInnateImmuneCellAgents )
            {
                // Get all local innate immune cell agents
                context.selectAgents(repast::SharedContext<VirusCellInteractionAgents>::LOCAL, theLocalInnateImmuneCells, 2, false);

                // Since we remove an innnate immune cell, but we need to keep their counts quite constant, we will create a new one at the dead one's place
                initialiseInnateImmuneCellAgent(currInnateImmuneCellAgendId++, true);
            }
        }
    }
    else if( theAgentType == 3 )
    {
        SpecialisedImmuneCellAgent* theSpecialisedImmuneCell = static_cast<SpecialisedImmuneCellAgent*>(theAgent);
        removeAgent = (theSpecialisedImmuneCell->getCellState() == theSpecialisedImmuneCell->dead);
    }

    if( removeAgent )
    {
        repast::RepastProcess::instance()->agentRemoved( theAgentId );
        context.removeAgent( theAgentId );
    }

}



/**********************
*   VirusCellModel::initSchedule - Function which schedules all events
**********************/
void VirusCellModel::initSchedule(repast::ScheduleRunner& runner)
{
	runner.scheduleEvent(1, 1, repast::Schedule::FunctorPtr(new repast::MethodFunctor<VirusCellModel> (this, &VirusCellModel::doSomething)));
	runner.scheduleEndEvent(repast::Schedule::FunctorPtr(new repast::MethodFunctor<VirusCellModel> (this, &VirusCellModel::recordResults)));
	runner.scheduleStop(stopAt);

    // Data collection
	runner.scheduleEvent(0.1, 1, repast::Schedule::FunctorPtr(new repast::MethodFunctor<repast::DataSet>(agentsData, &repast::DataSet::record)));
	runner.scheduleEvent(0.2, 3, repast::Schedule::FunctorPtr(new repast::MethodFunctor<repast::DataSet>(agentsData, &repast::DataSet::write)));
	runner.scheduleEndEvent(repast::Schedule::FunctorPtr(new repast::MethodFunctor<repast::DataSet>(agentsData, &repast::DataSet::write)));

    runner.scheduleEvent(1.3, 10, repast::Schedule::FunctorPtr(new repast::MethodFunctor<VirusCellModel> (this, &VirusCellModel::requestAgents)));
}

 

/**********************
*   VirusCellModel::recordResults - Function which schedules all events
**********************/
void VirusCellModel::recordResults()
{
    props->putProperty("RunNumber","1");
    std::vector<std::string> keyOrder;
    keyOrder.push_back("RunNumber");
    keyOrder.push_back("stop.at");
    props->writeToSVFile("./output/results.csv", keyOrder);
}





/******************************************
*   Agent Package Provider class
******************************************/

/**********************
*   VirusCellInteractionAgentsPackageProvider::VirusCellInteractionAgentsPackageProvider - Constructor for the package provider
**********************/
VirusCellInteractionAgentsPackageProvider::VirusCellInteractionAgentsPackageProvider(repast::SharedContext<VirusCellInteractionAgents>* contextPtr): 
agentsContext(contextPtr)
{ 

}



/**********************
*   VirusCellInteractionAgentsPackageProvider::providePackage - Function for providing an agent package
**********************/
void VirusCellInteractionAgentsPackageProvider::providePackage(VirusCellInteractionAgents * agent, std::vector<VirusCellInteractionAgentPackage>& out){
    repast::AgentId id = agent->getId();
    int agentType = id.agentType();

    VirusCellInteractionAgentPackage package;
    if (agentType == 0)
    {
        EpithelialCellAgent* theEpithelialCell = static_cast<EpithelialCellAgent*>(agent);

        repast::AgentId agentToModify = theEpithelialCell->getNeighbouringCellToModify();
        VirusCellInteractionAgentPackage cellPackage(id.id(), id.startingRank(), id.agentType() ,id.currentRank(), theEpithelialCell->getLifespan(), 
            theEpithelialCell->getAge(), theEpithelialCell->getInternalState(), theEpithelialCell->getExternalState(), theEpithelialCell->getInfectedLifespan(),theEpithelialCell->getTimeInfected(), 
            theEpithelialCell->getDivisionRate(), theEpithelialCell->getTimeSinceLastDivision(), theEpithelialCell->getReleaseDelay(), theEpithelialCell->getDisplayVirProteinsDelay(), 
            theEpithelialCell->getExtracellularReleaseProb(), theEpithelialCell->getCellToCellTransmissionProb(),
            -1, -1, -1, -1, -1, -1,
            theEpithelialCell->getTypeOfModifToNeighbCell(), agentToModify.id(), agentToModify.startingRank(), agentToModify.agentType(), agentToModify.currentRank());
        
        package = cellPackage;
    }
    else if( agentType == 1)
    {
        VirionAgent* theVirion = static_cast<VirionAgent*>(agent);
        VirusCellInteractionAgentPackage virionPackage(id.id(), id.startingRank(), id.agentType() ,id.currentRank(), theVirion->getLifespan(), 
            theVirion->getAge(), theVirion->getVirionState(), 
            -1, -1, -1, -1, -1, -1, -1, -1, -1,
            theVirion->getPenetrationProb(), theVirion->getClearanceProb(), theVirion->getClearanceProbScaler(),
             -1, -1, -1,
             -1, -1, -1, -1, -1);
        
        package = virionPackage;
    }
    else if( agentType == 2 )
    {
        InnateImmuneCellAgent* theInnateImmuneCell = static_cast<InnateImmuneCellAgent*>(agent);
        VirusCellInteractionAgentPackage innImuneCellPackage(id.id(), id.startingRank(), id.agentType() ,id.currentRank(), theInnateImmuneCell->getLifespan(), 
            theInnateImmuneCell->getAge(), theInnateImmuneCell->getCellState(), 
            -1 ,-1, -1, -1, -1, -1, -1, -1, -1, 
            -1, -1, -1,
            theInnateImmuneCell->getInfCellRecognitionProb(), theInnateImmuneCell->getInfCellEliminationProb(), theInnateImmuneCell->getSpecImmuneCellRecruitProb(),
            -1, -1, -1, -1, -1);
        
        package = innImuneCellPackage;
    }
    else
    {
        SpecialisedImmuneCellAgent* theSpecialisedImmuneCell = static_cast<SpecialisedImmuneCellAgent*>(agent);
        VirusCellInteractionAgentPackage spcImmuneCellPackage(id.id(), id.startingRank(), id.agentType() ,id.currentRank(), theSpecialisedImmuneCell->getLifespan(), 
            theSpecialisedImmuneCell->getAge(), theSpecialisedImmuneCell->getCellState(),
            -1, -1, -1, -1, -1, -1 ,-1, -1, -1, 
            -1, -1, -1,
            theSpecialisedImmuneCell->getInfCellRecognitionProb(), theSpecialisedImmuneCell->getInfCellEliminationProb(),
            -1, 
            -1, -1, -1, -1, -1);
        
        package = spcImmuneCellPackage;
    }

    out.push_back(package);
}



/**********************
*   VirusCellInteractionAgentsPackageProvider::provideContent - Function for providing an agent package's content
**********************/
void VirusCellInteractionAgentsPackageProvider::provideContent(repast::AgentRequest req, std::vector<VirusCellInteractionAgentPackage>& out){
    std::vector<repast::AgentId> ids = req.requestedAgents();
    for(size_t i = 0; i < ids.size(); i++){
        providePackage(agentsContext->getAgent(ids[i]), out);
    }
}



/******************************************
*   Agent Package Receiver class
******************************************/

/**********************
*   VirusCellInteractionAgentsPackageReceiver::VirusCellInteractionAgentsPackageReceiver - Constructor for the package receiver
**********************/
VirusCellInteractionAgentsPackageReceiver::VirusCellInteractionAgentsPackageReceiver(repast::SharedContext<VirusCellInteractionAgents>* contextPtr): 
agentsContext(contextPtr)
{

}



/**********************
*   VirusCellInteractionAgentsPackageReceiver::createAgent - Function for creating an agent from a received agent package
**********************/
VirusCellInteractionAgents * VirusCellInteractionAgentsPackageReceiver::createAgent(VirusCellInteractionAgentPackage package){
    repast::AgentId theAgentId(package.id, package.rank, package.type, package.currentRank);

    if(package.type == 0)
    {
        return new EpithelialCellAgent( theAgentId, package.lifespan, package.age, EpithelialCellAgent::InternalState(package.privateState), 
            EpithelialCellAgent::ExternalState(package.publicState), package.infectedLifespan, package.timeInfected, package.divisionRate, package.timeSinceLastDivision, package.releaseDelay, package.displayVirProteinsDelay,
            EpithelialCellAgent::NeighbouringCellModificationType(package.typeOfChangeToMakeToAgent), repast::AgentId(package.agentToEditID, package.agentToEditStartRank, package.agentToEditType, package.agentToEditCurrRank),
            package.extracellularReleaseProb, package.cellToCellTransmissionProb );
    }  
    else if(package.type == 1)
    {
        return new VirionAgent( theAgentId, package.lifespan, package.age, VirionAgent::VirionStates(package.privateState), package.penetrationProbability, package.clearanceProbability, package.clearanceProbScaler);
    }
    else if(package.type == 2)
    {
        return new InnateImmuneCellAgent( theAgentId, package.lifespan, package.age, InnateImmuneCellAgent::InnateImmuneCellStates(package.privateState), package.infectedCellRecognitionProb, package.infectedCellEliminationProb, package.specialisedImmuneCellRecruitProb );
    }
    else
    {
        return new SpecialisedImmuneCellAgent( theAgentId, package.lifespan, package.age, SpecialisedImmuneCellAgent::SpecialisedImmuneCellStates(package.privateState), package.infectedCellRecognitionProb, package.infectedCellEliminationProb );
    }
}



/**********************
*   VirusCellInteractionAgentsPackageReceiver::updateAgent - Function for updating an agent with data from a received agent package
**********************/
void VirusCellInteractionAgentsPackageReceiver::updateAgent(VirusCellInteractionAgentPackage package){
    repast::AgentId theAgentId(package.id, package.rank, package.type);
    VirusCellInteractionAgents * theAgent = agentsContext->getAgent(theAgentId);

    if( package.type == 0 )
    {
        EpithelialCellAgent* theEpithelialCell = static_cast<EpithelialCellAgent*>(theAgent);
        theEpithelialCell->set(package.currentRank, package.lifespan, package.age, EpithelialCellAgent::InternalState(package.privateState), 
           EpithelialCellAgent::ExternalState(package.publicState), package.infectedLifespan, package.timeInfected, package.divisionRate, package.timeSinceLastDivision, package.releaseDelay, package.displayVirProteinsDelay,
           EpithelialCellAgent::NeighbouringCellModificationType(package.typeOfChangeToMakeToAgent), repast::AgentId(package.agentToEditID, package.agentToEditStartRank, package.agentToEditType, package.agentToEditCurrRank),
           package.extracellularReleaseProb, package.cellToCellTransmissionProb);
    }
    else if( package.type == 1 )
    {
        VirionAgent* theVirion = static_cast<VirionAgent*>(theAgent);
        theVirion->set(package.currentRank, package.lifespan, package.age, VirionAgent::VirionStates(package.privateState), package.penetrationProbability, package.clearanceProbability, package.clearanceProbScaler );
    }
    else if( package.type == 2 )
    {
        InnateImmuneCellAgent* theInnateImmuneCell = static_cast<InnateImmuneCellAgent*>(theAgent);
        theInnateImmuneCell->set(package.currentRank, package.lifespan, package.age, InnateImmuneCellAgent::InnateImmuneCellStates(package.privateState), package.infectedCellRecognitionProb, package.infectedCellEliminationProb, package.specialisedImmuneCellRecruitProb );
    }
    else
    {
        SpecialisedImmuneCellAgent* theSpecialisedImmuneCell = static_cast<SpecialisedImmuneCellAgent*>(theAgent);
        theSpecialisedImmuneCell->set(package.currentRank, package.lifespan, package.age, SpecialisedImmuneCellAgent::SpecialisedImmuneCellStates(package.privateState), package.infectedCellRecognitionProb, package.infectedCellEliminationProb );
    }
}





/***********************************
********** DATA COLLECTION *********
***********************************/

/******************************************
* Data Source class for tracking the count of epithelial cells in the model
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
        if(currentEpithCell->getInternalState() != currentEpithCell->dead)
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

    // Count only the alive epithelial cells
    for( iter = theEpithelialCells.begin(); iter != theEpithelialCells.end(); ++iter )
    {

        EpithelialCellAgent* currentEpithCell = static_cast<EpithelialCellAgent*>(*iter);
        if(currentEpithCell->getInternalState() == currentEpithCell->infected)
        {
            ++infectedEpithCells;
        }

    } 
    return infectedEpithCells;
}