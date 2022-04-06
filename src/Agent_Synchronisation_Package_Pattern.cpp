// Agent_Sychronisation_Pattern.cpp

// Contains everything required for syncrhonising/moving agents across processes. 
// That includes the agent packages implementation, the package provider and receiver classes implementation
#include "Agent_Synchronisation_Package_Pattern.h"

// Include agent related files
#include "Epithelial_Cell_Agent.h"
#include "Virion_Agent.h"
#include "Innate_Immune_Cell.h"
#include "Specialised_Immune_Cell.h"


/*****************************
* Serializable Agent Package Data Struct
*****************************/
VirusCellInteractionAgentPackage::VirusCellInteractionAgentPackage(){ }



/**********************
* Constructor of the serializable package to be used for Epithelial Cell Agents
**********************/
VirusCellInteractionAgentPackage::VirusCellInteractionAgentPackage(int _id, int _rank, int _type, int _currentRank, double _lifespan, int _age, int _internalState, int _externalState, double _infectedLifespan, 
            int _timeInfected, double _divisionRate, int _timeSinceLastDivision, double _releaseDelay, double _displayVirProteinsDelay, 
            double _extracellularReleaseProb, double _cellToCellTransmissionProb,
            double _virionReleaseRate, int _countOfVirionsToRelease, double _virionReleaseRemainder,
            int _typeOfChangeToMakeToAgent, int _agentToEditID, int _agentToEditStartRank, int _agentToEditType, int _agentToEditCurrRank):
id(_id), 
rank(_rank), 
type(_type),
currentRank(_currentRank), 
lifespan(_lifespan), 
age(_age),
internalState(_internalState),

// Specific to epithelial cell agents
externalState(_externalState),
infectedLifespan(_infectedLifespan),
timeInfected(_timeInfected),
divisionRate(_divisionRate),
timeSinceLastDivision(_timeSinceLastDivision),
releaseDelay(_releaseDelay),
displayVirProteinsDelay(_displayVirProteinsDelay),
extracellularReleaseProb(_extracellularReleaseProb),
cellToCellTransmissionProb(_cellToCellTransmissionProb),
virionReleaseRate(_virionReleaseRate),
countOfVirionsToRelease(_countOfVirionsToRelease),
virionReleaseRemainder(_virionReleaseRemainder),

// Specific for handling buffer zone agent changes propagating back to original agents
typeOfChangeToMakeToAgent(_typeOfChangeToMakeToAgent),
agentToEditID(_agentToEditID),
agentToEditStartRank(_agentToEditStartRank),
agentToEditType(_agentToEditType),
agentToEditCurrRank(_agentToEditCurrRank)
{

}



/**********************
* Constructor of the serializable package to be used for Virion Agents
**********************/
VirusCellInteractionAgentPackage::VirusCellInteractionAgentPackage(int _id, int _rank, int _type, int _currentRank, double _lifespan, int _age, int _internalState,
            double _penetrationProbability, double _clearanceProbability, double _clearanceProbScaler):
id(_id), 
rank(_rank), 
type(_type),
currentRank(_currentRank), 
lifespan(_lifespan), 
age(_age),
internalState(_internalState),

// Specific to virion agents
penetrationProbability(_penetrationProbability),
clearanceProbability(_clearanceProbability),
clearanceProbScaler(_clearanceProbScaler)
{

}


/**********************
* Constructor of the serializable package to be used for Innate Immune Cell Agents
**********************/
VirusCellInteractionAgentPackage::VirusCellInteractionAgentPackage(int _id, int _rank, int _type, int _currentRank, double _lifespan, int _age, int _internalState, 
                                    double _infectedCellRecognitionProb, double _infectedCellEliminationProb,
                                    double _specialisedImmuneCellRecruitProb, double _specialisedImmuneCellRecruitRate,
                                    int _countOfSpecCellsToRecruit, double _specCellsRecruitRemainder,
                                    double _innateImmuneCellRecruitRate, double _countOfInnateCellsToRecruit,
                                    double _innateCellsRecruitRemainder):
id(_id), 
rank(_rank), 
type(_type),
currentRank(_currentRank), 
lifespan(_lifespan), 
age(_age),
internalState(_internalState),

// Specific to the innate immune cell agents
infectedCellRecognitionProb(_infectedCellRecognitionProb),
infectedCellEliminationProb(_infectedCellEliminationProb),
specialisedImmuneCellRecruitProb(_specialisedImmuneCellRecruitProb),
specialisedImmuneCellRecruitRate(_specialisedImmuneCellRecruitRate),
countOfSpecCellsToRecruit(_countOfSpecCellsToRecruit),
specCellsRecruitRemainder(_specCellsRecruitRemainder),
innateImmuneCellRecruitRate(_innateImmuneCellRecruitRate),
countOfInnateCellsToRecruit(_countOfInnateCellsToRecruit),
innateCellsRecruitRemainder(_innateCellsRecruitRemainder)
{

}



/**********************
* Constructor of the serializable package to be used for Specialised Immune Cell Agents
**********************/
VirusCellInteractionAgentPackage::VirusCellInteractionAgentPackage(int _id, int _rank, int _type, int _currentRank, double _lifespan, int _age, int _internalState, 
                                    double _infectedCellRecognitionProb, double _infectedCellEliminationProb,
                                    double _specialisedImmuneCellRecruitRate, int _countOfSpecCellsToRecruit, double _specCellsRecruitRemainder):
id(_id), 
rank(_rank), 
type(_type),
currentRank(_currentRank), 
lifespan(_lifespan), 
age(_age),
internalState(_internalState),

// Specific to the innate specialised cell agents
infectedCellRecognitionProb(_infectedCellRecognitionProb),
infectedCellEliminationProb(_infectedCellEliminationProb),
specialisedImmuneCellRecruitRate(_specialisedImmuneCellRecruitRate),
countOfSpecCellsToRecruit(_countOfSpecCellsToRecruit),
specCellsRecruitRemainder(_specCellsRecruitRemainder)
{

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
*   VirusCellInteractionAgentsPackageProvider::providePackage - Function for providing an agent package. 
*   Builds the package with all needed variables and passes it to the output vector
**********************/
void VirusCellInteractionAgentsPackageProvider::providePackage(VirusCellInteractionAgents * agent, std::vector<VirusCellInteractionAgentPackage>& out){
    repast::AgentId id = agent->getId();
    int agentType = id.agentType();

    // Build the package by adding the set of needed variables depending on the type of agent which needs to be synchronised.
    VirusCellInteractionAgentPackage package;
    if (agentType == 0)
    {
        EpithelialCellAgent* theEpithelialCell = static_cast<EpithelialCellAgent*>(agent);

        repast::AgentId agentToModify = theEpithelialCell->getNeighbouringCellToModify();
        VirusCellInteractionAgentPackage cellPackage(id.id(), id.startingRank(), id.agentType() ,id.currentRank(), theEpithelialCell->getLifespan(), 
            theEpithelialCell->getAge(), theEpithelialCell->getInternalState(), theEpithelialCell->getExternalState(), theEpithelialCell->getInfectedLifespan(),theEpithelialCell->getTimeInfected(), 
            theEpithelialCell->getDivisionRate(), theEpithelialCell->getTimeSinceLastDivision(), theEpithelialCell->getReleaseDelay(), theEpithelialCell->getDisplayVirProteinsDelay(), 
            theEpithelialCell->getExtracellularReleaseProb(), theEpithelialCell->getCellToCellTransmissionProb(), 
            theEpithelialCell->getVirionReleaseRate(), theEpithelialCell->getVirionCountToRelease(), theEpithelialCell->getVirionReleaseRemainder(),
            theEpithelialCell->getTypeOfModifToNeighbCell(), agentToModify.id(), agentToModify.startingRank(), agentToModify.agentType(), agentToModify.currentRank());
        
        package = cellPackage;
    }
    else if( agentType == 1)
    {
        VirionAgent* theVirion = static_cast<VirionAgent*>(agent);
        VirusCellInteractionAgentPackage virionPackage(id.id(), id.startingRank(), id.agentType() ,id.currentRank(), theVirion->getLifespan(), 
            theVirion->getAge(), theVirion->getVirionState(), 
            theVirion->getPenetrationProb(), theVirion->getClearanceProb(), theVirion->getClearanceProbScaler());
        
        package = virionPackage;
    }
    else if( agentType == 2 )
    {
        InnateImmuneCellAgent* theInnateImmuneCell = static_cast<InnateImmuneCellAgent*>(agent);
        VirusCellInteractionAgentPackage innImuneCellPackage(id.id(), id.startingRank(), id.agentType() ,id.currentRank(), theInnateImmuneCell->getLifespan(), 
            theInnateImmuneCell->getAge(), theInnateImmuneCell->getCellState(), 
            theInnateImmuneCell->getInfCellRecognitionProb(), theInnateImmuneCell->getInfCellEliminationProb(), theInnateImmuneCell->getSpecImmuneCellRecruitProb(),
            theInnateImmuneCell->getSpecialisedCellRecruitRate(), theInnateImmuneCell->getCountOfSpecialisedCellsToRecruit(), theInnateImmuneCell->getSpecialisedCellsRecruitRemainder(),
            theInnateImmuneCell->getInnateCellRecruitRate(), theInnateImmuneCell->getCountOfInnateCellsToRecruit(), theInnateImmuneCell->getInnateCellsRecruitRemainder());
        
        package = innImuneCellPackage;
    }
    else
    {
        SpecialisedImmuneCellAgent* theSpecialisedImmuneCell = static_cast<SpecialisedImmuneCellAgent*>(agent);
        VirusCellInteractionAgentPackage spcImmuneCellPackage(id.id(), id.startingRank(), id.agentType() ,id.currentRank(), theSpecialisedImmuneCell->getLifespan(), 
            theSpecialisedImmuneCell->getAge(), theSpecialisedImmuneCell->getCellState(),
            theSpecialisedImmuneCell->getInfCellRecognitionProb(), theSpecialisedImmuneCell->getInfCellEliminationProb(),
            theSpecialisedImmuneCell->getSpecialisedImmuneCellRecruitRateOfSpecCell(), theSpecialisedImmuneCell->getCountOfSpecCellsToRecruit(),
            theSpecialisedImmuneCell->getSpecCellsRecruitRemainder());
        
        package = spcImmuneCellPackage;
    }

    // Provide the package
    out.push_back(package);
}



/**********************
*   VirusCellInteractionAgentsPackageProvider::provideContent - Function for providing an all requested agents' packages
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

    // Create the correct agent type, using the appropriate package variable values.
    if(package.type == 0)
    {
        return new EpithelialCellAgent( theAgentId, package.lifespan, package.age, EpithelialCellAgent::InternalState(package.internalState), 
            EpithelialCellAgent::ExternalState(package.externalState), package.infectedLifespan, package.timeInfected, package.divisionRate, package.timeSinceLastDivision, package.releaseDelay, package.displayVirProteinsDelay,
            EpithelialCellAgent::NeighbouringCellModificationType(package.typeOfChangeToMakeToAgent), repast::AgentId(package.agentToEditID, package.agentToEditStartRank, package.agentToEditType, package.agentToEditCurrRank),
            package.extracellularReleaseProb, package.cellToCellTransmissionProb,
            package.virionReleaseRate, package.countOfVirionsToRelease, package.virionReleaseRemainder );
    }  
    else if(package.type == 1)
    {
        return new VirionAgent( theAgentId, package.lifespan, package.age, VirionAgent::VirionStates(package.internalState), package.penetrationProbability, package.clearanceProbability, package.clearanceProbScaler);
    }
    else if(package.type == 2)
    {
        return new InnateImmuneCellAgent( theAgentId, package.lifespan, package.age, InnateImmuneCellAgent::InnateImmuneCellStates(package.internalState), package.infectedCellRecognitionProb, package.infectedCellEliminationProb, package.specialisedImmuneCellRecruitProb,
                                        package.innateImmuneCellRecruitRate, package.specialisedImmuneCellRecruitRate, package.countOfInnateCellsToRecruit, 
                                        package.innateCellsRecruitRemainder, package.countOfSpecCellsToRecruit, package.specCellsRecruitRemainder );
    }
    else
    {
        return new SpecialisedImmuneCellAgent(theAgentId, package.lifespan, package.age, SpecialisedImmuneCellAgent::SpecialisedImmuneCellStates(package.internalState), 
                                            package.infectedCellRecognitionProb, package.infectedCellEliminationProb, package.specialisedImmuneCellRecruitRate, 
                                            package.countOfSpecCellsToRecruit, package.specCellsRecruitRemainder );
    }
}



/**********************
*   VirusCellInteractionAgentsPackageReceiver::updateAgent - Function for updating an agent with data from a received agent package
**********************/
void VirusCellInteractionAgentsPackageReceiver::updateAgent(VirusCellInteractionAgentPackage package){
    repast::AgentId theAgentId(package.id, package.rank, package.type);
    VirusCellInteractionAgents * theAgent = agentsContext->getAgent(theAgentId);

    // Update the correct agent type, using the appropriate package variable values.
    if( package.type == 0 )
    {
        EpithelialCellAgent* theEpithelialCell = static_cast<EpithelialCellAgent*>(theAgent);
        theEpithelialCell->set(package.currentRank, package.lifespan, package.age, EpithelialCellAgent::InternalState(package.internalState), 
           EpithelialCellAgent::ExternalState(package.externalState), package.infectedLifespan, package.timeInfected, package.divisionRate, package.timeSinceLastDivision, package.releaseDelay, package.displayVirProteinsDelay,
           EpithelialCellAgent::NeighbouringCellModificationType(package.typeOfChangeToMakeToAgent), repast::AgentId(package.agentToEditID, package.agentToEditStartRank, package.agentToEditType, package.agentToEditCurrRank),
           package.extracellularReleaseProb, package.cellToCellTransmissionProb,
           package.virionReleaseRate, package.countOfVirionsToRelease, package.virionReleaseRemainder );
    }
    else if( package.type == 1 )
    {
        VirionAgent* theVirion = static_cast<VirionAgent*>(theAgent);
        theVirion->set(package.currentRank, package.lifespan, package.age, VirionAgent::VirionStates(package.internalState), package.penetrationProbability, package.clearanceProbability, package.clearanceProbScaler );
    }
    else if( package.type == 2 )
    {
        InnateImmuneCellAgent* theInnateImmuneCell = static_cast<InnateImmuneCellAgent*>(theAgent);
        theInnateImmuneCell->set(package.currentRank, package.lifespan, package.age, InnateImmuneCellAgent::InnateImmuneCellStates(package.internalState), package.infectedCellRecognitionProb, package.infectedCellEliminationProb, 
        package.specialisedImmuneCellRecruitProb, package.innateImmuneCellRecruitRate, package.specialisedImmuneCellRecruitRate,
        package.countOfInnateCellsToRecruit, package.innateCellsRecruitRemainder, package.countOfSpecCellsToRecruit, package.specCellsRecruitRemainder );
    }
    else
    {
        SpecialisedImmuneCellAgent* theSpecialisedImmuneCell = static_cast<SpecialisedImmuneCellAgent*>(theAgent);
        theSpecialisedImmuneCell->set(package.currentRank, package.lifespan, package.age, SpecialisedImmuneCellAgent::SpecialisedImmuneCellStates(package.internalState), package.infectedCellRecognitionProb, 
                                    package.infectedCellEliminationProb, package.specialisedImmuneCellRecruitRate, package.countOfSpecCellsToRecruit, package.specCellsRecruitRemainder );
    }
}