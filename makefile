#************************************************************************************************************
#
# Repast HPC Tutorial Makefile
#
#************************************************************************************************************

include ./env

.PHONY: all
all : 


.PHONY: clean_output_files
clean_output_files:
	rm -f *.csv
	rm -f *.txt
	rm -f ./output/*.csv
	rm -f ./output/*.txt
	rm -f ./logs/*.*

.PHONY: clean_compiled_files
clean_compiled_files:
	rm -f *.exe
	rm -f ./bin/*.exe
	rm -f *.o
	rm -f ./object/*.o
	
.PHONY: clean
clean: clean_compiled_files clean_output_files remove_subdirectories
	rm -f *.cpp
	rm -f ./src/*.cpp
	rm -f *.props
	rm -f ./props/*.props
	rm -f ./include/*.h


.PHONY: Virus_Cell_Sim
Virus_Cell_Sim: clean_compiled_files
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) -I./include -c ./src/Virus_Cell_Main.cpp -o ./objects/Virus_Cell_Main.o
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) -I./include -c ./src/Virus_Cell_Model.cpp -o ./objects/Virus_Cell_Model.o
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) -I./include -c ./src/Virus_Cell_Agent.cpp -o ./objects/Virus_Cell_Agent.o
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) -I./include -c ./src/Epithelial_Cell_Agent.cpp -o ./objects/Epithelial_Cell_Agent.o
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) -I./include -c ./src/Virion_Agent.cpp -o ./objects/Virion_Agent.o
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) -I./include -c ./src/Innate_Immune_Cell.cpp -o ./objects/Innate_Immune_Cell.o
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) -I./include -c ./src/Specialised_Immune_Cell.cpp -o ./objects/Specialised_Immune_Cell.o
	$(MPICXX) $(BOOST_LIB_DIR) $(REPAST_HPC_LIB_DIR) -o ./bin/Virus_Cell_Model.exe  ./objects/Virus_Cell_Main.o ./objects/Virus_Cell_Model.o ./objects/Virus_Cell_Agent.o ./objects/Epithelial_Cell_Agent.o ./objects/Virion_Agent.o  ./objects/Innate_Immune_Cell.o ./objects/Specialised_Immune_Cell.o -O3 $(REPAST_HPC_LIB) $(BOOST_LIBS)