


# Variables you can change
# Are we in debug mode? (1= IC cylinder ; 2= Messy IC, relaxes to cylinder )
DEBUG = 0
#
# Append to flags
CPPFLAGS = -DDEBUG=$(DEBUG) 

# Code compiling related stuff
GCC = g++ 
CPPFLAGS += -O2 #-Wall
INCLUDES = -I ./eigen/

CPPFLAGS += $(INCLUDES)

OUTPUTNAME = FilamentCode

all: Filament_Main

Filament_Main : DCores_Main.o INIReader.o ini.o ErrorMessages.o TheState.o SetIC.o MagCylinder.o UpdateQ.o UpdateDQDPHI.o Integrals.o FindSS.o SolvePoisson.o SolveAmpere.o SolvePertPoisson.o SolvePoissonSOR.o SolveAmpereSOR.o PrintState.o EvalV.o calcDMDPHI.o 
	$(GCC) $(CPPFLAGS) -o $(OUTPUTNAME) $^
	rm ./*.o

clean:
	rm $(OUTPUTNAME) *.o
