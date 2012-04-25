all: KeccakTools

SOURCES = \
    Sources/duplex.cpp \
    Sources/genKATShortMsg.cpp \
    Sources/Keccak.cpp \
    Sources/KeccakCrunchyContest.cpp \
    Sources/Keccak-f.cpp \
    Sources/Keccak-f25LUT.cpp \
    Sources/Keccak-fAffineBases.cpp \
    Sources/Keccak-fCodeGen.cpp \
    Sources/Keccak-fDCEquations.cpp \
    Sources/Keccak-fDCLC.cpp \
    Sources/Keccak-fDisplay.cpp \
    Sources/Keccak-fEquations.cpp \
    Sources/Keccak-fParity.cpp \
    Sources/Keccak-fParityBounds.cpp \
    Sources/Keccak-fParts.cpp \
    Sources/Keccak-fPositions.cpp \
    Sources/Keccak-fPropagation.cpp \
    Sources/Keccak-fState.cpp \
    Sources/Keccak-fTrailCore3Rounds.cpp \
    Sources/Keccak-fTrailCoreInKernelAtC.cpp \
    Sources/Keccak-fTrailCoreParity.cpp \
    Sources/Keccak-fTrailCoreRows.cpp \
    Sources/Keccak-fTrailExtension.cpp \
    Sources/Keccak-fTrails.cpp \
    Sources/main.cpp \
    Sources/padding.cpp \
    Sources/progress.cpp \
    Sources/sponge.cpp \
    Sources/spongetree.cpp \
    Sources/transformations.cpp

HEADERS = \
    Sources/duplex.h \
    Sources/Keccak.h \
    Sources/KeccakCrunchyContest.h \
    Sources/Keccak-f.h \
    Sources/Keccak-f25LUT.h \
    Sources/Keccak-fAffineBases.h \
    Sources/Keccak-fCodeGen.h \
    Sources/Keccak-fDCEquations.h \
    Sources/Keccak-fDCLC.h \
    Sources/Keccak-fDisplay.h \
    Sources/Keccak-fEquations.h \
    Sources/Keccak-fParity.h \
    Sources/Keccak-fParityBounds.h \
    Sources/Keccak-fParts.h \
    Sources/Keccak-fPositions.h \
    Sources/Keccak-fPropagation.h \
    Sources/Keccak-fState.h \
    Sources/Keccak-fTrailCore3Rounds.h \
    Sources/Keccak-fTrailCoreInKernelAtC.h \
    Sources/Keccak-fTrailCoreParity.h \
    Sources/Keccak-fTrailCoreRows.h \
    Sources/Keccak-fTrailExtension.h \
    Sources/Keccak-fTrails.h \
    Sources/padding.h \
    Sources/progress.h \
    Sources/sponge.h \
    Sources/spongetree.h \
    Sources/types.h \
    Sources/transformations.h \
    Sources/translationsymmetry.h \
    Sources/types.h

BINDIR = bin

$(BINDIR):
	mkdir -p $(BINDIR)

OBJECTS = $(addprefix $(BINDIR)/, $(notdir $(patsubst %.cpp,%.o,$(SOURCES))))

CFLAGS = -O3 -g0

VPATH = Sources

INCLUDES = -ISources

$(BINDIR)/%.o:%.cpp $(HEADERS)
	$(CXX) $(INCLUDES) $(CFLAGS) -c $< -o $@

.PHONY: KeccakTools

KeccakTools: bin/KeccakTools

bin/KeccakTools:  $(BINDIR) $(OBJECTS)  $(HEADERS)
	$(CXX) $(CFLAGS) -o $@ $(OBJECTS)

clean:
	rm -rf bin/
