all: KeccakTools

SOURCES = \
    Sources/genKATShortMsg.cpp \
    Sources/Keccak.cpp \
    Sources/Keccak-f.cpp \
    Sources/Keccak-f25LUT.cpp \
    Sources/Keccak-fCodeGen.cpp \
    Sources/Keccak-fEquations.cpp \
    Sources/Keccak-fParts.cpp \
    Sources/main.cpp \
    Sources/sponge.cpp \
    Sources/transformations.cpp

HEADERS = \
    Sources/Keccak.h \
    Sources/Keccak-f.h \
    Sources/Keccak-f25LUT.h \
    Sources/Keccak-fCodeGen.h \
    Sources/Keccak-fEquations.h \
    Sources/Keccak-fParts.h \
    Sources/sponge.h \
    Sources/types.h \
    Sources/transformations.h

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
