all: KeccakTools

SOURCES=$(wildcard Sources/*.cpp)

BINDIR = bin

$(BINDIR):
	mkdir -p $(BINDIR)

OBJECTS = $(addprefix $(BINDIR)/, $(notdir $(patsubst %.cpp,%.o,$(SOURCES))))

CFLAGS = -O3 -g0 -Wreorder

VPATH = Sources

INCLUDES = -ISources

-include $(addsuffix .d, $(OBJECTS))

$(BINDIR)/%.o:%.cpp
	$(CXX) $(INCLUDES) $(CFLAGS) -c $< -o $@
	@$(CXX) $(INCLUDES) -MM $(CFLAGS) $< > $@.d.tmp
	@sed -e 's|.*:|$@:|' < $@.d.tmp > $@.d
	@rm $@.d.tmp

.PHONY: KeccakTools

KeccakTools: bin/KeccakTools

bin/KeccakTools:  $(BINDIR) $(OBJECTS)
	$(CXX) $(CFLAGS) -o $@ $(OBJECTS)

clean:
	rm -rf bin/
