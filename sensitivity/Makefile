PROGNAME      = analysis 
SOURCES       = analysis.cc simulation.cc

INCLUDES      = 
OBJECTS       = $(patsubst %.cc, %.o, $(SOURCES))
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTLIBS     := $(shell root-config --nonew --libs)
CFLAGS       += $(ROOTCFLAGS)
LIBS         += $(ROOTLIBS)
#  Not sure why Minuit isn't being included -- put in by hand
#  #
LIBS         += -lMinuit
LDFLAGS       = -O

#  Add RooFit and RooStats
#  #
LIBS        += -lRooFit -lRooFitCore -lRooStats

$(PROGNAME):    $(OBJECTS)
			g++ -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

%.o : %.cc $(INCLUDES)
		g++ ${CFLAGS} -c  -g -o $@ $<

test:
		@echo $(ROOTCFLAGS)

clean:	
		-rm -f ${PROGNAME} ${OBJECTS}
