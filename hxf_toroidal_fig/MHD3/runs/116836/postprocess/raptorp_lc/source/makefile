CC = h5cc  -I.  -I./src
CFLAGS =  -fopenmp  -std=c99  -O2 -lm -lgsl -Wall -Wno-unused-but-set-variable
LDFLAGS = -fopenmp -lm -lgsl -lgslcblas 

VPATH=./src:
CPATH=./src:

TARGET ?= RAPTOR
RAPTOR_NUM_FREQUENCIES ?= 4
RAPTOR_IMAGE_OUTPUT ?= 1

OBJDIR=build/$(TARGET)
BINDIR=bin
BIN=$(BINDIR)/$(TARGET)

CFLAGS += -DRAPTOR_NUM_FREQUENCIES=$(RAPTOR_NUM_FREQUENCIES) -DRAPTOR_IMAGE_OUTPUT=$(RAPTOR_IMAGE_OUTPUT)

SOURCES=main.c core.c io.c GRmath.c gr_integrator.c rte_integrator.c pol_rte_integrator.c metric.c pol_emission.c tetrad.c model.c constants.c camera.c pluto.c
OBJECTS := $(patsubst %.c,$(OBJDIR)/%.o,$(SOURCES))

all: build

build: create_directories $(BIN)

$(BIN): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

create_directories:
	@test -d $(OBJDIR) || mkdir -p $(OBJDIR)
	@test -d $(BINDIR) || mkdir -p $(BINDIR)


$(OBJECTS): $(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS)  -c $^ -o $@


RAPTOR_LC4_IMAGE:
	$(MAKE) build TARGET=$@ RAPTOR_NUM_FREQUENCIES=4 RAPTOR_IMAGE_OUTPUT=1

RAPTOR_LC4_SPEC:
	$(MAKE) build TARGET=$@ RAPTOR_NUM_FREQUENCIES=4 RAPTOR_IMAGE_OUTPUT=0

RAPTOR_SED50_IMAGE:
	$(MAKE) build TARGET=$@ RAPTOR_NUM_FREQUENCIES=50 RAPTOR_IMAGE_OUTPUT=1

RAPTOR_SED50_SPEC:
	$(MAKE) build TARGET=$@ RAPTOR_NUM_FREQUENCIES=50 RAPTOR_IMAGE_OUTPUT=0

clean:
	rm -rf build bin
