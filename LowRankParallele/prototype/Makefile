include Makefile.opts

SRC_DIR   = ../tools
BUILD_DIR = build
DATA_DIR  = data
EXEC=prototype
SOURCES = prototype.c pspblasinfo.c psmatgen.c stringTools.c pblasIOtools.c nameFile.c
OBJECTS   = $(SOURCES:%.c=$(BUILD_DIR)/%.o)
LIB       = $(LIBSCALAPACK) $(LIBLAPACK) $(LIBBLAS)
INC       = $(INCBLAS) $(INCTOOLS)

COMPILE = $(CC) $(CFLAGS) $(INC) -c
LINK = $(F77) $(LDFLAGS)

all: $(OBJECTS) prototype

%.o: %.c
	$(COMPILE) $*.c -o $@

$(BUILD_DIR)/%.o: %.c
	$(COMPILE) $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(COMPILE) $< -o $@

prototype: $(OBJECTS) 
	$(LINK) $(OBJECTS) $(LIB) -lm -o $@

clean:
	rm -rf $(EXEC) $(OBJECTS) $(DATA_DIR)/*.out $(DATA_DIR)/*SOL.dat
