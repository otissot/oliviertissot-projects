CC = g++
CCFLAGS = -Wall -g
LDFLAGS =
LIBRARY = -L/usr/X11R6/lib  -lglut -lGLU  -lGL
INCLUDES = -I/usr/X11R6/include

PartieI:
EXEC1 = PartieI
SRCS1 = partie_I.cpp mat_cr.cpp vect.cpp maillage.cpp assemblage.cpp
OBJS1 = partie_I.o mat_cr.o vect.o maillage.o assemblage.o

%.o:%.cpp
	$(CC) $(CCFLAGS) $(INCLUDES) -c $<

all: $(EXEC1)

$(EXEC1) : $(OBJS1)
	$(CC) $(CCFLAGS) $(OBJS1) $(LIBRARY) $(LDFLAGS) -o $@

PartieII:
EXEC2 = PartieII
SRCS2 = partie_II.cpp mat_cr.cpp vect.cpp maillage.cpp assemblage.cpp
OBJS2 = partie_II.o mat_cr.o vect.o maillage.o assemblage.o

%.o:%.cpp
	$(CC) $(CCFLAGS) $(INCLUDES) -c $<

all: $(EXEC2)

$(EXEC2) : $(OBJS2)
	$(CC) $(CCFLAGS) $(OBJS2) $(LIBRARY) $(LDFLAGS) -o $@

affichage:
EXEC3 = affichage
SRCS3 = affichage.cpp mat_cr.cpp vect.cpp maillage.cpp
OBJS3 = affichage.o mat_cr.o vect.o maillage.o

%.o:%.cpp
	$(CC) $(CCFLAGS) $(INCLUDES) -c $<

all: $(EXEC3)

$(EXEC3) : $(OBJS3)
	$(CC) $(CCFLAGS) $(OBJS3) $(LIBRARY) $(LDFLAGS) -o $@

clean:
	rm -f $(OBJS1) $(EXEC1) $(OBJS2) $(EXEC2) $(OBJS3) $(EXEC3)s
