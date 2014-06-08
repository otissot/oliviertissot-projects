# Le compilateur :
CC = g++

# Les opions de compilations (-g : le debuggueur, -Wall : les warnings, -O3 : optimisation, etc) :
#CCFLAGS = -Wall -g
CCFLAGS = -O3

# Pour pouvoir utiliser OpenGL :
OPENGL_LIB = -L/usr/X11R6/lib  -lglut -lGLU  -lGL
OPENGL_INC = -I/usr/X11R6/include

# A changer selon sa configuration :
UMFPACK_INC_INSTALL_DIRECTORY = /usr/local/include
UMFPACK_LIB_INSTALL_DIRECTORY = /usr/local/lib

# Pour pouvoir utiliser UMFPACK (on a besoin de BLAS) :
UMFPACK_INC = -I$(UMFPACK_INC_INSTALL_DIRECTORY)
# Pour pouvoir utiliser BLAS :
BLAS_LIB =  -lblas
UMFPACK_LIB = -L$(UMFPACK_LIB_INSTALL_DIRECTORY) -lumfpack -lcholmod -lamd -lcolamd -lsuitesparseconfig -lrt $(BLAS_LIB)

EXEC1 = glplotiso
SRCS1 = glplotiso.cpp Maillage.cpp
OBJS1 = glplotiso.o Maillage.o

EXEC2 = glplotisoP2
SRCS2 = glplotisoP2.cpp Maillage.cpp
OBJS2 = glplotisoP2.o Maillage.o

EXEC3 = poisson2d
SRCS3 = poisson2d.cpp Maillage.cpp MapMatrix.cpp Param.cpp MatriceElementaire.cpp QuadratureFormular.cpp
OBJS3 = poisson2d.o Maillage.o MapMatrix.o Param.o MatriceElementaire.o QuadratureFormular.o

EXEC4 = stokes2d
SRCS4 = stokes2d.cpp Maillage.cpp MapMatrix.cpp Param.cpp MatriceElementaire.cpp QuadratureFormular.cpp
OBJS4 = stokes2d.o Maillage.o MapMatrix.o Param.o MatriceElementaire.o QuadratureFormular.o

EXEC5 = NavierStokes2d
SRCS5 = NavierStokes2d.cpp Maillage.cpp MapMatrix.cpp Param.cpp MatriceElementaire.cpp SecondMembreNS.cpp QuadratureFormular.cpp
OBJS5 = NavierStokes2d.o Maillage.o MapMatrix.o Param.o MatriceElementaire.o SecondMembreNS.o QuadratureFormular.o


poisson2d.o:poisson2d.cpp
	$(CC) $(CCFLAGS) $(UMFPACK_INC) -c $<
stokes2d.o:stokes2d.cpp
	$(CC) $(CCFLAGS) $(UMFPACK_INC) -c $<
NavierStokes2d.o:NavierStokes2d.cpp
	$(CC) $(CCFLAGS) $(UMFPACK_INC) -c $<
glplotiso.o:glplotiso.cpp
	$(CC) $(CCFLAGS) $(OPENGL_INC) -c $<
glplotisoP2.o:glplotisoP2.cpp
	$(CC) $(CCFLAGS) $(OPENGL_INC) -c $<
%.o:%.cpp
	$(CC) $(CCFLAGS) -c $<

plot: $(EXEC1) $(EXEC2)

plotP1: $(EXEC1)
plotP2: $(EXEC2)

poisson2D: $(EXEC3)

stokes2D: $(EXEC4)

NS2D: $(EXEC5)

$(EXEC1) : $(OBJS1)
	$(CC) $(CCFLAGS) $(OBJS1) $(OPENGL_LIB) -o $@

$(EXEC2) : $(OBJS2)
	$(CC) $(CCFLAGS) $(OBJS2) $(OPENGL_LIB) -o $@

$(EXEC3) : $(OBJS3)
	$(CC) $(CCFLAGS) $(OBJS3) $(UMFPACK_LIB) -o $@

$(EXEC4) : $(OBJS4)
	$(CC) $(CCFLAGS) $(OBJS4) $(UMFPACK_LIB) -o $@

$(EXEC5) : $(OBJS5)
	$(CC) $(CCFLAGS) $(OBJS5) $(UMFPACK_LIB) -o $@


clean:
	@rm -f *.*~ *~ $(OBJS1) $(EXEC1) $(OBJS2) $(EXEC2) $(OBJS3) $(EXEC3) $(OBJS4) $(EXEC4) $(OBJS5) $(EXEC5)
	@ echo "Repertoire propre !"
