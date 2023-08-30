CC	= mpicxx
PETSC=$(PETSC_DIR)
PETSCLIB=$(PETSC)/lib
PETSCINC1=$(PETSC)/include
PETSCINC2=$(PETSC)/include
PETSCARCH=$(PETSC)/linux-prime
HYPRE=$(PETSC)/linux-prime


#PETSCLIB=/safl/software/aegean/petsc/3.2-p6-openmpi-1.6-gcc-4.7.0/lib
#PETSCINC1=/safl/software/aegean/petsc/3.2-p6-openmpi-1.6-gcc-4.7.0/include
#PETSCINC2=/safl/software/aegean/petsc/3.2-p6-openmpi-1.6-gcc-4.7.0/include

#TECINC=./#$(TEC360HOME)/include
#TECLIB=./#$(TEC360HOME)/lib

TECINC  = /safl/software/x86_64/tecplot/360_2009_R2/include
TECLIB  = /safl/software/x86_64/tecplot/360_2009_R2/lib


FFTWINC = /safl/software/x86_64/fftw/3.2.2/include
FFTWLIB = /safl/software/x86_64/fftw/3.2.2/lib



ACMLLIB=$(ACML)/lib

HYPRELIB = $(HYPRE)/lib
HYPREINC = $(HYPRE)/include

#########

LIBDIR=-L$(ACMLLIB) -L$(PETSCLIB) -L$(TECLIB) -L$(HYPRELIB) -L$(FFTWLIB)

LIBFLAG	=	-lpthread -lrt -ldl -lstdc++ \
		-lpetsc -lHYPRE -lgfortran   

SOURCEC	=	bcs.c bmv.c compgeom.c ibm.c ibm_io.c init.c \
		main.c metrics.c poisson.c rhs.c timeadvancing.c \
		timeadvancing1.c variables.c fsi.c implicitsolver.c\
		fsi_move.c solvers.c rhs2.c wallfunction.c \
		les.c k-omega.c distance.c level.c momentum.c poisson_hypre.c rotor_model.c wallmodel.c tmprt.c sediment_transport.c convection_diffusion.c 


OBJSC	=	bcs.o bmv.o compgeom.o ibm.o ibm_io.o init.o \
		main.o metrics.o poisson.o rhs.o timeadvancing.o \
		timeadvancing1.o variables.o fsi.o implicitsolver.o\
		fsi_move.o solvers.o rhs2.o wallfunction.o \
		les.o k-omega.o distance.o level.o momentum.o poisson_hypre.o rotor_model.o wallmodel.o tmprt.o sediment_transport.o convection_diffusion.o

CPPFLAGS =	-DNDEBUG -I$(PETSCINC1) -I$(PETSCINC2) -I$(TECINC) -I$(HYPREINC)  -I$(FFTWINC) -DTECIO=1 -O3 -no-pie 

ALLFLAGS =	-Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -g3  -L/mmfs1/home/trung.le/Petsc/petsc-3.2-p7/linux-prime/lib  -lpetsc -lX11 -lpthread -Wl,-rpath,/mmfs1/home/trung.le/Petsc/petsc-3.2-p7/linux-prime/lib -lHYPRE -Wl,-rpath,/cm/local/apps/gcc/10.2.0/lib/gcc/x86_64-linux-gnu/10.2.0 -Wl,-rpath,/cm/local/apps/gcc/10.2.0/lib64 -Wl,-rpath,/cm/local/apps/gcc/10.2.0/lib -Wl,-rpath,/mmfs1/apps/spack/0.16.1/linux-rhel8-zen2/gcc-10.2.0/hwloc-1.11.11-znktbn53x77rlqwk7jyl4aezf2hlaueh/lib -Wl,-rpath,/mmfs1/apps/spack/0.16.1/linux-rhel8-zen2/gcc-10.2.0/zlib-1.2.11-e4yw3y55ty527weaskhkdoabx5d7ptfc/lib -Wl,-rpath,/mmfs1/apps/spack/0.16.1/linux-rhel8-zen2/gcc-10.2.0/openmpi-3.1.6-sxqyko7cbnnwtdndem2gzrkgljwzrx6p/lib -lstdc++ -lflapack -lfblas -lm -L/mmfs1/apps/spack/0.16.1/linux-rhel8-zen2/gcc-10.2.0/hwloc-1.11.11-znktbn53x77rlqwk7jyl4aezf2hlaueh/lib -L/mmfs1/apps/spack/0.16.1/linux-rhel8-zen2/gcc-10.2.0/zlib-1.2.11-e4yw3y55ty527weaskhkdoabx5d7ptfc/lib -L/mmfs1/apps/spack/0.16.1/linux-rhel8-zen2/gcc-10.2.0/openmpi-3.1.6-sxqyko7cbnnwtdndem2gzrkgljwzrx6p/lib -L/cm/local/apps/gcc/10.2.0/lib/gcc/x86_64-linux-gnu/10.2.0 -L/cm/local/apps/gcc/10.2.0/lib64 -L/cm/local/apps/gcc/10.2.0/lib -ldl -lmpi -lgcc_s -lpthread -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lgfortran -lm -lgfortran -lm -lquadmath -lm -lstdc++ -ldl -lmpi -lgcc_s -lpthread -ldl 	
test:	$(OBJSC)
	$(CC) -o testt $(OBJSC) $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

run:	$(OBJSC)
	$(CC) -o testt $(OBJSC) $(ALLFLAGS)

data: data.o
	$(CC) -o data data.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) /safl/software/x86_64/tecplot/360_2009_R2/lib/tecio64.a

data05: data05.o
	$(CC) -o data05 data05.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) -ltecio64

itfcsearch: itfcsearch.o
	$(CC) -o itfcsearch itfcsearch.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

shear: shear.o

data1: data1.o
	$(CC) -o data1 data1.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) /safl/software/x86_64/tecplot/360_2009_R2/lib/tecio64.a

hill: hill.o
	$(CC) -o hill hill.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

xyz: xyz2Plot3d.o
	$(CC) -o xyz xyz2Plot3d.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

clean:
	rm *.o 
turbine: TurbineLoc.o
	$(CC) -o TurbineLoc TurbineLoc.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)
terrain: terrain.o
	$(CC) -o terrain terrain.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)
Tec2UCD: OSLTec2UCD.o
	$(CC) -o Tec2UCD OSLTec2UCD.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

GenerateInflow: GenerateInflow__.o
	$(CC) -o GenerateInflow GenerateInflow__.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) /safl/software/x86_64/tecplot/360_2009_R2/lib/tecio64.a

ReadKevin_3dhilldat: ReadKevin_3dhilldat.o
	$(CC) -o ReadKevin_3dhilldat ReadKevin_3dhilldat.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) /safl/software/x86_64/tecplot/360_2009_R2/lib/tecio64.a

UCD2Plot3d: ucd2Plot3d_OSL.o
	$(CC) -o UCD2Plot3d ucd2Plot3d_OSL.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) /safl/software/x86_64/tecplot/360_2009_R2/lib/tecio64.a

GenInflow2: GenInflow2.o
	$(CC) -o GenInflow2 GenInflow2.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) /safl/software/x86_64/tecplot/360_2009_R2/lib/tecio64.a


SPAnal: SPAnal.o
	$(CC) -o SPAnal SPAnal.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) /safl/software/x86_64/tecplot/360_2009_R2/lib/tecio64.a

GenInflow3: GenInflow3.o
	$(CC) -o GenInflow3 GenInflow3.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) /safl/software/x86_64/tecplot/360_2009_R2/lib/tecio64.a




