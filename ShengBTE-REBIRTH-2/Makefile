include arch.make

objects:=conductivity.o config.o data.o input.o iterations.o misc.o gruneisen.o phonon_routines.o integrals.o processes.o scaling.o symmetry.o wedgetc.o dos_routines.o ShengBTE.o

all: ../ShengBTE
../ShengBTE: $(objects)
	$(MPIFC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)
conductivity.o conductivity.mod : conductivity.f90 data.mod config.mod 
config.o config.mod : config.f90 symmetry.mod data.mod misc.mod 
data.o data.mod : data.f90 
dos_routines.o dos_routines.mod : dos_routines.f90 data.mod config.mod 
gruneisen.o gruneisen.mod : gruneisen.f90 input.mod config.mod data.mod misc.mod 
input.o input.mod : input.f90 data.mod config.mod 
integrals.o integrals.mod : integrals.f90 data.mod config.mod 
iterations.o iterations.mod : iterations.f90 wedgetc.mod config.mod data.mod 
misc.o misc.mod : misc.f90 data.mod 
phonon_routines.o phonon_routines.mod : phonon_routines.f90 input.mod config.mod data.mod misc.mod 
processes.o processes.mod : processes.f90 config.mod data.mod misc.mod 
scaling.o scaling.mod : scaling.f90 data.mod config.mod 
symmetry.o symmetry.mod : symmetry.f90 data.mod 
wedgetc.o wedgetc.mod : wedgetc.f90 misc.mod config.mod symmetry.mod
ShengBTE.o: integrals.mod gruneisen.mod dos_routines.mod phonon_routines.mod scaling.mod conductivity.mod iterations.mod wedgetc.mod processes.mod input.mod config.mod misc.mod data.mod 
%.o %.mod: %.f90
	$(MPIFC) $(FFLAGS) -c -o $*.o $<
	touch $*.mod
clean:
	rm -f ../ShengBTE *.mod *.o
