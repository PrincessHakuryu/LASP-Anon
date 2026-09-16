export FFLAGS = -O2 -Mpreprocess -acc -gpu=cc80,cc89  -Minfo=accel,mp -Minline -traceback -mp
export LDFLAGS=-L/home14/liangqm/spglib-2.6.0/build -lsymspg
export MPIFC=mpif90
MKL=$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group	\
$(MKLROOT)/lib/intel64/libmkl_intel_lp64.a				\
 $(MKLROOT)/lib/intel64/libmkl_sequential.a				\
 $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
export LAPACK=$(MKL)
export LIBS=$(LAPACK)
