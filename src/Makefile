.SUFFIXES: .o .i .f90 .f
.f90.o:
	$(CFT) -c $(FFLAGS) $*.f90 -o $*.o
.f.o:
	$(CFT) -c $(FFLAGS) $*.f -o $*.o

MAIN = BEF
#Source file
SRCS = main.f90

OBJS =$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCS)))

#Executable file
$(MAIN): $(OBJS)
	$(CFT) $(FFLAGS) $(LFLAGS) -o $(MAIN) $(OBJS)
