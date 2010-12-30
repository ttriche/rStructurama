LIBS     =  -lm -lgomp
CFLAGS	 =  -DOMP -fopenmp
CC       =   g++-4.5

OBJECTS  =   structurama.o hungarian.o interface.o integer.o iomanager.o MbBitfield.o machineInfo.o MbRandom.o mcmc.o model.o nexusfile.o observation.o partition.o samples.o stirling.o

PROGS    = st2

all:		$(PROGS)

st2:		$(OBJECTS)
		$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o st2
		
structurama.o:	structurama.cpp
		$(CC) $(CFLAGS) -c structurama.cpp

hungarian.o:	hungarian.cpp
		$(CC) $(CFLAGS) -c hungarian.cpp

interface.o:	interface.cpp
		$(CC) $(CFLAGS) -c interface.cpp

integer.o:	integer.cpp
		$(CC) $(CFLAGS) -c integer.cpp

iomanager.o:	iomanager.cpp
		$(CC) $(CFLAGS) -c iomanager.cpp

machineInfo.o:	machineInfo.cpp
		$(CC) $(CFLAGS) -c machineInfo.cpp

MbBitfield.o:	MbBitfield.cpp
		$(CC) $(CFLAGS) -c MbBitfield.cpp

MbRandom.o:	MbRandom.cpp
		$(CC) $(CFLAGS) -c MbRandom.cpp

mcmc.o:	mcmc.cpp
		$(CC) $(CFLAGS) -c mcmc.cpp

model.o:	model.cpp
		$(CC) $(CFLAGS) -c model.cpp

nexusfile.o:	nexusfile.cpp
		$(CC) $(CFLAGS) -c nexusfile.cpp

observation.o:	observation.cpp
		$(CC) $(CFLAGS) -c observation.cpp

partition.o:	partition.cpp
		$(CC) $(CFLAGS) -c partition.cpp

samples.o:	samples.cpp
		$(CC) $(CFLAGS) -c samples.cpp

stirling.o:	stirling.cpp
		$(CC) $(CFLAGS) -c stirling.cpp

clean:		
		rm -f *.o
