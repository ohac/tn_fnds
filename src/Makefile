OBJS = dio.o matlabfunctions.o platinum.o star.o synthesis.o test.o fft.o wavread.o
CXX = g++
CFLAGS = -O -msse3
INCLUDES = -I.
LIBS = -lfftw3

tn_fnds : $(OBJS) fft.h mmsystem.h windows.h wavread.h world.h
	$(CXX) -o $@ $(OBJS) $(LIBS) $(CFLAGS)

.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<

clean:
	rm tn_fnds $(OBJS)
