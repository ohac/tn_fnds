OBJS = dio.o matlabfunctions.o platinum.o star.o synthesis.o test.o fft.o wavread.o
CXX = g++
CFLAGS = -O
INCLUDES = -I.
LIBS = -lfftw3

tn_fnds : $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<

clean:
	rm tn_fnds $(OBJS)
