SRC= mpfiWrapper.cpp utilities.cpp affineArithmeticClass.cpp aaEnvironment.cpp  
CXXFLAGS = -MMD -g -Wall -I /usr/local/include
CXX=g++
LDFLAGS= -L/usr/local/lib -g  -lboost_timer -lboost_system -lmpfr -lmpfi

all: filter testAffine fersonExample eulerMaruyama ssExamples

testAffine: $(SRC:.cpp=.o) testAffine.o
	g++  -g -Wall -o $@ $^ $(LDFLAGS)

filter:  $(SRC:.cpp=.o) filter.o
	g++  -g -Wall -o $@ $^ $(LDFLAGS)

fersonExample: $(SRC:.cpp=.o) fersonExample.o
	g++  -g -Wall -o $@ $^ $(LDFLAGS)


eulerMaruyama: $(SRC:.cpp=.o) eulerMaruyama.o
	g++ -g -Wall -o $@ $^ $(LDFLAGS) 

ssExamples: $(SRC:.cpp=.o) ssExamples.o
	g++  -g -Wall -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o testAffine *.d filter fersonExample eulerMaruyama ssExamples

-include $(SRC:.cpp=.d)

