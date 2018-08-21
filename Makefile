#INCLUDES = -I/usr/local/include -I./include -I. # include from these three dirs. # ??? it does not inclue ./include/nlib.h   ...why?
#LIBS = -L./libmisc -lmisc # ??? what is library?

CXX = g++ # compiler
#CXXFLAGS = -std=c++0x -Wall -g -O2 -pipe 
CXXFLAGS =  -Wall -g -O2 -pipe 
TARGET = a.out # the name of the program you want to make
SRCS = main.cpp majorminor.cpp pdbnishi.cpp tranishi.cpp centnishi.cpp inpnishi.cpp math_nishi.cpp# the names of all source files
DEPS = nlib.h math_nishi.h # header; a dependency of *.o files
OBJS := $(SRCS:.cpp=.o) # change the suffix from cpp to o

$(TARGET): $(OBJS)
	g++ -Wall -o $@ $(OBJS)
$(OBJS): $(DEPS)	# enable rebuilding with changing header
#main.o: main.cpp
#	g++ -Wall -c main.cpp
#sub.o: sub.cpp
#	g++ -Wall -c sub.cpp
clean:
	rm -f $(TARGET) $(OBJS)

