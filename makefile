CC	= g++
CPPFLAGS= -Wno-deprecated -O2 -c
LDFLAGS	= -O2
SOURCES	= main.cpp Util.cpp BiCover.cpp Distance.cpp Graph.cpp GraphUtil.cpp 
OBJECTS	= $(SOURCES:.cpp=.o)
EXECUTABLE=distance

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

rm:
	rm -f *.o
clean:
	rm -f *.o
