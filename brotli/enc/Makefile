#brotli/enc

include ../../shared.mk

OBJS = backward_references.o block_splitter.o encode.o entropy_encode.o histogram.o literal_cost.o prefix.o

all : $(OBJS)

encode : $(OBJS)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(OBJS) main.cc -o encode

clean :
	rm -f $(OBJS) $(SO)

