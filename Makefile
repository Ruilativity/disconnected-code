CHROMA=/work/06332/zhangrui/stampede2/bin/chroma-double_s4/bin/
CONFIG=$(CHROMA)chroma-config
CXX=$(shell $(CONFIG) --cxx)
CXXFLAGS=$(shell $(CONFIG) --cxxflags) -I. -g
LDFLAGS=-L/work/06332/zhangrui/stampede2/bin/chroma-double_s4/lib $(shell $(CONFIG) --ldflags) -g
LIBS=$(shell $(CONFIG) --libs) 

GSL=/opt/apps/gcc7_1/gsl/2.3/bin/
GSL_CONFIG=$(GSL)gsl-config
GSL_CFLAGS=$(shell $(GSL_CONFIG) --cflags) 
GSL_LIBS=$(shell $(GSL_CONFIG) --libs)

QLA=/work/06332/zhangrui/stampede2/bin/qla/bin/
QLA_CONFIG=$(QLA)qla-config
QLA_CFLAGS=$(shell $(QLA_CONFIG) --cflags) 
QLA_LIBS=$(shell $(QLA_CONFIG) --libs)
QLA_LDFLAGS=$(shell $(QLA_CONFIG) --ldflags)

QDPC=/work/06332/zhangrui/stampede2/bin/qdp/bin/
QDPC_CONFIG=$(QDPC)qdp-config
QDPC_CFLAGS=$(shell $(QDPC_CONFIG) --cflags) 
QDPC_LIBS=$(shell $(QDPC_CONFIG) --libs)
QDPC_LDFLAGS=$(shell $(QDPC_CONFIG) --ldflags)

SRCS=$(wildcard *.cc)
OBJS=$(SRCS:.cc=.o)
EXES=$(SRCS:.cc=)

all: $(OBJS) $(EXES)

%: %.o
	$(CXX) -o $@ $(CXXFLAGS) $< $(LDFLAGS) $(LIBS) $(GSL_LIBS) $(QDPC_LIBS) $(QDPC_LDFLAGS) $(QDPC_LIBS) $(QDPC_LDFLAGS) $(LDFLAGS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(GSL_CFLAGS) $(QLA_CFLAGS) $(QDPC_CFLAGS) -c $<

clean:
	rm -rf $(OBJS) $(EXES) *~
