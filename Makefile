CHROMA=${HOME}/bin/pdfchroma/bin/
CONFIG=$(CHROMA)chroma-config
CXX=$(shell $(CONFIG) --cxx)
CXXFLAGS=$(shell $(CONFIG) --cxxflags) -I. -I${HOME}/bin/qopqdp/include -g 
LDFLAGS=-L${HOME}/bin/pdfchroma/lib $(shell $(CONFIG) --ldflags) -L${HOME}/bin/qopqdp/lib -g
LIBS=$(shell $(CONFIG) --libs)   -lwilsonmg -lqopqdp

GSL=/global/common/sw/cray/cnl7/haswell/gsl/2.5/gcc/8.2.0/sr445ay/bin/
GSL_CONFIG=$(GSL)gsl-config
GSL_CFLAGS=$(shell $(GSL_CONFIG) --cflags) 
GSL_LIBS=$(shell $(GSL_CONFIG) --libs)

QLA=${HOME}/bin/qla/bin/
QLA_CONFIG=$(QLA)qla-config
QLA_CFLAGS=$(shell $(QLA_CONFIG) --cflags) 
QLA_LIBS=$(shell $(QLA_CONFIG) --libs)
QLA_LDFLAGS=$(shell $(QLA_CONFIG) --ldflags)

QDPC=${HOME}/bin/qdp/bin/
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
