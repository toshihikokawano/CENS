LDFLAGS	=	-lm # -g
CPPFLAGS	=	-O3 -Wall -Wextra
CPP	=	g++
CXX	=	g++
RM      =	rm

OBJS	= cens.o censgamma.o censstat.o ensdfread.o riplread.o \
		 outxml.o outripl.o outstat.o masstable.o \
		 polysq.o polycalc.o \
		 cfgread.o

PROG	= cens

all:		$(PROG)

$(PROG):	$(OBJS)
		$(CPP) $(OBJS) $(LDFLAGS) -o $(PROG)

clean:
		$(RM) -f $(OBJS) $(PROG)

# g++ -E -MM -w *.cpp
cens.o: cens.cpp cens.h ensdf.h terminate.h elements.h cfgread.h
censgamma.o: censgamma.cpp cens.h ensdf.h terminate.h
censstat.o: censstat.cpp cens.h ensdf.h polysq.h
cfgread.o: cfgread.cpp cfgread.h
ensdfread.o: ensdfread.cpp cens.h ensdf.h terminate.h elements.h physicalconstant.h
masstable.o: masstable.cpp masstable.h masstable_audi2012_frdm2012.h
outripl.o: outripl.cpp cens.h ensdf.h elements.h masstable.h physicalconstant.h
outstat.o: outstat.cpp cens.h ensdf.h polysq.h
outxml.o: outxml.cpp cens.h ensdf.h xmltag.h
polycalc.o: polycalc.cpp polysq.h
polysq.o: polysq.cpp physicalconstant.h polysq.h
riplread.o: riplread.cpp cens.h ensdf.h terminate.h
