SHELL=/bin/sh

OBJS=calcmass.o isotopes.o parse.o multiloop.o

CC       = gcc
INCLIB   = 
CPPFLAGS = 
CFLAGS   = -O3 -g -DLINUX -D_NO_PROTO $(INCLIB)

COMMAND=calcmass

$(COMMAND):  $(OBJS)
	$(CC) -o $(COMMAND) $(OBJS) $(SYSLIB)

$(OBJS):    isotopes.h

clean:
	rm -f *.o *.l


