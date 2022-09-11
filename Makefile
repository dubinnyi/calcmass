SHELL=/bin/sh

OBJS=calcmass.o isotopes.o parse.o multiloop.o

CC       = gcc
INCLIB   = 
CPPFLAGS = 
LDFLAGS  = -lm
CFLAGS   = -O3 -g -DLINUX -D_NO_PROTO $(INCLIB)

COMMAND=calcmass

$(COMMAND):  $(OBJS)
	$(CC) $(LDFLAGS) -o $(COMMAND) $(OBJS) $(SYSLIB)

$(OBJS):    isotopes.h

clean:
	rm -f *.o *.l


