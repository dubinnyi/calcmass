# @file Makefile
# @brief Build configuration for the calcmass mass-spectrometry CLI utility.
#
# Targets:
#   calcmass (default) - Compiles and links the executable with strict C11 warnings.
#   clean              - Removes the compiled executable and intermediate object files.
#
# Note: $(LDLIBS) (-lm) is placed after $(OBJS) to satisfy linkers with --as-needed.

SHELL    = /bin/sh

OBJS     = calcmass.o isotopes.o parse.o multiloop.o

CC       = gcc
INCLIB   =
CPPFLAGS =
CFLAGS   = -std=c11 -O3 -g -Wall -Wextra -Wpedantic -Wshadow -DLINUX $(INCLIB)
LDFLAGS  =
LDLIBS   = -lm

COMMAND  = calcmass

$(COMMAND): $(OBJS)
	$(CC) $(LDFLAGS) -o $(COMMAND) $(OBJS) $(LDLIBS) $(SYSLIB)

$(OBJS): isotopes.h

clean:
	rm -f $(COMMAND) *.o *.l

.PHONY: clean


