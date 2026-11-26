# @file Makefile
# @brief Build, test, and install configuration for the calcmass mass-spectrometry CLI utility.
#
# Targets:
#   calcmass (default) - Compiles and links the executable with strict C11 warnings.
#   test               - Builds and runs C unit tests and CLI integration tests.
#   debug              - Builds with -O0 -g and Address/UndefinedBehavior sanitizers.
#   install            - Installs the executable to $(DESTDIR)$(PREFIX)/bin.
#   uninstall          - Removes the installed executable from $(DESTDIR)$(PREFIX)/bin.
#   clean              - Removes compiled executables and intermediate object files.
#
# Note: $(LDLIBS) (-lm) is placed after $(OBJS) to satisfy linkers with --as-needed.

SHELL    = /bin/sh

OBJS     = calcmass.o isotopes.o parse.o multiloop.o
LIB_OBJS = isotopes.o parse.o multiloop.o

CC       ?= gcc
PREFIX   ?= /usr/local
INCLIB   =
CPPFLAGS =
CFLAGS   = -std=c11 -O3 -g -Wall -Wextra -Wpedantic -Wshadow -DLINUX $(INCLIB)
LDFLAGS  =
LDLIBS   = -lm

COMMAND  = calcmass
UNIT_BIN = tests/test_unit

$(COMMAND): $(OBJS)
	$(CC) $(LDFLAGS) -o $(COMMAND) $(OBJS) $(LDLIBS) $(SYSLIB)

$(OBJS): isotopes.h

$(UNIT_BIN): tests/test_unit.c $(LIB_OBJS) isotopes.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $(UNIT_BIN) tests/test_unit.c $(LIB_OBJS) $(LDLIBS) $(SYSLIB)

test: $(COMMAND) $(UNIT_BIN)
	./$(UNIT_BIN)
	sh tests/test_cli.sh ./$(COMMAND)

debug:
	$(MAKE) clean
	$(MAKE) CFLAGS="-std=c11 -O0 -g -Wall -Wextra -Wpedantic -Wshadow -fsanitize=address,undefined -DLINUX $(INCLIB)" LDFLAGS="-fsanitize=address,undefined" $(COMMAND)

install: $(COMMAND)
	install -d "$(DESTDIR)$(PREFIX)/bin"
	install -m 755 $(COMMAND) "$(DESTDIR)$(PREFIX)/bin/$(COMMAND)"

uninstall:
	rm -f "$(DESTDIR)$(PREFIX)/bin/$(COMMAND)"

clean:
	rm -f $(COMMAND) $(UNIT_BIN) *.o *.l tests/*.o

.PHONY: test debug install uninstall clean
