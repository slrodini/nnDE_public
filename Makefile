
BASE := $(PWD)

SRCDIR=$(BASE)/src
OBJDIR=$(BASE)/obj
INCDIR=$(BASE)/inc
BINDIR=$(BASE)/bin


CC=gcc
CFLAGS = -O3 -Wall -I$(BASE)/inc -O3 
LDFLAGS = -lm -lgsl -lgslcblas
PRODIR=$(BASE)/runs
TESTPRODIR=$(BASE)/tests
PROGRAM=generic_pot
TESTPROGRAM = generic_pot

CSRC=$(wildcard $(SRCDIR)/*.c)
#OBJS = $(addprefix $(OBJDIR)/,$(CSRC:.c=.o))
#The upper one ok if CSRC has no path in front!
OBJS=$(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(CSRC))

all: $(BINDIR)/$(PROGRAM)
test: $(BINDIR)/$(TESTPROGRAM)

$(BINDIR)/$(PROGRAM): $(OBJS) $(OBJDIR)/$(PROGRAM).o
	$(CC) $(CFLAGS) $(OBJS) $(OBJDIR)/$(PROGRAM).o -o $@ $(LDFLAGS)

$(BINDIR)/$(TESTPROGRAM): $(OBJS) $(OBJDIR)/$(TESTPROGRAM).o
	$(CC) $(CFLAGS) $(OBJS) $(OBJDIR)/$(TESTPROGRAM).o -o $@ $(LDFLAGS)

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

$(OBJDIR)/$(PROGRAM).o : $(PRODIR)/$(PROGRAM).c
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

$(OBJDIR)/$(TESTPROGRAM).o : $(TESTPRODIR)/$(TESTPROGRAM).c
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

clean:
	rm -rf $(BINDIR)/* $(OBJDIR)/* && clear
