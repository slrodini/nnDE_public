
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
PROGRAM=generic_pot_D2
PROGRAMD1=generic_pot_D1
TESTPROGRAM = test_der_dim

CSRC=$(wildcard $(SRCDIR)/*.c)
#OBJS = $(addprefix $(OBJDIR)/,$(CSRC:.c=.o))
#The upper one ok if CSRC has no path in front!
OBJS=$(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(CSRC))

all: $(BINDIR)/$(PROGRAM) $(BINDIR)/$(PROGRAMD1)
d2der: $(BINDIR)/$(PROGRAM)
d1der: $(BINDIR)/$(PROGRAMD1)
test: $(BINDIR)/$(TESTPROGRAM)

$(BINDIR)/$(PROGRAM): $(OBJS) $(OBJDIR)/$(PROGRAM).o
	$(CC) $(CFLAGS) $(OBJS) $(OBJDIR)/$(PROGRAM).o -o $@ $(LDFLAGS)

$(BINDIR)/$(PROGRAMD1): $(OBJS) $(OBJDIR)/$(PROGRAMD1).o
	$(CC) $(CFLAGS) $(OBJS) $(OBJDIR)/$(PROGRAMD1).o -o $@ $(LDFLAGS)

$(BINDIR)/$(TESTPROGRAM): $(OBJS) $(OBJDIR)/$(TESTPROGRAM).o
	$(CC) $(CFLAGS) $(OBJS) $(OBJDIR)/$(TESTPROGRAM).o -o $@ $(LDFLAGS)

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

$(OBJDIR)/$(PROGRAM).o : $(PRODIR)/$(PROGRAM).c
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

$(OBJDIR)/$(PROGRAMD1).o : $(PRODIR)/$(PROGRAMD1).c
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

$(OBJDIR)/$(TESTPROGRAM).o : $(TESTPRODIR)/$(TESTPROGRAM).c
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

clean:
	rm -rf $(BINDIR)/* $(OBJDIR)/* && clear
