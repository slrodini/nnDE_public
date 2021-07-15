
BASE := $(PWD)

SRCDIR=$(BASE)/src
OBJDIR=$(BASE)/obj
INCDIR=$(BASE)/inc
BINDIR=$(BASE)/bin

CC=gcc
CFLAGS = -O3 -Wall -I$(BASE)/inc
PROGRAM=e6502

CSRC=$(wildcard $(SRCDIR)/*.c)
#OBJS = $(addprefix $(OBJDIR)/,$(CSRC:.c=.o))
#The upper one ok if CSRC has no path in front!
OBJS=$(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(CSRC))

all: $(BINDIR)/$(PROGRAM)

$(BINDIR)/$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@

$(OBJDIR)/%.o : $(SRCDIR)/%.c $(INCDIR)/%.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BINDIR)/* $(OBJDIR)/* && clear
