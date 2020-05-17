SHELL := /bin/bash

CC        := gcc
CFLAGS    := -std=gnu99 -static -Wall -g

BIN     := bin
SRC     := src workspace
INCLUDE := include

LIBRARIES   := -lgsl -lgslcblas -lm

EXECUTABLE  := odesolvers

SOURCEDIRS  := $(shell find $(SRC) -maxdepth 0 -type d)
INCLUDEDIRS := $(shell find $(INCLUDE) -maxdepth 0 -type d)

CINCLUDES   := $(patsubst %,-I %, $(INCLUDEDIRS:%/=%))
SOURCES     := $(wildcard $(patsubst %, %/*.c, $(SOURCEDIRS)))


all: $(BIN)/$(EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(SOURCES)
	$(CC) $(CFLAGS) $(CINCLUDES) $^ -o $@ $(LIBRARIES)

run: all
	./$(BIN)/$(EXECUTABLE)

debug:
	gdb ./$(BIN)/$(EXECUTABLE)

.PHONY: clean

clean:
	-$(RM) $(BIN)/$(EXECUTABLE)
