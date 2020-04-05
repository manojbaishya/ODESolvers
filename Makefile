CC        := gcc
CFLAGS  := -std=gnu99 -Wall -g -fstack-protector

BIN     := bin
SRC     := src
INCLUDE := include ~/.local/include/vcpkg/installed/x64-linux/include
LIB     := lib ~/.local/include/vcpkg/installed/x64-linux/lib
OBJ     := objs
IO      := iodata

LIBRARIES   := -lm -lgslcblas -lgsl

EXECUTABLE  := main

SOURCEDIRS  := $(shell find $(SRC) -maxdepth 0 -type d)
OBJDIRS     := $(shell find $(OBJ) -maxdepth 0 -type d)
INCLUDEDIRS := $(shell find $(INCLUDE) -maxdepth 0 -type d)
LIBDIRS     := $(shell find $(LIB) -maxdepth 0 -type d)
IODIR       := $(shell find $(IO) -maxdepth 0 -type d)

CINCLUDES   := $(patsubst %,-I%, $(INCLUDEDIRS:%/=%))
CLIBS       := $(patsubst %,-L%, $(LIBDIRS:%/=%))

SOURCES     := $(wildcard $(patsubst %, %/*.c, $(SOURCEDIRS)))
OBJECTS     := $(patsubst $(SOURCEDIRS)/%.c, $(OBJDIRS)/%.o, $(SOURCES))
IODATA      := $(wildcard $(patsubst %, %/*, $(IODIR:%/=%)))


all: $(BIN)/$(EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(OBJECTS)
	$(CC) -o $@ $^ $(CLIBS) $(LIBRARIES)
	@printf "\n================= Built All Executables ===================\n"

$(OBJECTS): $(OBJDIRS)/%.o: $(SOURCEDIRS)/%.c
	$(CC) $(CFLAGS) $(CINCLUDES) -c -o $@ $<

run: all
	@printf "\n---------------------- Running File -----------------------\n"
	./$(BIN)/$(EXECUTABLE)

debug:
	gdb ./$(BIN)/$(EXECUTABLE)

.PHONY: clean

clean:
	@printf "\n***************** Erasing Residual Files ******************\n\n"
	-$(RM) $(BIN)/$(EXECUTABLE)
	-$(RM) $(OBJECTS)
ifneq ($(IODATA),)
	rm -rf $(IODATA)
endif
	@printf "\n************************ Success **************************\n\n"

