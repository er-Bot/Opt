COMP = g++
FLAGS = -std=c++0x -Wall

LIB_TARGET = nlo.a
TST_TARGET = nlo

SRC = $(shell find ./src -name "*.cpp")
OBJ = $(patsubst ./src/%.cpp, ./obj/%.o, $(SRC))

# run lib and test
.phony: all
all: lib test

.phony: dirs
dirs:
	@mkdir -p obj lib

# first make folders obj and lib
.phony: lib
lib: dirs $(OBJ)
	@mkdir -p ./lib
	@echo [AR] ./lib/$(LIB_TARGET)
	@ar -crs ./lib/$(LIB_TARGET) $(OBJ)

# create the .o object files
obj/%.o: src/%.cpp
	@echo [COMP] $@
	@$(COMP) -c $< -o $@ $(FLAGS) -Iinclude

.phony: test
test:
	@mkdir -p bin
	@echo [*] Compiling the test ..
	@$(COMP) ./main.cpp ./lib/$(LIB_TARGET) -o bin/$(TST_TARGET) $(CCFLAGS) -Iinclude
	
.phony: clean
clean:
	@rm -r obj lib bin
	