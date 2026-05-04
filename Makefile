CXX := g++

TARGET := fluid_sim
BINDIR := bin
OBJDIR := build/obj
CLIOBJDIR := build/obj-cli

SRC := main.cpp aux.cpp core-sim-functions.cpp display-functions.cpp initializations.cpp imgui/imgui.cpp imgui/imgui_draw.cpp imgui/imgui_tables.cpp imgui/imgui_widgets.cpp imgui-sfml/imgui-SFML.cpp

OBJ := $(SRC:%.cpp=$(OBJDIR)/%.o)
DEP := $(OBJ:.o=.d)

SFML_CFLAGS ?= $(shell pkg-config --cflags sfml-graphics 2>/dev/null)
SFML_LIBS ?= $(shell pkg-config --libs sfml-graphics 2>/dev/null)

ifeq ($(strip $(SFML_LIBS)),)
	SFML_LIBS := -lsfml-graphics -lsfml-window -lsfml-system
endif

CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra
ENABLE_OPENMP ?= 1
ifeq ($(ENABLE_OPENMP),1)
	OMPFLAGS := -fopenmp
else
	OMPFLAGS :=
endif
BASE_CPPFLAGS := -I. -Iimgui -Iimgui-sfml -DIMGUI_USER_CONFIG='"imconfig-SFML.h"'
CPPFLAGS += $(BASE_CPPFLAGS)
CXXFLAGS += $(OMPFLAGS)

GUI_CPPFLAGS := $(CPPFLAGS) $(SFML_CFLAGS) -DENABLE_SFML_RENDERING
GUI_LDFLAGS := $(SFML_LIBS) -lGL $(OMPFLAGS)
CLI_CPPFLAGS := $(CPPFLAGS)
CLI_LDFLAGS := $(OMPFLAGS)

.PHONY: all build run clean help

all: build

build: $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(OBJ)
	@mkdir -p $(BINDIR)
	$(CXX) $(OBJ) -o $@ $(GUI_LDFLAGS)


# Headless CLI solver target (no runtime rendering)
CLI_TARGET := fluid_sim_cli
CLI_SRC := cli_solver.cpp aux.cpp core-sim-functions.cpp display-functions.cpp initializations.cpp
CLI_OBJ := $(CLI_SRC:%.cpp=$(CLIOBJDIR)/%.o)

$(BINDIR)/$(CLI_TARGET): $(CLI_OBJ)
	@mkdir -p $(BINDIR)
	$(CXX) $(CLI_OBJ) -o $@ $(CLI_LDFLAGS)

$(CLIOBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CLI_CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

.PHONY: cli
cli: $(BINDIR)/$(CLI_TARGET)


$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(GUI_CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

run: build
	./$(BINDIR)/$(TARGET)

help:
	@echo "Available targets:"
	@echo "  make build                 - Build the project (default)"
	@echo "  make run                   - Build and run the simulation"
	@echo "  make clean                 - Clean build artifacts"
	@echo ""
	@echo "Build options:"
	@echo "  make ENABLE_OPENMP=0       - Build without OpenMP (compatible with systems without OpenMP)"
	@echo "  make ENABLE_OPENMP=1       - Build with OpenMP enabled (default, faster on multi-core)"

clean:
	rm -rf build $(BINDIR)

-include $(DEP)
