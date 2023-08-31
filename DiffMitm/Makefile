#Polynomial used to define the finite field GF(2^n) -
# AES Field
GF_POLY := 0x11b
# F_2
#GF_POLY := 0x02
# Piccolo
#GF_POLY := 0x13

GRBPATH = /opt/gurobi
USERDEFINES := -DGF_POLY=$(GF_POLY)

LDFLAGS := -L$(GRBPATH)/linux64/lib -lgurobi_c++ -lgurobi95 -ltbb -lm
#LDFLAGS := -L/Library/gurobi952/macos_universal2/lib -lgurobi_c++ -lgurobi95

CXX      := -g++
CXXFLAGS := -Wall -Wextra -std=c++11
#LDFLAGS  := -L/usr/lib -lstdc++ -lm
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/apps
TARGET   := program
INCLUDE  := -Iinclude/
SRC      :=                      \
   $(wildcard src/*.cpp)         \

OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES \
         := $(OBJECTS:.o=.d)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(USERDEFINES) $(INCLUDE) -c $< -MMD -o $@

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(USERDEFINES) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

-include $(DEPENDENCIES)

.PHONY: all build clean debug release info

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O3 -march=native -DNDEBUG
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*

info:
	@echo "[*] Application dir: ${APP_DIR}     "
	@echo "[*] Object dir:      ${OBJ_DIR}     "
	@echo "[*] Sources:         ${SRC}         "
	@echo "[*] Objects:         ${OBJECTS}     "
	@echo "[*] Dependencies:    ${DEPENDENCIES}"
