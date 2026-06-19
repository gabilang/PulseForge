# BioSignalProcessor — plain Makefile (no cmake required).
CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Iinclude
BUILD    := build

# Some macOS Command Line Tools installs ship an incomplete toolchain-local
# libc++; fall back to the SDK's headers if <complex> isn't otherwise found.
ifeq ($(shell uname),Darwin)
  SDK_CXX := $(shell xcrun --show-sdk-path 2>/dev/null)/usr/include/c++/v1
  ifneq ($(wildcard $(SDK_CXX)/complex),)
    CXXFLAGS += -isystem $(SDK_CXX)
  endif
endif

LIB_SRCS := $(wildcard src/*.cpp src/io/*.cpp src/bio/*.cpp src/viz/*.cpp)
LIB_OBJS := $(patsubst src/%.cpp,$(BUILD)/%.o,$(LIB_SRCS))

APP_BIN   := $(BUILD)/pulseforge
TEST_BIN  := $(BUILD)/bsp_tests

.PHONY: all app test clean run demo

all: app test

app: $(APP_BIN)

$(APP_BIN): $(LIB_OBJS) $(BUILD)/app/main.o
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "Built $(APP_BIN)"

$(TEST_BIN): $(LIB_OBJS) $(BUILD)/tests/test_main.o
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "Built $(TEST_BIN)"

test: $(TEST_BIN)
	./$(TEST_BIN)

run: $(APP_BIN)
	./$(APP_BIN) demo

demo: run

# Generic compile rule for any .cpp under src/, app/, tests/.
$(BUILD)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD)/app/%.o: app/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD)/tests/%.o: tests/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD)
