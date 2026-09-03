CXX ?= c++
CPPFLAGS ?=
CXXFLAGS ?= -O3 -Wall -Wextra
override CXXFLAGS += -std=c++17

BUILD_DIR ?= build/make
WITH_RACCESS ?= 0
RACCESS_INCLUDE_DIR ?=

CORE_NAMES := FileReader accessibility api beam_inside_outside debug debugging \
	energy_linearcapr io legacy_energy_param lincapr

ifeq ($(WITH_RACCESS),1)
ifeq ($(strip $(RACCESS_INCLUDE_DIR)),)
$(error WITH_RACCESS=1 requires RACCESS_INCLUDE_DIR=/path/to/raccess/src)
endif
CORE_NAMES += energy_raccess
override CPPFLAGS += -DLINEARRACCESS_WITH_RACCESS=1 -I$(RACCESS_INCLUDE_DIR)
override CPPFLAGS += -D_LIBCPP_ENABLE_CXX17_REMOVED_BINDERS
override CPPFLAGS += -D_LIBCPP_ENABLE_CXX17_REMOVED_RANDOM_SHUFFLE
override CPPFLAGS += -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION
endif

override CPPFLAGS += -Iinclude -Isrc

CORE_OBJECTS := $(addprefix $(BUILD_DIR)/,$(addsuffix .o,$(CORE_NAMES)))
CLI_OBJECT := $(BUILD_DIR)/linracc_main.o
TEST_OBJECT := $(BUILD_DIR)/core_test.o
DEPS := $(CORE_OBJECTS:.o=.d) $(CLI_OBJECT:.o=.d) $(TEST_OBJECT:.o=.d)

.PHONY: all test clean

all: $(BUILD_DIR)/LinRacc

test: $(BUILD_DIR)/LinRacc $(BUILD_DIR)/linearraccess_core_test
	$(BUILD_DIR)/linearraccess_core_test
	python3 tests/cli_smoke_test.py --executable $(BUILD_DIR)/LinRacc \
		$(if $(filter 1,$(WITH_RACCESS)),--with-raccess,)

$(BUILD_DIR)/LinRacc: $(CLI_OBJECT) $(CORE_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/linearraccess_core_test: $(TEST_OBJECT) $(CORE_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/core_test.o: tests/core_test.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/%.o: src/%.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	@if [ -n "$(BUILD_DIR)" ] && [ "$(BUILD_DIR)" != "/" ] && [ "$(BUILD_DIR)" != "." ]; then \
		rm -rf "$(BUILD_DIR)"; \
	fi

-include $(DEPS)
