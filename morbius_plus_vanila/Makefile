.PHONY: all test clean

CXX := g++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra -pedantic -pthread
TARGET := morbius_plus_vanila
SOURCES := main.cpp Utility.cpp SeedInitialization.cpp GibbsPipeline.cpp Result.cpp Refinement.cpp RefinementOutput.cpp
HEADERS := MorbiusPlus.h Refinement.h

all: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

test: $(TARGET)
	bash test/run_test.sh

clean:
	rm -f $(TARGET)
	rm -rf test/output test/output_san
