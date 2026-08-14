CXX := g++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra -pedantic
TARGET := synthethic_dataset_generator
SOURCES := main.cpp Utility.cpp DatasetPreparation.cpp DatasetGeneration.cpp
OBJECTS := $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)

%.o: %.cpp SyntheticDatasetGenerator.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

test: $(TARGET)
	bash test/run_test.sh

clean:
	rm -f $(TARGET) $(OBJECTS)
	rm -rf test/input test/output

.PHONY: all test clean
