CXX=g++
CFLAGS=-c
CXXFLAGS=-std=c++17


BUILD_DIR=build
SOURCE_DIR=src

SOURCES=$(wildcard $(SOURCE_DIR)/*.cpp)
OBJECTS=$(patsubst $(SOURCE_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
MAIN=test/test.cpp
MATRIX_LIB=matrixlib.a
EXECUTABLE=test

all: $(BUILD_DIR) $(MATRIX_LIB) $(EXECUTABLE)

$(BUILD_DIR):
	mkdir -p $@

$(MATRIX_LIB): $(BUILD_DIR)
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(OPTIMIZE) $(SOURCES) -o $(OBJECTS)
	ar rc $(MATRIX_LIB) $(OBJECTS)

$(EXECUTABLE): $(MATRIX_LIB) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -I./ $(MAIN) $< -o $(BUILD_DIR)/$@

run:
	./$(BUILD_DIR)/$(EXECUTABLE) $(N)

clean:
	rm -rf $(BUILD_DIR) *.a