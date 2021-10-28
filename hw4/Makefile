CC=gcc
CFLAGS=-O3 -fopenmp -lcblas

BUILD_DIR=build
SOURCE_DIR=src

SOURCES=$(wildcard $(SOURCE_DIR)/*.c)
MAIN=main.c
EXECUTABLE=main

all: $(EXECUTABLE)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(EXECUTABLE): $(BUILD_DIR) $(SOURCES)
	$(CC) $(MAIN) $(SOURCES) -o $(BUILD_DIR)/$@ $(CFLAGS)

test: all
	./$(BUILD_DIR)/$(EXECUTABLE) $(N) $(m)

clean:
	rm -rf $(BUILD_DIR)