CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

SRC = main
OUT = assembler

all: $(OUT)

$(OUT): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(OUT) $(SRC)

clean:
	rm -f $(OUT)
