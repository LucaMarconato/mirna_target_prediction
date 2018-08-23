INCLUDE = -I./dependencies/

SRCS = main.cpp mirna.cpp matching.cpp

build: $(SRCS)
	g++ -std=c++11 $(SRCS) $(INCLUDE)

run: build
	./a.out

