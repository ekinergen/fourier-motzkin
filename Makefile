all: main.cpp
	mkdir bin
	g++ -Wall -Werror -Wextra -Wpedantic -o bin/main main.cpp
clean:
	rm -r bin
