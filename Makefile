all: debug

debug: clean
	gcc -g -D DEBUG -Wall -lm taco.c -o taco

release: clean
	gcc -O3 -lm taco.c -o taco
	
clean:
	rm -rf taco
