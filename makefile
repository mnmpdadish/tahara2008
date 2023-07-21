
all: 
	cd pfapack && $(MAKE)
	gcc -Wall -O2 test.c -lm -llapack -lblas -lgfortran ./pfapack/libpfapack.a

clean:
	rm a.out
