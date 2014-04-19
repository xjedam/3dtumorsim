CFLAGS = -I"freeglut\include" -L"freeglut\lib" -lfreeglut -lopengl32

3dtumorsim: 3dtumorsim.c lattice.o
	gcc 3dtumorsim.c lattice.o -o bin/3dtumorsim $(CFLAGS)

lattice.o: lattice.c
	gcc -c lattice.c $(CFLAGS)
  
clean:
	rm *.o

