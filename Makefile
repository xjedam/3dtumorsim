CFLAGS = -I"freeglut\include" -L"freeglut\lib" -lfreeglut -lopengl32

3dtumorsim: 3dtumorsim.c lattice.o potts.o
	gcc 3dtumorsim.c lattice.o potts.o -o bin/3dtumorsim $(CFLAGS)

lattice.o: lattice.c
	gcc -c lattice.c $(CFLAGS)
  
potts.o: potts.c
	gcc -c potts.c $(CFLAGS)
  
clean:
	rm *.o

