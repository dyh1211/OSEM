FLAGS=-g -O6
LIBS=-lm
OBJECTS=fdtd.o  fdtd_sparameters.o fdtd_abc.o fdtd_mesh.o fdtd_updates.o  fdtd_farfield.o  fdtd_pml.o


fdtd:	$(OBJECTS)
	g++ -o $@ $(OBJECTS) $(LIBS) $(FLAGS)

%.o:	%.c
	g++ $(FLAGS) -o $@ -c $<

clean:
	rm -f *.o fdtd *.c~ *.h~ *.cpp~ *.hpp~ Makefile~
