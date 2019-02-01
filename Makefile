ex4: ex4.o Bitmap.o
	gcc -o $@ $^
ex3: ex3.o Bitmap.o
	gcc -o $@ $^
ex2: ex2.o Bitmap.o
	gcc -o $@ $^
ex1: ex1.o Bitmap.o
	gcc -o $@ $^
.c.o:
	gcc -c $<
