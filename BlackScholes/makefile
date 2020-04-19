CC = g++ # compiler
.PHONY: all compile clean
all: compile run

# Cpp file names
file1 = Parabolic
file2 = Function1D
file3 = Function2D
file4 = InitialU
file5 = ExactU
file6 = Driver

compile: ${file1}.cpp ${file2}.cpp ${file3}.cpp ${file4}.cpp ${file5}.cpp ${file6}.cpp
	${CC} -g -c -o ${file1}.o ${file1}.cpp
	${CC} -g -c -o ${file2}.o ${file2}.cpp
	${CC} -g -c -o ${file3}.o ${file3}.cpp
	${CC} -g -c -o ${file4}.o -std=c++11 ${file4}.cpp # Dervied class
	${CC} -g -c -o ${file5}.o -std=c++11 ${file5}.cpp
	${CC} -g -c -o ${file6}.o ${file6}.cpp
	${CC} -g -o Execute ${file1}.o ${file2}.o ${file3}.o ${file4}.o ${file5}.o ${file6}.o

run: Execute
	./Execute
	rm *.o Execute

clean:
	rm *.o Execute
