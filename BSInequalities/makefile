CC = g++ # compiler
.PHONY: all compile clean
all: compile run

# Cpp file names
file1 = Driver
file2 = AbstractFunctions
file3 = Functions
file4 = Option


compile: ${file1}.cpp ${file2}.cpp ${file3}.cpp ${file4}.cpp
	${CC} -g -c -o ${file1}.o ${file1}.cpp
	${CC} -g -c -o ${file2}.o ${file2}.cpp
	${CC} -g -c -o ${file3}.o -std=c++11 ${file3}.cpp
	${CC} -g -c -o ${file4}.o ${file4}.cpp
	${CC} -g -o Execute ${file1}.o ${file2}.o ${file3}.o ${file4}.o

run: Execute
	./Execute
	rm *.o Execute

clean:
	rm *.o Execute
