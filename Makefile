EXEC_NAME=FASE
CC=g++

# CFLAGS= -Wall -Wno-write-strings -O0 -g
#CFLAGS= -w -Wno-write-strings -O3 -g
CFLAGS= -w -Wno-write-strings -O3 -std=c++11

SRC =                   	\
	Fase.cpp 		\
	Timer.cpp		\
	GTrie.cpp		\
	GraphMatrix.cpp		\
	Graphlets5.cpp		\
	GraphList.cpp		\
	GraphUtils.cpp		\
	Isomorphism.cpp		\
	Random.cpp		\
	nauty/nauty.c		\
	nauty/nautil.c		\
	nauty/naugraph.c 	\
	FaseMain.cpp

OBJ =  ${SRC:.cpp=.o}

#------------------------------------------------------------

all: ${EXEC_NAME}

${EXEC_NAME}: ${OBJ}
	${CC} ${CFLAGS} ${CLIBS} -o ${EXEC_NAME} ${OBJ}

%.o: %.cpp
	${CC} ${CFLAGS} -c -o $@ $+

clean:
	rm ${EXEC_NAME} *.o *~ *# -rf
