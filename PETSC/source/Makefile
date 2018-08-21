ALL: heat_petsc

CFLAGS		 = -I${PETSC_DIR}/../include -I${PETSC_DIR}/include
FFLAGS		 = 
CPPFLAGS         =
FPPFLAGS         =

include ${PETSC_DIR}/lib/petsc/conf/petscvariables
include ${PETSC_DIR}/lib/petsc/conf/petscrules

heat_petsc: heat_petsc.o
	${CC} -o heat_petsc heat_petsc.o ${PETSC_SNES_LIB}
	${RM} heat_petsc.o
