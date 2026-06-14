
ifeq ($(OS),Windows_NT)
	LIBFLAGS=-L../../target/x86_64-pc-windows-gnu/release  -lseuif97  
	EXEDIR=../../target/x86_64-pc-windows-gnu/release/
else
	UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        LIBFLAGS= -L../../target/release -Wl,-rpath=../../target/release  -lseuif97 -lm
		EXEDIR=./
    endif
endif

all: 
	gcc demo.c -o$(EXEDIR)demo  $(LIBFLAGS)
	$(EXEDIR)demo



 