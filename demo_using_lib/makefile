
ifeq ($(OS),Windows_NT)
	LIBFLAGES= -L../target/release/ -lseuif97  
	EXEDIR=../target/release/
else
	UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        LIBFLAGES= -L../target/release -Wl,-rpath=../target/release  -lseuif97 -lm
		EXEDIR=./
    endif
endif

all: 
	gcc demo.c -o$(EXEDIR)demo  $(LIBFLAGES)
	$(EXEDIR)demo



 