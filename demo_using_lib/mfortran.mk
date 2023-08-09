
ifeq ($(OS),Windows_NT)
	LIBFLAGES= -L../target/release/ -lif97  
	EXEDIR=../target/release/
else
	UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        LIBFLAGES= -L../target/release -Wl,-rpath=../target/release  -lif97 -lm
		EXEDIR=./
    endif
endif

all: 
	gfortran -fno-underscoring demo.f08 -o$(EXEDIR)demo  $(LIBFLAGES)
	$(EXEDIR)demo



 