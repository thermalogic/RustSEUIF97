ifeq ($(OS),Windows_NT)
	LDFLAGS= -L../target/release/ -lseuif97  
	EXEDIR=../target/release/
	CGO_SET=set CGO_LDFLAGS=$(LDFLAGS)&&
else
	UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        LDFLAGS= -L../target/release -Wl,-rpath,../target/release -lseuif97 -lm
		EXEDIR=./
		CGO_SET=CGO_LDFLAGS="$(LDFLAGS)"
    endif
endif

all: 
	$(CGO_SET) go build -o $(EXEDIR)demo demo.go
	$(EXEDIR)demo