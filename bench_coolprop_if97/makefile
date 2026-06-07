
ifeq ($(OS),Windows_NT)
	LIBFLAGS= -L../target/x86_64-pc-windows-gnu/release/ -lseuif97  
	EXEDIR=../target/x86_64-pc-windows-gnu/release/
else
	UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        LIBFLAGS= -L../target/release -Wl,-rpath=../target/release  -lseuif97 -lm	
		EXEDIR=./
    endif
endif

all: 
	g++ bench_if97.cpp -o$(EXEDIR)bench_if97 $(LIBFLAGS) 
	$(EXEDIR)bench_if97

run: 
	$(EXEDIR)bench_if97

save: 
	g++ bench_if97.cpp -o$(EXEDIR)bench_if97 $(LIBFLAGS) 
	$(EXEDIR)bench_if97 > benchmark_results.txt

clean:
	rm -f $(EXEDIR)bench_if97 benchmark_results.txt