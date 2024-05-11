CC = g++
CPPFLAGS = -O3 -I/home/jovyan/lapackpp/include -I/home/jovyan/include -I.
LDFLAGS = -L/home/jovyan/lapackpp/lib -L/home/jovyan/lib -L.
LDLIBS = -llapack -lblas -llapackpp -fopenmp -lintegrals

DEPENDENCIES_DIR = dependencies
INCLUDE = include

CPP_FILES = $(shell grep -l -L 'int main' *.cpp)
OBJECT = $(CPP_FILES:.cpp=.o)
EXECUTABLES = $(shell grep -l 'int main' *.cpp | sed 's/.cpp//')
DEPENDENCIES = $(OBJECT:.o=.d)

%.d: %.cpp
	@set -e; rm -f $@; \
    $(CC) -M $(CPPFLAGS) $< > $@.$$$$; \
    sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
    rm -f $@.$$$$; \
	mkdir -p $(DEPENDENCIES_DIR); \
	mv $@ $(DEPENDENCIES_DIR)

all: $(EXECUTABLES)
	@set -e; mkdir -p $(INCLUDE); mv $(OBJECT) $(INCLUDE)

$(EXECUTABLES): $(OBJECT)

include $(DEPENDENCIES)

clean:
	rm -r $(DEPENDENCIES_DIR) $(INCLUDE); rm $(EXECUTABLES)

install: 

.PHONY: clean install