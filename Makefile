CC      = g++-9 
CFLAGS  = -g -Wall -fopenmp 
EXEC    = dragoncello
LIB     = libdragoncello.a
AR      = ar
ARFLAGS = rcs

GSL_DIR = /usr/local/Cellar/gsl/2.5

INCDIR += -I$(GSL_DIR)/include
LIBDIR += -L$(GSL_DIR)/lib -lgsl -lgslcblas

OBJS += algorithm.o buildgrids.o dragoncello.o run2D.o 

all: $(EXEC) $(LIB)

$(EXEC): $(OBJS) main.o
	$(CC) $(CFLAGS) $(INCDIR) -o $@ $^ $(LIBDIR)

%.o: %.cpp constants.h
	$(CC) $(CFLAGS) $(INCDIR) -c -o $@ $< 	

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

.PHONY: clean

clean:
	@rm -vf *.o
	@rm -vf $(EXEC)
	@rm -vf $(LIB)
