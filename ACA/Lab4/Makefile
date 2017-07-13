BIN = smt
SRC = smt.c
OBJ = $(SRC:.c=.o)
CFLAGS = -O2 -g -Wall -std=c99 -D_GNU_SOURCE
LIBS = -lrt -lpthread

all: $(BIN)

$(BIN):$(OBJ)
	gcc $(CFLAGS) $(OBJ) -o $(BIN) $(LIBS)

%.o:%.c
	gcc $(CFLAGS) -c $<

clean:
	rm -f $(BIN) $(OBJ)

