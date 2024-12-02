ARCH := $(shell arch)
ifeq ($(ARCH),x86_64)
CFLAGS=-Wall -O3 -static
LIBS=-lhts -lpthread -ldeflate -lisal -lgsl -lgslcblas -lcairo -lfontconfig -lpixman-1 -lfreetype -lexpat -lzip -lbz2 -llzma -lpng -lm -ldl -lz
else
CFLAGS=-Wall -O3 -march=native
LIBS=-lhts -lisal -lgsl -lgslcblas -lm -lcairo
endif

CC=cc

SRC = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(SRC))

bamgd: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

test: bamgd
	./$^ -i ../dat/573.bam -o ../dat/573.png -r ../dat/573.fa.gz -d

clean:
	rm -f bamgd *.o ../dat/573.png
