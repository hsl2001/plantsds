CC = gcc
CFLAGS = -O3 -march=native -Wall -Wextra
LDFLAGS = -lm -lpthread -lz
TARGET = segtrace
SRCS = segtrace.c klib/kthread.c

.PHONY: clean static

$(TARGET): $(SRCS) segtrace.h
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)

all: $(SRCS) segtrace.h
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)

static: $(SRCS) segtrace.h
	$(CC) $(CFLAGS) -static -o $(TARGET) $(SRCS) $(LDFLAGS)

clean:
	rm -f $(TARGET)