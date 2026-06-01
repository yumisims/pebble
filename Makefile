CC      = cc
CFLAGS  = -std=c11 -Wall -Wextra -Wpedantic -Iinclude
LDFLAGS =
ifeq ($(shell uname -s),Linux)
LDFLAGS += -static
endif
LDLIBS  = -lm

PEBBLE_BIGWIG ?= 0

BUILD_DIR = build
LIB       = $(BUILD_DIR)/libpebble.a
CLI       = $(BUILD_DIR)/pebble
TESTS     = $(BUILD_DIR)/pebble_tests
IO_TESTS  = $(BUILD_DIR)/pebble_io_tests

LIB_SRCS  = src/pebble.c src/pebble_io.c

ifeq ($(PEBBLE_BIGWIG),1)
CFLAGS += -DPEBBLE_BIGWIG -DNOCURL
LIBBIGWIG_DIR = third_party/libBigWig
LIBBIGWIG = $(BUILD_DIR)/libBigWig.a
LIBBIGWIG_SRCS = io.c bwValues.c bwRead.c bwStats.c bwWrite.c
LIBBIGWIG_OBJS = $(addprefix $(BUILD_DIR)/libbigwig/,$(LIBBIGWIG_SRCS:.c=.o))
LIB_SRCS += src/pebble_bigwig.c
ifeq ($(shell uname -s),Linux)
ZLIB_A := $(firstword $(wildcard /usr/lib/x86_64-linux-gnu/libz.a /usr/lib/libz.a))
ifneq ($(ZLIB_A),)
LDLIBS += $(ZLIB_A)
else
$(warning libz.a not found; linking zlib dynamically — omitting -static)
LDFLAGS := $(filter-out -static,$(LDFLAGS))
LDLIBS += -lz
endif
else
LDLIBS += -lz
endif
else
LIBBIGWIG =
LIBBIGWIG_DIR =
LIB_SRCS += src/pebble_bigwig_stub.c
endif

LIB_OBJS  = $(LIB_SRCS:src/%.c=$(BUILD_DIR)/%.o)

.PHONY: all clean test demo example libbigwig

all: $(CLI) $(TESTS) $(IO_TESTS)

ifeq ($(PEBBLE_BIGWIG),1)
all: libbigwig
libbigwig: $(LIBBIGWIG)
endif

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

ifeq ($(PEBBLE_BIGWIG),1)
$(BUILD_DIR)/libbigwig:
	mkdir -p $(BUILD_DIR)/libbigwig

$(BUILD_DIR)/libbigwig/%.o: $(LIBBIGWIG_DIR)/%.c | $(BUILD_DIR)/libbigwig
	$(CC) $(CFLAGS) -I$(LIBBIGWIG_DIR) -Wno-sign-compare -c $< -o $@

$(LIBBIGWIG): $(LIBBIGWIG_OBJS)
	ar rcs $@ $^
endif

$(BUILD_DIR)/%.o: src/%.c | $(BUILD_DIR)
ifeq ($(PEBBLE_BIGWIG),1)
	$(CC) $(CFLAGS) -I$(LIBBIGWIG_DIR) -c $< -o $@
else
	$(CC) $(CFLAGS) -c $< -o $@
endif

$(LIB): $(LIB_OBJS)
	ar rcs $@ $^

$(CLI): $(LIB) $(LIBBIGWIG) src/main.c
ifeq ($(PEBBLE_BIGWIG),1)
	$(CC) $(CFLAGS) -I$(LIBBIGWIG_DIR) src/main.c -o $@ $(LDFLAGS) $(LIB) $(LIBBIGWIG) $(LDLIBS)
else
	$(CC) $(CFLAGS) src/main.c -o $@ $(LDFLAGS) $(LIB) $(LDLIBS)
endif

$(TESTS): tests/test_pebble.c $(LIB) $(LIBBIGWIG)
ifeq ($(PEBBLE_BIGWIG),1)
	$(CC) $(CFLAGS) -I$(LIBBIGWIG_DIR) tests/test_pebble.c -o $@ $(LDFLAGS) $(LIB) $(LIBBIGWIG) $(LDLIBS)
else
	$(CC) $(CFLAGS) tests/test_pebble.c -o $@ $(LDFLAGS) $(LIB) $(LDLIBS)
endif

$(IO_TESTS): tests/test_io.c $(LIB) $(LIBBIGWIG)
ifeq ($(PEBBLE_BIGWIG),1)
	$(CC) $(CFLAGS) -I$(LIBBIGWIG_DIR) tests/test_io.c -o $@ $(LDFLAGS) $(LIB) $(LIBBIGWIG) $(LDLIBS)
else
	$(CC) $(CFLAGS) tests/test_io.c -o $@ $(LDFLAGS) $(LIB) $(LDLIBS)
endif

test: $(TESTS) $(IO_TESTS)
	$(TESTS)
	$(IO_TESTS)

demo: $(CLI)
	$(CLI) --demo

example: $(CLI)
	$(CLI) -i examples/mock_scaffold.bedgraph -o $(BUILD_DIR)/mock_scaffold.smoothed.bedgraph
	cat $(BUILD_DIR)/mock_scaffold.smoothed.bedgraph

clean:
	rm -rf $(BUILD_DIR)
