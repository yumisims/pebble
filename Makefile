CC      = cc
CFLAGS  = -std=c11 -Wall -Wextra -Wpedantic -Iinclude
LDFLAGS =
LDLIBS  = -lm

BUILD_DIR = build
LIB       = $(BUILD_DIR)/libpebble.a
CLI       = $(BUILD_DIR)/pebble
TESTS     = $(BUILD_DIR)/pebble_tests
IO_TESTS  = $(BUILD_DIR)/pebble_io_tests

LIB_SRCS  = src/pebble.c src/pebble_io.c
LIB_OBJS  = $(LIB_SRCS:src/%.c=$(BUILD_DIR)/%.o)

.PHONY: all clean test demo example

all: $(CLI) $(TESTS) $(IO_TESTS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: src/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(LIB): $(LIB_OBJS)
	ar rcs $@ $^

$(CLI): $(LIB) src/main.c
	$(CC) $(CFLAGS) src/main.c -o $@ $(LDFLAGS) $(LIB) $(LDLIBS)

$(TESTS): tests/test_pebble.c $(LIB)
	$(CC) $(CFLAGS) tests/test_pebble.c -o $@ $(LDFLAGS) $(LIB) $(LDLIBS)

$(IO_TESTS): tests/test_io.c $(LIB)
	$(CC) $(CFLAGS) tests/test_io.c -o $@ $(LDFLAGS) $(LIB) $(LDLIBS)

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
