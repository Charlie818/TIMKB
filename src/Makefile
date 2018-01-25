DIR=../src

CC=g++
CFLAGS=-O3
# CFLAGS=-g3 -Wall
LDFLAGS=-lm 
TAGS=main
target_tests=tic.o tim.o tests.o
target_offline=tic.o tim.o offline.o
target_experiments=tic.o tim.o main.o
target_pmia_no_topic_spread=tic.o tim.o pmia_no_topic_spread.o

all:$(TAGS)

main: $(target_experiments)
		$(CC) $(CFLAGS) $(target_experiments) -o main $(LDFLAGS)
main.o: $(DIR)/main.cpp
		$(CC) $(CFLAGS) -o main.o -c $(DIR)/main.cpp

tests: $(target_tests)
		$(CC) $(CFLAGS) $(target_tests) -o tests $(LDFLAGS)
tests.o: $(DIR)/tests.cpp
		$(CC) $(CFLAGS) -o tests.o -c $(DIR)/tests.cpp

pmia_no_topic_spread: $(target_pmia_no_topic_spread)
		$(CC) $(CFLAGS) $(target_pmia_no_topic_spread) -o pmia_no_topic_spread $(LDFLAGS)
pmia_no_topic_spread.o: $(DIR)/pmia_no_topic_spread.cpp
		$(CC) $(CFLAGS) -o pmia_no_topic_spread.o -c $(DIR)/pmia_no_topic_spread.cpp

offline: $(target_offline)
		$(CC) $(CFLAGS) $(target_offline) -o offline $(LDFLAGS)
offline.o: $(DIR)/offline.cpp
		$(CC) $(CFLAGS) -o offline.o -c $(DIR)/offline.cpp

tim.o: $(DIR)/tim.cpp
		$(CC) $(CFLAGS) -o tim.o -c $(DIR)/tim.cpp
tic.o: $(DIR)/tic.cpp
		$(CC) $(CFLAGS) -o tic.o -c $(DIR)/tic.cpp
clean:
	rm -rf *.o $(TAGS)
