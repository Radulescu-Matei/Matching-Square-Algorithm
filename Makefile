build: matching_square.c
	gcc matching_square.c helpers.c -o matching_square -lm -lpthread -Wall -Wextra
clean:
	rm -rf matching_square