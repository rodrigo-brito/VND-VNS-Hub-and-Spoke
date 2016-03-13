all:
		g++ -w src/main.cpp -o bin/main
		@echo "\n--\nCompilação finalizada...\ndigite 'make run' para executar\n--"
run:
		./bin/main
