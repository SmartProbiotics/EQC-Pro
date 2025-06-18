CC = g++
CFLAGS = -O3 -std=c++11
OBJDIR = .obj
SRCS = Graph.cpp
OBJS = $(OBJDIR)/Graph.o

all: EQC EQC-TD EQC-NH

# ========= Standard EQC =========
EQC: main.cpp $(SRCS)
	$(CC) $(CFLAGS) -DNO_PROGRESS_BAR -DNDEBUG -c Graph.cpp -o $(OBJDIR)/Graph.o
	$(CC) $(CFLAGS) -DNO_PROGRESS_BAR -DNDEBUG -c main.cpp -o $(OBJDIR)/main.o
	$(CC) $(OBJDIR)/main.o $(OBJDIR)/Graph.o -o EQC

# ========= Top-down variant =========
EQC-TD: EQC-TD.cpp $(SRCS)
	$(CC) $(CFLAGS) -DNO_PROGRESS_BAR -DNDEBUG -c Graph.cpp -o $(OBJDIR)/Graph.o
	$(CC) $(CFLAGS) -DNO_PROGRESS_BAR -DNDEBUG -c EQC-TD.cpp -o $(OBJDIR)/EQC-TD.o
	$(CC) $(OBJDIR)/EQC-TD.o $(OBJDIR)/Graph.o -o EQC-TD

# ========= Bottom-up no-heuristic variant =========
EQC-NH: main.cpp $(SRCS)
	$(CC) $(CFLAGS) -DNO_PROGRESS_BAR -DNDEBUG -DNHEU -c Graph.cpp -o $(OBJDIR)/Graph.o
	$(CC) $(CFLAGS) -DNO_PROGRESS_BAR -DNDEBUG -DNHEU -c main.cpp -o $(OBJDIR)/main.o
	$(CC) $(OBJDIR)/main.o $(OBJDIR)/Graph.o -o EQC-NH

clean:
	rm -f $(OBJDIR)/*.o EQC EQC-TD EQC-NH

analysis:
	mkdir -p .obj
	g++ -pg -O3 -c -I. -std=c++11 -o .obj/main.o main.cpp
	g++ -pg -O3 -c -I. -std=c++11 -o .obj/Graph.o Graph.cpp
	g++ -pg -O3 .obj/main.o .obj/Graph.o -o EQC 
	rm .obj/*.o
test:
	./EQC ~/Program/DATA/realworld_graphs/soc/soc-LiveMocha.bin 0.99
	# ./EQC ./datasets/bio-celegans.bin 0.8
	gprof ./EQC gmon.out > analysis.txt
debug:
	g++ -g -std=c++11 -o EQC main.cpp Graph.cpp