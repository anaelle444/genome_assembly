# Compilateur et options
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

# Nom de l'exécutable
OUT = assembler

# Fichiers sources (.cpp)
SOURCES = main.cpp \
          kmer_extract.cpp \
          calcul_arcs.cpp \
          graphe_bruijn.cpp \
          chemin_eulerien.cpp

# Fichiers objets (.o) générés à partir des sources
OBJECTS = $(SOURCES:.cpp=.o)

# Fichiers headers (.hpp)
HEADERS = kmer_extract.hpp \
          calcul_arcs.hpp \
          graphe_bruijn.hpp \
          chemin_eulerien.hpp

# Règle par défaut : compile l'exécutable
all: $(OUT)

# Règle pour créer l'exécutable à partir des objets
$(OUT): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(OUT) $(OBJECTS)

# Règle pour compiler chaque fichier .cpp en .o
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Règle pour nettoyer les fichiers compilés
clean:
	rm -f $(OUT) $(OBJECTS)

# Règle pour tout recompiler depuis zéro
rebuild: clean all

# Déclaration des règles qui ne sont pas des fichiers
.PHONY: all clean rebuild
