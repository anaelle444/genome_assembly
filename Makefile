# Compilateur et options
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -Iinclude

# Nom de l'exécutable
OUT = assembler

# Dossiers
SRC_DIR = src
INC_DIR = include
OBJ_DIR = obj

# Fichiers sources (.cpp) dans src/
SOURCES = $(SRC_DIR)/main.cpp \
          $(SRC_DIR)/kmer_extract.cpp \
          $(SRC_DIR)/calcul_arcs.cpp \
          $(SRC_DIR)/graphe_bruijn.cpp \
          $(SRC_DIR)/chemin_eulerien.cpp

# Fichiers objets (.o) générés dans obj/
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

# Fichiers headers (.hpp) dans include/
HEADERS = $(INC_DIR)/kmer_extract.hpp \
          $(INC_DIR)/calcul_arcs.hpp \
          $(INC_DIR)/graphe_bruijn.hpp \
          $(INC_DIR)/chemin_eulerien.hpp

# Règle par défaut : compile l'exécutable
all: $(OBJ_DIR) $(OUT)

# Créer le dossier obj/ s'il n'existe pas
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Règle pour créer l'exécutable à partir des objets
$(OUT): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(OUT) $(OBJECTS)

# Règle pour compiler chaque fichier .cpp en .o
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Règle pour nettoyer les fichiers compilés
clean:
	rm -f $(OUT) $(OBJ_DIR)/*.o
	rm -rf $(OBJ_DIR)

# Règle pour tout recompiler depuis zéro
rebuild: clean all

# Déclaration des règles qui ne sont pas des fichiers
.PHONY: all clean rebuild
