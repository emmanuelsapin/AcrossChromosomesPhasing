# Makefile pour le Générateur de Données Aléatoires

CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall
TARGET = GenerateurDonneesAleatoires
SOURCE = GenerateurDonneesAleatoires.cpp

# Cible par défaut
all: $(TARGET)

# Compilation
$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCE)

# Nettoyage
clean:
	rm -f $(TARGET) $(TARGET).exe

# Installation (copie dans un répertoire système, optionnel)
install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

# Aide
help:
	@echo "Makefile pour GenerateurDonneesAleatoires"
	@echo ""
	@echo "Cibles disponibles :"
	@echo "  make          - Compile le générateur"
	@echo "  make clean    - Supprime les fichiers compilés"
	@echo "  make install  - Installe dans /usr/local/bin (nécessite sudo)"
	@echo "  make help     - Affiche cette aide"

.PHONY: all clean install help
