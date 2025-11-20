//
//  graphe_bruijn.hpp
//  
//
//  Created by Anaelle Ji-Seun Joo on 20/11/2025.
//

#ifndef graphe_bruijn_hpp
#define graphe_bruijn_hpp

#include <string>
#include <vector>
#include <utility>  // Pour std::pair

// Structure représentant un nœud du graphe
struct Noeud {
    int id;
    std::string kmer;
    std::vector<int> successeurs;  // Liste des indices des nœuds successeurs
    
    Noeud() : id(-1), kmer("") {}
    Noeud(int i, const std::string& k) : id(i), kmer(k) {}
};

// Classe représentant un graphe de De Bruijn
class GrapheBruijn {
private:
    std::vector<Noeud> noeuds;
    
public:
    GrapheBruijn();
    GrapheBruijn(int taille);
    
    // Crée un nœud avec un identifiant et un k-mer
    void creerNoeud(int id, const std::string& kmer);
    
    // Ajoute un arc du nœud source vers le nœud destination
    void ajouterArc(int source, int destination);
    
    // Vérifie si un nœud possède au moins un arc sortant
    bool possèdeArcSortant(int noeud) const;
    
    // Retire et retourne un successeur du nœud (pour algorithme eulérien)
    int retirerUnSuccesseur(int noeud);
    
    // Accesseurs
    int nombreNoeuds() const { return noeuds.size(); }
    const Noeud& getNoeud(int id) const { return noeuds[id]; }
    const std::vector<Noeud>& getNoeuds() const { return noeuds; }
};

// Construit le graphe de De Bruijn
// L : liste des k-mers
// A : liste de paires d'indices de k-mers (arcs)
// Retourne : Graphe orienté de De Bruijn
GrapheBruijn grapheBruijn(const std::vector<std::string>& L, 
                          const std::vector<std::pair<int, int>>& A);

#endif /* graphe_bruijn_hpp */
