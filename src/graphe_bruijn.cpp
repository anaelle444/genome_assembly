//
//  graphe_bruijn.cpp
//  
//
//  Created by Anaelle Ji-Seun Joo on 20/11/2025.
//

#include "graphe_bruijn.hpp"
#include <vector>
#include <string>
#include <utility>

using namespace std;

// Constructeurs de la classe GrapheBruijn
GrapheBruijn::GrapheBruijn() {}

GrapheBruijn::GrapheBruijn(int taille) {
    noeuds.reserve(taille);
}

// Crée un nœud avec un identifiant et un k-mer
void GrapheBruijn::creerNoeud(int id, const string& kmer) {
    noeuds.push_back(Noeud(id, kmer));
}

// Ajoute un arc du nœud source vers le nœud destination
void GrapheBruijn::ajouterArc(int source, int destination) {
    if (source >= 0 && source < (int)noeuds.size()) {
        noeuds[source].successeurs.push_back(destination);
    }
}

// Vérifie si un nœud possède au moins un arc sortant
bool GrapheBruijn::possèdeArcSortant(int noeud) const {
    if (noeud >= 0 && noeud < (int)noeuds.size()) {
        return !noeuds[noeud].successeurs.empty();
    }
    return false;
}

// Retire et retourne un successeur du nœud
int GrapheBruijn::retirerUnSuccesseur(int noeud) {
    if (noeud >= 0 && noeud < (int)noeuds.size() && !noeuds[noeud].successeurs.empty()) {
        int successeur = noeuds[noeud].successeurs.back();
        noeuds[noeud].successeurs.pop_back();
        return successeur;
    }
    return -1;
}

// GrapheBruijn
// Entrée : L : liste des k-mers de longueur n
//          A : liste de paires d'indices de k-mers de longueur m
// Sortie : T : Graphe orienté
GrapheBruijn grapheBruijn(const vector<string>& L, 
                          const vector<pair<int, int>>& A) {
    int n = L.size();
    
    // Création d'un graphe vide avec n nœuds pré-alloués
    GrapheBruijn T(n);
    
    // Création de tous les nœuds en une seule passe
    for (int i = 0; i < n; i++) {
        T.creerNoeud(i, L[i]);
    }
    
    // Ajout des arcs
    for (const auto& arc : A) {
        T.ajouterArc(arc.first, arc.second);
    }
    
    return T;
}
