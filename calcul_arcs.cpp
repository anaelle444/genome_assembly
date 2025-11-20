//
//  calcul_arcs.cpp
//  
//
//  Created by Anaelle Ji-Seun Joo on 20/11/2025.
//

#include "calcul_arcs.hpp"
#include <algorithm>
#include <string>
#include <vector>
#include <utility>

using namespace std;

// Recherche binaire d'un k-mer dans une liste triée
// Retourne l'indice du k-mer, ou -1 s'il n'est pas trouvé
int rechercher(const string& kmer, const vector<string>& L) {
    auto it = lower_bound(L.begin(), L.end(), kmer);
    
    if (it != L.end() && *it == kmer) {
        return distance(L.begin(), it);
    }
    return -1;
}

// CalculArcs
// Entrée : L : liste de k-mers (ordre alphabétique), k : entier
// Sortie : A : liste des arcs (paires d'indices des k-mers)
vector<pair<int, int>> calculArcs(const vector<string>& L, int k) {
    vector<pair<int, int>> A;  // Liste des arcs
    char bases[] = {'A', 'C', 'G', 'T'};
    
    // Pour chaque k-mer dans L
    for (size_t i = 0; i < L.size(); i++) {
        // Récupérer le suffixe de longueur k-1 du k-mer courant
        string suffixe = L[i].substr(1, k - 1);
        
        // Tester chaque base pour voir si elle forme un k-mer existant
        for (char X : bases) {
            // Construire le k-mer suivant potentiel
            string kmer_suivant = suffixe + X;
            
            // Chercher si ce k-mer existe dans L
            int j = rechercher(kmer_suivant, L);
            
            // Si le k-mer existe, ajouter l'arc
            if (j != -1) {
                A.push_back(make_pair(i, j));
            }
        }
    }
    
    return A;
}
