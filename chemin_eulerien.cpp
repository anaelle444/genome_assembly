//
//  chemin_eulerien.cpp
//  
//
//  Created by Anaelle Ji-Seun Joo on 20/11/2025.
//

#include "chemin_eulerien.hpp"
#include "graphe_bruijn.hpp"
#include <stack>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

// CheminEulerienEtAssemblage
// Entrée : T : graphe orienté de Bruijn eulérien
//          L : liste des k-mers
//          k : entier (taille des k-mers)
// Sortie : S : séquence assemblée
string cheminEulerienEtAssemblage(GrapheBruijn& T, 
                                  const vector<string>& L, 
                                  int k) {
    stack<int> pile;
    vector<int> P;  // Chemin eulérien
    
    // Commencer avec le premier nœud du graphe (nœud 0)
    int v = 0;
    pile.push(v);
    
    // Algorithme de Hierholzer pour trouver un chemin eulérien
    while (!pile.empty()) {
        v = pile.top();
        
        if (T.possèdeArcSortant(v)) {
            // Si le nœud a un arc sortant, suivre cet arc
            int u = T.retirerUnSuccesseur(v);
            pile.push(u);
        } else {
            // Sinon, ajouter le nœud au chemin et le dépiler
            pile.pop();
            P.push_back(v);
        }
    }
    
    // Inverser le chemin (car il a été construit à l'envers)
    reverse(P.begin(), P.end());
    
    // Reconstruction de la séquence à partir du chemin
    if (P.empty()) {
        return "";
    }
    
    // Commencer avec le premier k-mer complet
    string S = L[P[0]];
    
    // Ajouter le dernier caractère de chaque k-mer suivant
    for (size_t i = 1; i < P.size(); i++) {
        S += L[P[i]][k - 1];
    }
    
    return S;
}
