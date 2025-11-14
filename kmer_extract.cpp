//
//  kmer_extract.cpp
//  
//
//  Created by Anaelle Ji-Seun Joo on 14/11/2025.
//

#include "kmer_extract.hpp"
#include <vector>

using namespace std;

//trier
std::vector<string> trier(std::vector<string> unsorted_list){
}


//kmerextract
//Entrée :
//k un entier
//F un ensemble de mots avec un alphabet A de taille n
//Sortie : une liste triée de k-mers issus des mots de F
//Variables :
//L une liste de mots sur l’alphabet A
//i,j un entier
//
//début
//L <- liste vide
//pour i de 0 a n-1 faire :
//    pour j de 0 a |F[i]|-k faire :
//        ajouter (F[i][j:j+k-1],L)
//    fin pour
//fin pour
//trier(L)
//retourner L
//fin
