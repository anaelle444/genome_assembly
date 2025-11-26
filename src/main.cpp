////

//  main.cpp//  main.cpp

//  Assembleur de génome basé sur le graphe de De Bruijn//  

////

//  Created by Anaelle Ji-Seun Joo on 14/11/2025.//  Created by Anaelle Ji-Seun Joo on 14/11/2025.

////



#include "kmer_extract.hpp"
#include "calcul_arcs.hpp"
#include "graphe_bruijn.hpp"
#include "chemin_eulerien.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <sys/resource.h>
#include <sys/stat.h>
#include <iomanip>

using namespace std;

// Fonction pour créer un dossier s'il n'existe pas
void creerDossier(const string& chemin) {
    struct stat info;
    if (stat(chemin.c_str(), &info) != 0) {
        // Le dossier n'existe pas, le créer
        mkdir(chemin.c_str(), 0755);
    }
}

// Fonction pour obtenir la mémoire utilisée en MB
double getMemoryUsage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    // ru_maxrss est en kilobytes sur Linux, convertir en MB
    return usage.ru_maxrss / 1024.0;
}

// Fonction pour formater le temps en heures:minutes:secondes
string formatTime(double seconds) {
    int hours = (int)(seconds / 3600);
    int minutes = (int)((seconds - hours * 3600) / 60);
    double secs = seconds - hours * 3600 - minutes * 60;
    
    if (hours > 0) {
        char buffer[50];
        sprintf(buffer, "%dh %dm %.2fs", hours, minutes, secs);
        return string(buffer);
    } else if (minutes > 0) {
        char buffer[50];
        sprintf(buffer, "%dm %.2fs", minutes, secs);
        return string(buffer);
    } else {
        char buffer[50];
        sprintf(buffer, "%.2fs", secs);
        return string(buffer);
    }
}

// Fonction pour lire les séquences depuis un fichier FASTQ
vector<string> lireFastq(const string& nomFichier) {
    vector<string> sequences;
    ifstream fichier(nomFichier);
    
    if (!fichier.is_open()) {
        cerr << "🙈 Erreur : impossible d'ouvrir le fichier " << nomFichier << endl;
        return sequences;
    }
    
    string ligne;
    int numLigne = 0;
    
    while (getline(fichier, ligne)) {
        numLigne++;
        // Dans un fichier FASTQ, la séquence est sur la ligne 2 de chaque bloc de 4 lignes
        if (numLigne % 4 == 2) {
            sequences.push_back(ligne);
        }
    }
    
    fichier.close();
    return sequences;
}

// Fonction pour lire les séquences depuis un fichier FASTA
vector<string> lireFasta(const string& nomFichier) {
    vector<string> sequences;
    ifstream fichier(nomFichier);
    
    if (!fichier.is_open()) {
        cerr << "🙈 Erreur : impossible d'ouvrir le fichier " << nomFichier << endl;
        return sequences;
    }
    
    string ligne;
    string sequenceCourante = "";
    
    while (getline(fichier, ligne)) {
        if (ligne.empty()) continue;
        
        if (ligne[0] == '>') {
            // Nouvelle séquence
            if (!sequenceCourante.empty()) {
                sequences.push_back(sequenceCourante);
                sequenceCourante = "";
            }
        } else {
            // Continuation de la séquence
            sequenceCourante += ligne;
        }
    }
    
    // Ajouter la dernière séquence
    if (!sequenceCourante.empty()) {
        sequences.push_back(sequenceCourante);
    }
    
    fichier.close();
    return sequences;
}

// Fonction pour écrire la séquence assemblée dans un fichier FASTA
void ecrireFasta(const string& nomFichier, const string& sequence, const string& nom = "sequence_assemblee") {
    ofstream fichier(nomFichier);
    
    if (!fichier.is_open()) {
        cerr << "🙈 Erreur : impossible d'écrire dans le fichier " << nomFichier << endl;
        return;
    }
    
    fichier << ">" << nom << endl;
    
    // Écrire la séquence avec 80 caractères par ligne
    for (size_t i = 0; i < sequence.length(); i += 80) {
        fichier << sequence.substr(i, 80) << endl;
    }
    
    fichier.close();
    cout << "🎀 Séquence assemblée écrite dans " << nomFichier << endl;
}

int main(int argc, char* argv[]) {
    // Démarrage du chronomètre
    clock_t tempsDebut = clock();
    double memoireDebut = getMemoryUsage();
    
    cout << "=== Assembleur de génome - Graphe de De Bruijn ===" << endl << endl;
    
    // Paramètres par défaut
    string fichierEntree = "reads.fastq.fq";
    string fichierSortie = "out.fa";
    int k = 21;  // Taille des k-mers
    
    // Lecture des arguments
    if (argc > 1) {
        fichierEntree = argv[1];
    }
    if (argc > 2) {
        k = stoi(argv[2]);
    }
    if (argc > 3) {
        fichierSortie = argv[3];
    }
    
    cout << "Paramètres :" << endl;
    cout << "  Fichier d'entrée : " << fichierEntree << endl;
    cout << "  Taille des k-mers (k) : " << k << endl;
    cout << "  Fichier de sortie : " << fichierSortie << endl << endl;
    
    // Créer le dossier de résultats
    string dossierResultats = "resultats";
    creerDossier(dossierResultats);
    cout << "📁 Dossier de résultats : " << dossierResultats << "/" << endl << endl;
    
    // Étape 1 : Lecture des séquences
    cout << "Étape 1 : Lecture des séquences..." << endl;
    clock_t temps1 = clock();
    vector<string> sequences;
    
    // Déterminer le type de fichier
    if (fichierEntree.find(".fastq") != string::npos || 
        fichierEntree.find(".fq") != string::npos) {
        sequences = lireFastq(fichierEntree);
    } else {
        sequences = lireFasta(fichierEntree);
    }
    
    if (sequences.empty()) {
        cerr << "🙈 Erreur : aucune séquence lue" << endl;
        return 1;
    }
    
    cout << "  " << sequences.size() << " séquences lues" << endl;
    double temps1Ecoule = (double)(clock() - temps1) / CLOCKS_PER_SEC;
    cout << "  ⏱️  Temps : " << formatTime(temps1Ecoule) << endl << endl;
    
    // Étape 2 : Extraction des k-mers
    cout << "Étape 2 : Extraction des k-mers..." << endl;
    clock_t temps2 = clock();
    vector<string> kmers = kmerExtract(k, sequences);
    cout << "  " << kmers.size() << " k-mers extraits et triés" << endl;
    double temps2Ecoule = (double)(clock() - temps2) / CLOCKS_PER_SEC;
    cout << "  ⏱️  Temps : " << formatTime(temps2Ecoule) << endl;
    
    // Écriture des k-mers dans un fichier intermédiaire FASTA
    string fichierKmers = dossierResultats + "/kmers_sorted.fasta";
    cout << "  Écriture des k-mers triés dans " << fichierKmers << "..." << endl;
    ofstream fichierK(fichierKmers);
    if (fichierK.is_open()) {
        for (size_t i = 0; i < kmers.size(); i++) {
            fichierK << ">kmer_" << (i + 1) << endl;
            fichierK << kmers[i] << endl;
        }
        fichierK.close();
        cout << "  ✅ K-mers sauvegardés dans " << fichierKmers << endl;
    } else {
        cerr << "  ⚠️  Avertissement : impossible d'écrire le fichier " << fichierKmers << endl;
    }
    
    // Écriture des k-mers dans un fichier TSV avec index
    string fichierKmersTSV = dossierResultats + "/kmers_sorted.tsv";
    ofstream fichierKTSV(fichierKmersTSV);
    if (fichierKTSV.is_open()) {
        fichierKTSV << "Index\tKmer\tPrefixe\tSuffixe" << endl;
        for (size_t i = 0; i < kmers.size(); i++) {
            string prefixe = kmers[i].substr(0, k-1);
            string suffixe = kmers[i].substr(1, k-1);
            fichierKTSV << i << "\t" << kmers[i] << "\t" << prefixe << "\t" << suffixe << endl;
        }
        fichierKTSV.close();
        cout << "  ✅ K-mers avec index sauvegardés dans " << fichierKmersTSV << endl;
    }
    cout << endl;
    
    // Étape 3 : Calcul des arcs
    cout << "Étape 3 : Calcul des arcs du graphe..." << endl;
    clock_t temps3 = clock();
    vector<pair<int, int>> arcs = calculArcs(kmers, k);
    cout << "  " << arcs.size() << " arcs calculés" << endl;
    double temps3Ecoule = (double)(clock() - temps3) / CLOCKS_PER_SEC;
    cout << "  ⏱️  Temps : " << formatTime(temps3Ecoule) << endl;
    
    // Écriture des arcs dans un fichier TSV
    string fichierArcs = dossierResultats + "/arcs.tsv";
    cout << "  Écriture des arcs dans " << fichierArcs << "..." << endl;
    ofstream fichierA(fichierArcs);
    if (fichierA.is_open()) {
        fichierA << "Source_Index\tDestination_Index\tSource_Kmer\tDestination_Kmer\tChevauchement" << endl;
        for (const auto& arc : arcs) {
            int source = arc.first;
            int dest = arc.second;
            string chevauchement = kmers[source].substr(1, k-1); // suffixe de source = préfixe de dest
            fichierA << source << "\t" << dest << "\t" 
                     << kmers[source] << "\t" << kmers[dest] << "\t"
                     << chevauchement << endl;
        }
        fichierA.close();
        cout << "  ✅ Arcs sauvegardés dans " << fichierArcs << endl;
    } else {
        cerr << "  ⚠️  Avertissement : impossible d'écrire le fichier " << fichierArcs << endl;
    }
    cout << endl;
    
    // Étape 4 : Construction du graphe de De Bruijn
    cout << "Étape 4 : Construction du graphe de De Bruijn..." << endl;
    clock_t temps4 = clock();
    GrapheBruijn graphe = grapheBruijn(kmers, arcs);
    cout << "  Graphe construit avec " << graphe.nombreNoeuds() << " nœuds" << endl;
    double temps4Ecoule = (double)(clock() - temps4) / CLOCKS_PER_SEC;
    cout << "  ⏱️  Temps : " << formatTime(temps4Ecoule) << endl;
    
    // Écriture du graphe de De Bruijn dans un fichier TXT (format lisible)
    string fichierGrapheTXT = dossierResultats + "/graphe_debruijn.txt";
    cout << "  Écriture du graphe dans " << fichierGrapheTXT << "..." << endl;
    ofstream fichierGTXT(fichierGrapheTXT);
    if (fichierGTXT.is_open()) {
        fichierGTXT << "=== Graphe de De Bruijn ===" << endl;
        fichierGTXT << "Nombre de nœuds : " << graphe.nombreNoeuds() << endl;
        fichierGTXT << "Nombre d'arcs : " << arcs.size() << endl;
        fichierGTXT << "Taille des k-mers : " << k << endl << endl;
        
        fichierGTXT << "=== Liste des nœuds et leurs successeurs ===" << endl << endl;
        
        const vector<Noeud>& noeuds = graphe.getNoeuds();
        for (size_t i = 0; i < noeuds.size(); i++) {
            fichierGTXT << "Nœud " << i << " : " << noeuds[i].kmer << endl;
            fichierGTXT << "  Successeurs (" << noeuds[i].successeurs.size() << ") : ";
            if (noeuds[i].successeurs.empty()) {
                fichierGTXT << "aucun";
            } else {
                for (size_t j = 0; j < noeuds[i].successeurs.size(); j++) {
                    int succ = noeuds[i].successeurs[j];
                    fichierGTXT << succ << " (" << noeuds[succ].kmer << ")";
                    if (j < noeuds[i].successeurs.size() - 1) {
                        fichierGTXT << ", ";
                    }
                }
            }
            fichierGTXT << endl << endl;
        }
        
        fichierGTXT.close();
        cout << "  ✅ Graphe sauvegardé dans " << fichierGrapheTXT << endl;
    } else {
        cerr << "  ⚠️  Avertissement : impossible d'écrire le fichier " << fichierGrapheTXT << endl;
    }
    
    // Écriture du graphe au format DOT pour visualisation avec Graphviz
    string fichierGrapheDOT = dossierResultats + "/graphe_debruijn.dot";
    cout << "  Écriture du graphe au format DOT dans " << fichierGrapheDOT << "..." << endl;
    ofstream fichierGDOT(fichierGrapheDOT);
    if (fichierGDOT.is_open()) {
        fichierGDOT << "digraph DeBruijnGraph {" << endl;
        fichierGDOT << "  rankdir=LR;" << endl;
        fichierGDOT << "  node [shape=circle, fontsize=10];" << endl;
        fichierGDOT << "  edge [fontsize=8];" << endl << endl;
        
        // Ajouter les nœuds
        const vector<Noeud>& noeuds = graphe.getNoeuds();
        for (size_t i = 0; i < noeuds.size(); i++) {
            fichierGDOT << "  " << i << " [label=\"" << i << "\\n" << noeuds[i].kmer << "\"];" << endl;
        }
        fichierGDOT << endl;
        
        // Ajouter les arcs
        for (const auto& arc : arcs) {
            fichierGDOT << "  " << arc.first << " -> " << arc.second << ";" << endl;
        }
        
        fichierGDOT << "}" << endl;
        fichierGDOT.close();
        cout << "  ✅ Graphe DOT sauvegardé dans " << fichierGrapheDOT << endl;
        cout << "  💡 Pour visualiser : dot -Tpng " << fichierGrapheDOT << " -o graphe.png" << endl;
    } else {
        cerr << "  ⚠️  Avertissement : impossible d'écrire le fichier " << fichierGrapheDOT << endl;
    }
    cout << endl;
    
    // Étape 5 : Recherche du chemin eulérien et assemblage
    cout << "Étape 5 : Recherche du chemin eulérien et assemblage..." << endl;
    clock_t temps5 = clock();
    string sequenceAssemblee = cheminEulerienEtAssemblage(graphe, kmers, k);
    cout << "  Séquence assemblée : " << sequenceAssemblee.length() << " bases" << endl;
    double temps5Ecoule = (double)(clock() - temps5) / CLOCKS_PER_SEC;
    cout << "  ⏱️  Temps : " << formatTime(temps5Ecoule) << endl;
    
    // Écriture du chemin eulérien dans un fichier
    string fichierChemin = dossierResultats + "/chemin_eulerien.txt";
    cout << "  Écriture du chemin eulérien dans " << fichierChemin << "..." << endl;
    ofstream fichierC(fichierChemin);
    if (fichierC.is_open()) {
        fichierC << "=== Chemin Eulérien ===" << endl;
        fichierC << "Longueur de la séquence assemblée : " << sequenceAssemblee.length() << " bases" << endl << endl;
        fichierC << "Séquence assemblée :" << endl;
        fichierC << sequenceAssemblee << endl;
        fichierC.close();
        cout << "  ✅ Chemin eulérien sauvegardé dans " << fichierChemin << endl;
    }
    cout << endl;
    
    // Étape 6 : Écriture du résultat
    cout << "Étape 6 : Écriture du résultat..." << endl;
    // Ajouter le chemin du dossier résultats si le fichier de sortie n'a pas de chemin
    string fichierSortieFinal = fichierSortie;
    if (fichierSortie.find('/') == string::npos) {
        fichierSortieFinal = dossierResultats + "/" + fichierSortie;
    }
    ecrireFasta(fichierSortieFinal, sequenceAssemblee);
    
    // Calcul des statistiques finales
    clock_t tempsFin = clock();
    double tempsTotal = (double)(tempsFin - tempsDebut) / CLOCKS_PER_SEC;
    double memoireFin = getMemoryUsage();
    double memoireUtilisee = memoireFin - memoireDebut;
    
    cout << endl << "=== Assemblage terminé avec succès ===" << endl;
    cout << endl << "📊 STATISTIQUES D'EXÉCUTION" << endl;
    cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" << endl;
    cout << "⏱️  Temps total d'exécution : " << formatTime(tempsTotal) << endl;
    cout << "💾 Mémoire utilisée : " << fixed << setprecision(2) << memoireUtilisee << " MB" << endl;
    cout << "💾 Mémoire maximale : " << fixed << setprecision(2) << memoireFin << " MB" << endl;
    cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" << endl;
    
    cout << endl << "📁 FICHIERS GÉNÉRÉS" << endl;
    cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" << endl;
    cout << "  ✅ " << dossierResultats << "/kmers_sorted.fasta - K-mers triés (FASTA)" << endl;
    cout << "  ✅ " << dossierResultats << "/kmers_sorted.tsv - K-mers avec index (TSV)" << endl;
    cout << "  ✅ " << dossierResultats << "/arcs.tsv - Liste des arcs du graphe (TSV)" << endl;
    cout << "  ✅ " << dossierResultats << "/graphe_debruijn.txt - Description du graphe" << endl;
    cout << "  ✅ " << dossierResultats << "/graphe_debruijn.dot - Graphe pour visualisation" << endl;
    cout << "  ✅ " << dossierResultats << "/chemin_eulerien.txt - Chemin eulérien et séquence" << endl;
    cout << "  ✅ " << fichierSortieFinal << " - Séquence assemblée finale (FASTA)" << endl;
    cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" << endl;
    
    cout << endl << "💡 VISUALISATION DU GRAPHE" << endl;
    cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" << endl;
    cout << "Pour générer une image du graphe, utilisez Graphviz :" << endl;
    cout << "  $ dot -Tpng " << dossierResultats << "/graphe_debruijn.dot -o " << dossierResultats << "/graphe.png" << endl;
    cout << "  $ dot -Tsvg " << dossierResultats << "/graphe_debruijn.dot -o " << dossierResultats << "/graphe.svg" << endl;
    cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" << endl;
    
    cout << endl << "✨ Thank you for trusting us with your genome assembly ✨" << endl;
    
    return 0;
}
