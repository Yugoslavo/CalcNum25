Guide de lecture du projet et ligne directrice

Ce document accompagne mon code et explique, question par question, les choix effectués ainsi que la logique générale du projet. Il complète la présentation orale.

QUESTION 1 – Génération de la matrice A (prob.c)

La majeure partie du travail pour cette question se trouve dans le fichier prob.c, où j’ai adapté le code pour gérer la géométrie fournie dans l’énoncé.

J’ai choisi d’utiliser un système de mapping (similaire à un dictionnaire en Python) afin d’associer chaque point intérieur de la grille à son indice lexicographique. Cette approche rend le code plus clair, flexible et évite les cas particuliers.

L’objectif était également de conserver une flexibilité totale : si les paramètres de la fonction in_hole() sont modifiés pour définir une autre géométrie ou un autre niveau de raffinement, la matrice CSR est automatiquement recalculée sans autre changement.

Pour faciliter le débogage, j’ai ajouté la fonction write_csr_arrays dans debug.c, permettant de visualiser la structure CSR de la matrice générée.

QUESTION 2 – Calcul du résidu relatif

Cette question a été traitée dans main.c via la fonction residual_ratio.

La fonction calcule efficacement le résidu relatif ||AΦ − β²Φ||₂ / ||Φ||₂ en optimisant les accès mémoire dans la structure CSR, ce qui est essentiel pour respecter les critères d’efficacité demandés dans l’énoncé.

L’implémentation reste concise et intuitive.

QUESTION 3 – Calcul du facteur critique alpha_crit et homothétie

Pour cette partie, j’utilise le solveur PRIMME pour calculer la plus petite valeur propre de la matrice A associée à la géométrie de base.

Comme demandé dans l’énoncé, j’ai introduit un facteur d’homothétie alpha, appliqué directement dans la génération de la matrice :

Au début, alpha = 1 pour déterminer la valeur propre β² associée à ma géométrie brute.

J’en déduis ensuite la valeur alpha_critique, c’est-à-dire le facteur d’homothétie permettant d’obtenir un réacteur stationnaire.

Toutes mes fonctions (notamment dans prob.c) acceptent ce facteur alpha comme paramètre, ce qui me permet de recalculer facilement la matrice pour n’importe quelle géométrie homothétique.

Grâce à cette approche, j’ai pu vérifier que pour alpha = alpha_crit, la plus petite valeur propre vaut bien 4, avec un résidu relatif inférieur à 1e-9, confirmant l’exactitude du calcul.

QUESTION 4 – Visualisation graphique du flux stationnaire

Pour cette question, j’ai utilisé gnuplot via des pipes ("gp") comme suggéré dans l’énoncé.

La difficulté principale était de rendre le flux phi compatible avec le format attendu par gnuplot. Ce traitement est encapsulé dans plotflux.c, où le flux est reconstruit sur la grille complète, y compris au bord, comme demandé.

J’ai vérifié que le mode propre fondamental respecte bien les conditions de Dirichlet.

QUESTION 5 – Résolution temporelle par Euler progressive

La fonction dans eulerprog.c implémente la méthode d’Euler progressive pour résoudre le problème de Cauchy discrétisé.

Le pas de temps et le nombre d’itérations sont choisis de manière à rendre l’évolution du flux clairement visible, notamment lorsque l’amplitude du flux double.

J’ai représenté :

un cas surcritique (alpha = alpha_crit × 1.01),

un cas souscritique (alpha = alpha_crit × 0.99).

J’ai choisi de tracer ces deux cas statiquement plutôt qu’en animation, pour une présentation plus lisible et moins coûteuse en mémoire.

QUESTION 6 – Solveur alternatif (ARPACK / JADAMILU)

J’ai d’abord implémenté une interface avec ARPACK, qui fonctionne correctement pour ma matrice.

Ensuite, j’ai également intégré JADAMILU, un solveur spécialisé pour matrices creuses symétriques. Finalement, j’ai choisi d’utiliser JADAMILU pour cette question car il exploite mieux la structure creuse, converge rapidement et ne nécessite que la partie triangulaire supérieure de la matrice.

La comparaison avec PRIMME montre une excellente cohérence des résultats.
ARPACK reste pleinement fonctionnel dans mon projet, mais n’est plus utilisé comme solveur principal.