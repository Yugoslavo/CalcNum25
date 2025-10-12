#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "csr_io.h"

/*J'ai crée cette fonction pour voir a quoi la matrice en CSR ressemble
Et puis si j'ai le temps je vais créer une autre fonction pour comparer
*/
int write_csr_arrays(const char *prefix,int n,const int *ia, const int *ja, const double *a){
    char f_ia[256], f_ja[256],f_a[256];
    snprintf(f_ia,sizeof f_ia, "%s_ia.txt",prefix);
    snprintf(f_ja,sizeof f_ja, "%s_ja.txt",prefix);
    snprintf(f_a,sizeof f_a, "%s_a.txt", prefix);

    FILE *fia = fopen(f_ia,"w"); 
    FILE *fja = fopen(f_ja, "w");
    FILE *fa  = fopen(f_a , "w");
    if (!fia || !fja || !fa){
        if (fia) fclose(fia); if (fja) fclose(fja); if (fa) fclose(fa);
        fprintf(stderr, "Erreur: impossible d'ouvrir les fichiers de sortie.\n");
        return 0;
    }

  // ia: n+1 entiers
  for (int i=0; i<=n; ++i) fprintf(fia, "%d\n", ia[i]);

  // nnz = ia[n]
  int nnz = ia[n];

  // ja: nnz entiers
  for (int k=0; k<nnz; ++k) fprintf(fja, "%d\n", ja[k]);

  // a: nnz doubles (format scientifique pour sécurité)
  for (int k=0; k<nnz; ++k) fprintf(fa, "%.15f\n", a[k]);

  fclose(fia); fclose(fja); fclose(fa);
  return 1;
}