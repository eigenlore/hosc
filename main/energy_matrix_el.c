
/*******************************************************************************
 *
 * File energy_matrix_el.c
 *
 * Questo codice effettua il calcolo del gap energetico e dell'elemento di
 * matrice. Stampa su file le grandezze al variare di t, e le grandezze
 * mediate finali, con associato l'errore.
 * Riceve da riga di comando il seed della simulazione.
 * Stampa a schermo la percentuale di sweep già effettuati, per
 * monitorare l'avanzamento del programma.
 * Stampa su file le estrazioni del correlatore binnato.
 *
 * Author: Lorenzo Tasca
 *
 *******************************************************************************/

#define MAIN_PROGRAM
#define ENERGY_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "harmonic_oscillator.h"
#include "random.h"

int main(int argc, char *argv[])
{
    int t, seed, i;
    double **c_binned, *c_bar, **cluster_c, *deltaE, **cluster_deltaE, *error_deltaE, *matrix_el, **cluster_matrix_el, *error_matrix_el, *mean_deltaE, *mean_matrix_el;
    FILE *fd;
    char file_name[100]; /*nome del file.dat in cui viene stampato l'output*/

    if (argc < 2)
    {
        printf("ERRORE: inserire un seed da riga di comando.\n");
        exit(EXIT_FAILURE);
    }

    termalisation();
    seed = atoi(argv[1]); /*passo il seed da riga di comando*/

    /*calcolo di tutte le grandezze fisiche*/
    c_binned = eval_c_binned(seed);
    printf("Analisi dei dati in corso\n");
    c_bar = eval_c_bar(c_binned);
    cluster_c = eval_cluster_c(c_bar, c_binned);
    sprintf(file_name, "data_files/energy/cmedieprova_N%d.dat", N);


    fd = fopen(file_name, "w");
    for (t = 0; t < N; t++)
        fprintf(fd, "%f\n", c_bar[t]);
    fprintf(fd, "\n");

    fclose(fd);

    /*stampa su file di cbinned*/ 
     sprintf(file_name, "data_files/energy/cbinnateprova_N%d.dat", N);
     fd = fopen(file_name, "w");
     for (i = 0; i < N_BIN; i++)
     {
         for (t = 0; t < N; t++)
             fprintf(fd, "%f\t", c_binned[t][i]);
         fprintf(fd, "\n");
     }
     fclose(fd);
    free(c_binned);

    deltaE = eval_deltaE(c_bar);
    matrix_el = eval_matrix_el(c_bar, deltaE);
    free(c_bar);
    cluster_deltaE = eval_cluster_deltaE(cluster_c);
    cluster_matrix_el = eval_cluster_matrix_el(cluster_c, cluster_deltaE);
    free(cluster_c);
    error_deltaE = eval_error_jackknife(deltaE, cluster_deltaE);
    error_matrix_el = eval_error_jackknife(matrix_el, cluster_matrix_el);
    mean_deltaE = eval_mean_x(deltaE, error_deltaE, cluster_deltaE);
    mean_matrix_el = eval_mean_x(matrix_el, error_matrix_el, cluster_matrix_el);

    /*stampa su file di DeltaE al variare di t*/
    sprintf(file_name, "data_files/energy/deltaEprova_N%d.dat", N);
    fd = fopen(file_name, "w");
    for (t = 0; t < N; t++)
        fprintf(fd, "%f\t%f\n", deltaE[t], error_deltaE[t]);
    fclose(fd);

    /*stampa su file elemento di matrice al variare di t*/
    sprintf(file_name, "data_files/energy/matrix_elprova_N%d.dat", N);
    fd = fopen(file_name, "w");
    for (t = 0; t < N; t++)
        fprintf(fd, "%f\t%f\n", matrix_el[t], error_matrix_el[t]);
    fclose(fd);

    /*stampa su file di DeltaE medio*/
    sprintf(file_name, "data_files/energy/mean_deltaEprova_N%d.dat", N);
    fd = fopen(file_name, "w");
    fprintf(fd, "%.15f\t%.15f\n", mean_deltaE[0], mean_deltaE[1]);
    fclose(fd);

    /*stampa su file elemento di matrice medio*/
    sprintf(file_name, "data_files/energy/mean_matrix_elprova_N%d.dat", N);
    fd = fopen(file_name, "w");
    fprintf(fd, "%.15f\t%.15f\n", mean_matrix_el[0], mean_matrix_el[1]);
    fclose(fd);
    
    return 0;
}
