
/*******************************************************************************
 *
 * File decorrelation_time.c
 *
 * Stampa su file la funzione di autocorrelazione Gamma(t,t_M). Crea un nuovo
 * file per ogni t fino a T_P_MAX, e stampa Gamma al variare di t_M fino a
 * T_M_MAX. Effettua poi il binnaggio e ricalcola Gamma come sopra.
 * Per il calcolo di Gamma effettua N_SWEEP sweep. Per il binnaggio utilizza
 * bin di grandezza D_BIN.
 * Stampa periodicamente a schermo il numero di sweep già effettuati, per
 * monitorare l'avanzamento del programma.
 *
 * Author: Lorenzo Tasca
 *
 *******************************************************************************/

#define MAIN_PROGRAM
#define DECORRELATION_TIME_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "harmonic_oscillator.h"
#include "random.h"

int main(int argc, char *argv[])
{
    int t, t_M;
    double **gamma, **c, **c_binned, **gamma_binned;
    FILE *fd;
    char file_name[50]; /*nome del file in cui viene stampato l'output*/

    if (N_SWEEP != D_BIN * N_BIN)
    {
        printf("ERRORE: impostare N_BIN e D_BIN tali che N_BIN*D_BIN=N_SWEEP.\n");
        exit(EXIT_FAILURE);
    }

    termalisation();

    c = eval_c(10112000);     /*calcolo del correlatore*/
    gamma = eval_gamma(c, 0); /*calcolo di Gamma non binnata*/

    for (t = 1; t < T_P_MAX; t++) /*stampa su file*/
    {
        sprintf(file_name, "data_files/gamma/gamma_t%d_N%d_prova.dat", t, N);
        fd = fopen(file_name, "w");
        for (t_M = 0; t_M < T_M_MAX; t_M++)
            fprintf(fd, "%d\t%f\n", t_M, gamma[t][t_M] / gamma[t][0]);
        fclose(fd);
    }

    c_binned = bin_c(c);                    /*binnaggio dei correlatori*/
    gamma_binned = eval_gamma(c_binned, 1); /*calcolo di Gamma binnata*/

    for (t = 1; t < T_P_MAX; t++) /*stampa su file*/
    {
        sprintf(file_name, "data_files/gamma/gamma_binned_t%d_N%d_prova.dat", t, N);
        fd = fopen(file_name, "w");
        for (t_M = 0; t_M < T_M_MAX; t_M++)
            fprintf(fd, "%d\t%f\n", t_M, gamma_binned[t][t_M] / gamma_binned[t][0]);
        fclose(fd);
    }

    return 0;
}
