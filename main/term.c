
/*******************************************************************************
 *
 * File term.c
 *
 * Il programma riceve da riga di comando la configurazione iniziale da cui
 * partire per termalizzare (0 a freddo, oppure 1 a caldo). Successivamente
 * esegue N_TERM sweep, stampa l'accettanza media e stampa su file l'azione
 * a ogni sweep.
 *
 * Author: Lorenzo Tasca
 *
 *******************************************************************************/

#define MAIN_PROGRAM
#define TERM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "harmonic_oscillator.h"
#include "random.h"

int main(int argc, char *argv[])
{
    int i;
    int heat;        /*0 a freddo, 1 a caldo*/
    double acc_rate; /*accettanza*/
    FILE *fd;
    char file_name[100]; /*nome del file in cui viene stampato l'output*/

    if (argc < 2)
    {
        printf("ERRORE: specificare da riga di comando se si vuole partire da caldo o freddo (0 freddo, 1 caldo).\n");
        exit(EXIT_FAILURE);
    }

    heat = atoi(argv[1]);
    acc_rate = 0;

    if (heat == 0)
    {
        for (i = 0; i < N; i++)
            xx[i] = 0;                                                                        /*a freddo tutte le coordinate nulle*/
        sprintf(file_name, "data_files/thermalisation/thermalisation_cold_N%d_prova.dat", N); /*per evitare di sovrascrivere i dati al run del programma c'è l'attributo 'prova'*/
        fd = fopen(file_name, "w");
    }
    else
    {
        rlxd_init(1, 3122000);
        ranlxd(xx, N); /*a caldo scelgo le coordinate random*/
        for (i = 0; i < N; i++)
            xx[i] *= 2;
        sprintf(file_name, "data_files/thermalisation/thermalisation_hot_N%d_prova.dat", N);
        fd = fopen(file_name, "w");
    }

    for (i = 0; i < N_TERM; i++)
    {
        fprintf(fd, "%d\t%f\n", i, action());
        acc_rate += sweep(3122000 + i);
    }

    fclose(fd);
    acc_rate /= (double)N_TERM;
    printf("L'accettanza media è: %f\n", acc_rate);

    return 0;
}
