
/*******************************************************************************
*
* File global.h
*
* Global parameters and arrays
*
* N numero di siti reticolari
* M massa dell'oscillatore per unità di passo reticolare
* W frequenza dell'oscillatore per unità di passo reticolare
* DELTA larghezza con cui vengono generati i numeri piatti nel Metropolis
* N_TERM numero di sweep da saltare per termalizzare la catena di Markov
* N_SWEEP numero di sweep totali della catena nel calcolo della funzione
*         di autocorrelazione
* T_M_MAX nei plot della funzione di autocorrelazione plottiamo fino a T_M_MAX
* T_P_MAX nei plot della funzione di autocorrelazione plottiamo la funzioni
*         per valori di t tra 1 e T_P_MAX
* D_BIN dimensione dei bin (circa 10 volte il tempo di autocorrelazione)
* N_BIN numero totale di bin
* ERR_REL_MAX errore relativo massimo accettato nel media in t dei valori di 
* energia ed elemento di matrice
*
*******************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

/*I valori attuali sono solo valori per testare il funzionamento dei programmi*/
#define N 64
#define M 1.0
#define W 1.0
#define DELTA 2.0
#define N_TERM 200
#define N_SWEEP 120000
#define T_M_MAX 15
#define T_P_MAX 10
#define D_BIN 100
#define N_BIN 10000
#define ERR_REL_MAX 0.1

#if defined MAIN_PROGRAM
  #define EXTERN
#else
  #define EXTERN extern
#endif

EXTERN double xx[N];

#undef EXTERN

#endif

