#ifndef HARMONIC_OSCILLATOR_H
#define HARMONIC_OSCILLATOR_H
#ifndef HARMONIC_OSCILLATOR_C

/*Calcola l'azione dell'oscillatore armonico sulla configurazione atturale del Metropolis*/
extern double action();

/*Calcola la variazione dell'azione data una variazione deltax della k-esima entrata della configurazione*/
extern double delta_action(int, double);

/*Effettua uno sweep dato il seed del generatore pseudorandom passato come input. Aggiorna tutte le coordinate
della configurazione secondo l'algoritmo del Metropolis. Inoltre calcola l'accettanza dello sweep e la restituisce
come output*/
extern double sweep(int);

/*Inizializza il vettore a caldo o freddo e lo termalizza effettuando N_TERM sweep*/
extern void termalisation();

/*Calcola il correlatore c^k(t) effettuando N_SWEEP. Restituisce una matrice c[t_fisico][t_markoviano]*/
extern double **eval_c(int);

/*Calcola la media nei bin di c*/
extern double **bin_c(double **);

/*Calcola la funzione di autocorrelazione Gamma[t_fisico][t_markoviano]. Riceve in input il correlatore
c^k(t) e un booleano isBinned per stabilire se i dati sono stati precedentemente binnati.*/
extern double **eval_gamma(double **, int);

/*Calcola da zero il correlatore c binnato, ovvero mediato sui bin. Il risultato della chiamata di questa funzione è equivalente
a chiamare in successione le funzioni eval_c e bin_c, tuttavia qua stiamo evitando di allocare un vettore grosso D_BIN*N_BIN,
rendendo l'esecuzione più leggera sulla memoria.*/
extern double **eval_c_binned(int);

/*Calcola la media sui bin del correlatore binnato*/
extern double *eval_c_bar(double **);

/*Calcolo arcocoseno iperbolico*/
extern double acosh(double);

/*Funzione accessoria per il calcolo dell'energia*/
extern double f(double, double, double);

/*Funzione accessoria per il calcolo dell'elemento di matrice*/
extern double g(double, int, double);

/*Calcolo dei cluster Jackknife del correlatore*/
extern double **eval_cluster_c(double *, double **);

/*Calcolo del gap energetico a partire dal correlatore*/
extern double *eval_deltaE(double *);

/*Calcolo dei cluster Jackknife del gap energetico*/
extern double **eval_cluster_deltaE(double **);

/*Calcolo dell'elemento di matrice a partire dal correlatore*/
extern double *eval_matrix_el(double *, double *);

/*Calcolo dei cluster Jackknife dell'elemento di matrice*/
extern double **eval_cluster_matrix_el(double **, double **);

/*Calcola la varianza Jackknife di una variabile aleatoria x dati i suoi cluster*/
extern double *eval_error_jackknife(double *, double **);

/*Effettua la media su t della variabile x(t) (x può essere deltaE o matrix_el), e propaga l'errore utilizzando i cluster Jackknife di x.
La media viene effettuata per i valori di x che hanno un errore relativo minore di ERR_REL_MAX. Restituisce un vettore la cui prima
entrata è il valore medio e la seconda il suo errore*/
extern double *eval_mean_x(double *, double *, double **);

#endif
#endif
