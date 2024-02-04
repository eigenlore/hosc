
/*******************************************************************************
 *
 * File harmonic_oscillator.c
 *
 * Contains all the functions necessary to implement the simulation of
 * quantum harmonic oscillator on lattice.
 *
 * All the functions are listed in "harmonic_oscillator.h" in the directory include.
 *
 * Author: Lorenzo Tasca
 *
 *******************************************************************************/

#define HARMONIC_OSCILLATOR_C
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "random.h"

/*Calcola l'azione dell'oscillatore armonico sulla configurazione atturale del Metropolis*/
double action()
{
	int i;
	double S;

	S = 0;

	for (i = 0; i < N; i++)
		S += (xx[(i + 1) % N] - xx[i]) * (xx[(i + 1) % N] - xx[i]) + W * W * xx[i] * xx[i];
	S *= M / 2;

	return S;
}

/*Calcola la variazione dell'azione data una variazione deltax della k-esima entrata della configurazione*/
double delta_action(int k, double deltax)
{
	double deltaS;
	deltaS = deltax * M / 2 * (-2 * xx[(k + 1) % N] + (2 * W * W + 4) * xx[k] - 2 * xx[(k - 1 + N) % N] + 2 * deltax + W * W * deltax);
	return deltaS;
}

/*Effettua uno sweep dato il seed del generatore pseudorandom passato come input. Aggiorna tutte le coordinate
della configurazione secondo l'algoritmo del Metropolis. Inoltre calcola l'accettanza dello sweep e la restituisce
come output*/
double sweep(int seed)
{
	double r[2 * N], deltaS, deltax, accettanza;
	int i;

	deltax = 0;
	accettanza = 0;

	rlxd_init(1, seed);
	ranlxd(r, 2 * N);
	for (i = 0; i < N; i++)
	{
		deltax = 2 * DELTA * (r[i] - 0.5); /*generazione piatta della nuova coordinata proposta*/
		deltaS = delta_action(i, deltax);
		if (deltaS < 0 || r[N + i] <= exp(-deltaS)) /*se deltaS<0 accetto con probabilità 1, altrimenti con probabilità e^-DeltaS*/
		{
			xx[i] = xx[i] + deltax;
			accettanza += 1.0;
		}
	}
	accettanza /= (double)N;
	return accettanza;
}

/*Inizializza il vettore a caldo o freddo e lo termalizza effettuando N_TERM sweep*/
void termalisation()
{
	int i, heat = 0; /*0 a freddo, 1 a caldo*/

	if (heat == 0)
	{
		for (i = 0; i < N; i++)
			xx[i] = 0; /*a freddo*/
	}
	else
	{
		rlxd_init(1, 1592000);
		ranlxd(xx, N); /*a caldo*/
	}

	for (i = 0; i < N_TERM; i++)
		sweep(3122000 + i); /*sweep per termalizzare*/

	return;
}

/*Calcola il correlatore c^k(t) effettuando N_SWEEP. Restituisce una matrice c[t_fisico][t_markoviano]*/
double **eval_c(int seed)
{
	double **c, *p;
	int i, k, t;

	c = (double **)malloc(N * sizeof(double *));
	p = (double *)malloc(N_SWEEP * N * sizeof(double));
	for (t = 0; t < N; t++)
		c[t] = p + t * N_SWEEP;

	/*fine dichiarazioni*/

	termalisation();

	for (i = 0; i < N_SWEEP; i++)
	{
		if ((i % (N_SWEEP / 10)) == 0)
			printf("%d sweep effettuati su %d\n", i, N_SWEEP); /*stampa periodicamente a schermo il numero di sweep già effettuati per monitorare l'esecuzione del programma*/
		for (t = 0; t < N; t++)
		{
			c[t][i] = 0;
			for (k = 0; k < N; k++)
				c[t][i] += xx[k] * xx[(k + t) % N];
			c[t][i] /= (double)N;
		}
		sweep(seed + i);
	}

	return c;
}

/*Calcola la media nei bin di c*/
double **bin_c(double **c)
{
	int t, k, i;
	double **c_binned, *p;
	c_binned = (double **)malloc(N * sizeof(double *));
	p = (double *)malloc(N_BIN * N * sizeof(double));
	for (t = 0; t < N; t++)
		c_binned[t] = p + t * N_BIN;

	/*fine dichiarazioni*/

	for (t = 0; t < N; t++)
	{
		for (k = 0; k < N_BIN; k++)
		{
			c_binned[t][k] = 0;
			for (i = k * D_BIN; i < (k + 1) * D_BIN; i++)
				c_binned[t][k] += c[t][i];
			c_binned[t][k] /= (double)D_BIN;
		}
	}

	return c_binned;
}

/*Calcola la funzione di autocorrelazione Gamma[t_fisico][t_markoviano]. Riceve in input il correlatore
c^k(t) e un booleano isBinned per stabilire se i dati sono stati precedentemente binnati.*/
double **eval_gamma(double **c, int isBinned)
{
	int i, t_M, t, N_TOT;
	double **gamma, sum1, sum2, *p;
	if (isBinned == 0) /*dati non binnati*/
		N_TOT = N_SWEEP;
	else /*dati binnati*/
		N_TOT = N_BIN;
	gamma = (double **)malloc(T_P_MAX * sizeof(double *));
	p = (double *)malloc(T_M_MAX * T_P_MAX * sizeof(double));
	for (t = 0; t < T_P_MAX; t++)
		gamma[t] = p + t * T_M_MAX;
	sum1 = 0;
	sum2 = 0;

	/*fine dichiarazioni*/

	for (t = 0; t < T_P_MAX; t++)
	{
		for (t_M = 0; t_M < T_M_MAX; t_M++)
		{
			for (i = 0; i < N_TOT - T_M_MAX; i++)
			{
				sum1 += c[t][i] * c[t][i + t_M];
				sum2 += c[t][i];
			}
			sum1 /= (double)(N_TOT - T_M_MAX);
			sum2 /= (double)(N_TOT - T_M_MAX);
			gamma[t][t_M] = sum1 - sum2 * sum2;
			sum1 = 0;
			sum2 = 0;
		}
	}

	return gamma;
}

/*Calcola da zero il correlatore c binnato, ovvero mediato sui bin. Il risultato della chiamata di questa funzione è equivalente
a chiamare in successione le funzioni eval_c e bin_c, tuttavia qua stiamo evitando di allocare un vettore grosso D_BIN*N_BIN,
rendendo l'esecuzione più leggera sulla memoria.*/
double **eval_c_binned(int seed)
{
	int t, k, i, count;
	double **c_binned, *p, c;

	c_binned = (double **)malloc(N * sizeof(double *));
	p = (double *)malloc(N_BIN * N * sizeof(double));
	for (t = 0; t < N; t++)
		c_binned[t] = p + t * N_BIN;
	for (t = 0; t < N; t++)
		for (k = 0; k < N_BIN; k++)
			c_binned[t][k] = 0;
	count = 1;

	/*fine dichiarazioni*/

	for (k = 0; k < D_BIN * N_BIN; k++)
	{
		if ((k % (D_BIN * N_BIN / 100)) == 0)
		{
			printf("Completamento al %d%%.\r", count); /*stampa periodicamente a schermo la percentuale di sweep effettuati per monitorare l'esecuzione del programma*/
			fflush(stdout);
			count++;
		}
		sweep(seed + k);
		for (t = 0; t < N; t++)
		{
			c = 0;
			for (i = 0; i < N; i++)
				c += xx[i] * xx[(i + t) % N];
			c /= (double)N;
			c_binned[t][(k - k % D_BIN) / D_BIN] += c / (double)D_BIN; 
		}
	}
	printf("\n");
	return c_binned;
}

/*Calcola la media sui bin del correlatore binnato*/
double *eval_c_bar(double **c_binned)
{
	int t, k;
	double *c_bar;

	c_bar = (double *)malloc(N * sizeof(double));

	/*fine dichiarazioni*/

	for (t = 0; t < N; t++)
	{
		c_bar[t] = 0;
		for (k = 0; k < N_BIN; k++)
		{
			c_bar[t] += c_binned[t][k];
		}
		c_bar[t] /= N_BIN;
	}

	return c_bar;
}

/*Calcolo arcocoseno iperbolico*/
double acosh(double x)
{
	return log(x + sqrt(x * x - 1));
}

/*Funzione accessoria per il calcolo dell'energia*/
double f(double a, double b, double c)
{
	return acosh((a + b) / (2.0 * c));
}

/*Funzione accessoria per il calcolo dell'elemento di matrice*/
double g(double c, int t, double deltaE)
{
	return sqrt(c * exp(t * deltaE));
}

/*Calcolo dei cluster Jackknife del correlatore*/
double **eval_cluster_c(double *c_bar, double **c_binned)
{
	int t, k;
	double **cluster_c, *p;

	cluster_c = (double **)malloc(N * sizeof(double *));
	p = (double *)malloc(N_BIN * N * sizeof(double));
	for (t = 0; t < N; t++)
		cluster_c[t] = p + t * N_BIN;

	/*fine dichiarazioni*/

	/*Calcolo dei cluster delle variabili primarie*/
	for (t = 0; t < N; t++)
	{
		for (k = 0; k < N_BIN; k++)
			cluster_c[t][k] = c_bar[t] - (c_binned[t][k] - c_bar[t]) / ((double)(N_BIN - 1));
	}

	return cluster_c;
}

/*Calcolo del gap energetico a partire dal correlatore*/
double *eval_deltaE(double *c_bar)
{
	int t;
	double *deltaE;

	deltaE = (double *)malloc(N * sizeof(double));

	/*fine dichiarazioni*/

	for (t = 0; t < N; t++)
		deltaE[t] = f(c_bar[(t + 1) % N], c_bar[(t - 1 + N) % N], c_bar[t]);

	return deltaE;
}

/*Calcolo dei cluster Jackknife del gap energetico*/
double **eval_cluster_deltaE(double **cluster_c)
{
	int t, k;
	double **cluster_deltaE, *p;

	cluster_deltaE = (double **)malloc(N * sizeof(double *));
	p = (double *)malloc(N_BIN * N * sizeof(double));
	for (t = 0; t < N; t++)
		cluster_deltaE[t] = p + t * N_BIN;

	/*fine dichiarazioni*/

	for (t = 0; t < N; t++)
	{
		for (k = 0; k < N_BIN; k++)
		{
			cluster_deltaE[t][k] = f(cluster_c[(t + 1) % N][k], cluster_c[(t - 1 + N) % N][k], cluster_c[t][k]);
		}
	}

	return cluster_deltaE;
}

/*Calcolo dell'elemento di matrice a partire dal correlatore*/
double *eval_matrix_el(double *c_bar, double *deltaE)
{
	int t;
	double *matrix_el;

	matrix_el = (double *)malloc(N * sizeof(double));

	/*fine dichiarazioni*/

	for (t = 0; t < N; t++)
		matrix_el[t] = g(c_bar[t], t, deltaE[t]);

	return matrix_el;
}

/*Calcolo dei cluster Jackknife dell'elemento di matrice*/
double **eval_cluster_matrix_el(double **cluster_c, double **cluster_deltaE)
{
	int t, k;
	double **cluster_matrix_el, *p;

	cluster_matrix_el = (double **)malloc(N * sizeof(double *));
	p = (double *)malloc(N_BIN * N * sizeof(double));
	for (t = 0; t < N; t++)
		cluster_matrix_el[t] = p + t * N_BIN;

	/*fine dichiarazioni*/

	for (t = 0; t < N; t++)
	{
		for (k = 0; k < N_BIN; k++)
		{
			cluster_matrix_el[t][k] = g(cluster_c[t][k], t, cluster_deltaE[t][k]);
		}
	}

	return cluster_matrix_el;
}

/*Calcola la varianza Jackknife di una variabile aleatoria x dati i suoi cluster*/
double *eval_error_jackknife(double *x, double **cluster_x)
{
	int t, k;
	double *error_x;

	error_x = (double *)malloc(N * sizeof(double));

	/*fine dichiarazioni*/

	for (t = 0; t < N; t++)
	{
		error_x[t] = 0;
		for (k = 0; k < N_BIN; k++)
			error_x[t] += (cluster_x[t][k] - x[t]) * (cluster_x[t][k] - x[t]);
		error_x[t] *= (double)(N_BIN - 1);
		error_x[t] /= (double)N_BIN;
		error_x[t] = sqrt(error_x[t]);
	}

	return error_x;
}

/*Effettua la media su t della variabile x(t) (x può essere deltaE o matrix_el), e propaga l'errore utilizzando i cluster Jackknife di x.
La media viene effettuata per i valori di x che hanno un errore relativo minore di ERR_REL_MAX. Restituisce un vettore la cui prima
entrata è il valore medio e la seconda il suo errore*/
double *eval_mean_x(double *x, double *error_x, double **cluster_x)
{
	int t, k;
	double *mean_x, *cluster_mean_x, err_rel, sum_weights;

	sum_weights = 0;
	mean_x = (double *)malloc(2 * sizeof(double)); /*mean_x[0] è il valore medio, mean_x[1] il suo errore*/
	cluster_mean_x = (double *)malloc(N_BIN * sizeof(double));
	mean_x[0] = 0;
	err_rel = 0;
	for (k = 0; k < N_BIN; k++)
		cluster_mean_x[k] = 0;

	/*fine dichiarazioni*/

	for (t = 1; err_rel < ERR_REL_MAX; t++)
	{
		sum_weights += 1 / (error_x[t] * error_x[t]);
		mean_x[0] += x[t] / (error_x[t] * error_x[t]); /*media pesata*/
		for (k = 0; k < N_BIN; k++)
			cluster_mean_x[k] += cluster_x[t][k] / (error_x[t] * error_x[t]);
		err_rel = error_x[t + 1] / x[t + 1];
	}

	mean_x[0] /= sum_weights;
	for (k = 0; k < N_BIN; k++)
		cluster_mean_x[k] /= sum_weights;

	mean_x[1] = 0;
	for (k = 0; k < N_BIN; k++)
		mean_x[1] += (cluster_mean_x[k] - mean_x[0]) * (cluster_mean_x[k] - mean_x[0]); /*varianza Jaccknife*/
	mean_x[1] *= (double)(N_BIN - 1);
	mean_x[1] /= (double)N_BIN;
	mean_x[1] = sqrt(mean_x[1]);

	return mean_x;
}
