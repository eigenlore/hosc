
Spiegazione del contenuto dei file dati

=============================================================================

Directory thermalisation:

 - thermalisation_hot_N*.dat contiene l'azione stampata a ogni sweep prima 
   della termalizzazione, partendo da caldo con N=*.
   File generato da term.c

 - thermalisation_cold_N*.dat contiene l'azione stampata a ogni sweep prima   
   della termalizzazione, partendo da freddo con N=*
   File generato da term.c

 - accettanza_N*.dat contiene l'accettanza media con N=*, al variare di Delta
   Per N!=64 è riportata solo quella per Delta=1. L'accettanza è calcolata
   partendo da freddo, ma quella a caldo è sostanzialmente la stessa. 

=============================================================================

Directory gamma:

 - gamma_t*_N?.dat contiene la funzione Gamma al variare di t_M valutata a t=*, 
   con N=?. gamma_binned_t*_N?.dat analogo ma dopo il binnaggio.
   File generati da decorrelation_time.c  
   Di seguito sono riportati i valori usati per le simulazioni:

	- N=32		 D_BIN=20	N_SWEEP=1200000		T_M_MAX=15
	- N=64		 D_BIN=40	N_SWEEP=1200000		T_M_MAX=15
	- N=128		 D_BIN=90	N_SWEEP=1350000		T_M_MAX=25
	- N=256		 D_BIN=250	N_SWEEP=1250000		T_M_MAX=35
	- N=512		 D_BIN=800	N_SWEEP=800000		T_M_MAX=50	
	- N=1024	 D_BIN=2500	N_SWEEP=250000		T_M_MAX=100

 - tau_correl_N*.dat contiene il tempo di decorrelazione di Gamma(t,t_M), al
   variare di t, per N=*. 

=============================================================================
   
Directory energy:

 - I parametri usati nelle simulazioni sono 

	N	SEED		BIN		DURATA SIMULAZIONE	
	==========================================================
	32	4122023         5 milioni	24 min
	64	15092000        2,5 milioni	32 min
	128	13121989        1 milione	53 min
	256	13062017	500 mila	215 min
	512	6122023		50 mila		200 min
	1024	20589712	10 mila		418 min

 - deltaE_N*.dat contiene l'energia per unità di passo reticolare, al variare
   di t, calcolata per N=*. 

 - matrix_el_N*.dat è analogo ma con l'elemento di matrice.

 - mean_deltaE_N*.dat contiene l'energia per unità di passo reticolare media, 
   con il relativo errore.

 - mean_matrix_el_N*.dat analogo ma con l'elemento di matrice.
