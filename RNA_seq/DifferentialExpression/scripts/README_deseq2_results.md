# Come interpretare i risultati di DESeq2

Questa nota descrive i concetti dietro l'analisi di espressione differenziale e la tabella dei risultati. Non è una spiegazione riga per riga dello script.

## La domanda dell'analisi

Per ogni gene, DESeq2 stima se l'espressione differisce tra le condizioni sperimentali. Nell'analisi corrente il confronto è:

```text
tumore vs normale
```

La tabella principale è `de/results/deseq2_all_genes.csv`: una riga per ogni gene che ha superato il filtro di espressione bassa.

## Colonne della tabella

| Colonna | Significato pratico |
|---|---|
| `gene_id` | Identificativo Ensembl stabile del gene. |
| `gene_symbol` | Nome breve del gene, utile per la lettura e la ricerca bibliografica. |
| `baseMean` | Media dell'espressione **normalizzata** del gene su tutti i campioni. Indica il livello generale di espressione, non la differenza tra le condizioni. |
| `log2FoldChange` | Differenza stimata tra le condizioni, su scala log2. Nel confronto tumore vs normale: positivo = maggiore nel tumore; negativo = maggiore nel normale. `+1` significa circa 2 volte maggiore, `-1` circa 2 volte minore. |
| `lfcSE` | Errore standard del `log2FoldChange`: misura l'incertezza della differenza stimata. Più è piccolo, più la stima è precisa. |
| `stat` | Statistica del test usata per calcolare il p-value. Non è normalmente la prima colonna da interpretare biologicamente. |
| `pvalue` | Evidenza statistica non corretta per un singolo gene. |
| `padj` | p-value corretto per i molti geni testati (FDR/BH). È il criterio statistico principale; una soglia frequente è `padj < 0.05`. |

## `stat`, `pvalue` e `padj`: dal fold change alla significatività

Per ogni gene, il test parte dall'ipotesi nulla:

```text
non esiste una differenza media di espressione tra le condizioni
log2FoldChange = 0
```

### `stat`

La statistica del test mette in rapporto la dimensione dell'effetto con la sua incertezza. In modo approssimato:

```text
stat = log2FoldChange / lfcSE
```

Un valore vicino a zero indica che l'effetto è piccolo rispetto all'incertezza. Un valore molto positivo è evidenza di maggiore espressione nella prima condizione del contrasto; un valore molto negativo indica maggiore espressione nella seconda. In genere non si applica una soglia direttamente a `stat`: `pvalue` e `padj` sono più utili per la decisione statistica.

### `pvalue`

Il `pvalue` traduce `stat` nella domanda: **se il gene non avesse davvero alcuna differenza tra condizioni, quanto sarebbe raro osservare un risultato almeno così estremo?** Un valore basso è evidenza contro l'ipotesi di nessuna differenza. Non significa che sia la probabilità che l'ipotesi nulla sia vera.

Non è corretto selezionare direttamente `pvalue < 0.05` quando si testano migliaia di geni. Per un gene senza vera differenza, una soglia di 0.05 lascia circa il 5% di probabilità di un falso positivo dovuto al caso. Se tutti i 16.918 geni testati fossero invariati, ci si aspetterebbero in media:

```text
16.918 × 0.05 ≈ 846 falsi positivi
```

Il numero effettivo non deve essere esattamente 846 in una singola analisi, ma il problema dei molti test rimane.

### `padj`

`padj` è il `pvalue` corretto con il metodo Benjamini–Hochberg (BH), che controlla il False Discovery Rate (FDR). Usando `padj < 0.05`, tra i geni dichiarati significativi ci si aspetta che circa il 5% sia costituito da falsi positivi in media. Non identifica quali geni siano falsi positivi e non garantisce esattamente il 5% in ogni analisi.

La sequenza concettuale è:

```text
log2FoldChange + lfcSE → stat → pvalue → padj
```

## Perché `lfcSE` basso e `padj` alto possono coesistere

Sì, è possibile e ha un significato importante. Un `lfcSE` basso dice che il fold change è stimato con precisione; non dice che il fold change è diverso da zero o biologicamente grande.

```text
log2FoldChange = 0.03
lfcSE = 0.02
```

In questo esempio la stima è precisa, ma è anche quasi zero. La statistica è circa `0.03 / 0.02 = 1.5`, quindi il `pvalue` e probabilmente il `padj` saranno alti. L'interpretazione è: “la differenza è stata stimata bene e appare molto piccola”.

Un valore basso di `lfcSE` contribuisce a ottenere un `padj` basso soltanto se anche il valore assoluto di `log2FoldChange` è sufficientemente lontano da zero. Inoltre, `padj` dipende dalla correzione effettuata sull'insieme dei geni testati, non solo dal singolo gene.

## `baseMean`: perché serve

DESeq2 normalizza i conteggi per rendere confrontabili campioni con profondità di sequenziamento diverse. Il `baseMean` è la media di questi conteggi normalizzati.

```text
baseMean = media dei conteggi normalizzati del gene nei campioni
```

Un `baseMean` alto **non rende automaticamente corretto** un fold change. Indica però che il gene è stato osservato molte volte e quindi il rumore dovuto a pochi read casuali tende a pesare meno.

Per esempio, il passaggio da 1 a 2 read è un raddoppio (`log2FC = 1`), ma un singolo read in più o in meno può cambiare radicalmente il risultato. Il passaggio da 1.000 a 2.000 read è sempre un raddoppio, ma piccole fluttuazioni casuali incidono molto meno sul rapporto.

Anche un gene molto espresso può essere incerto se i campioni biologici mostrano differenze incoerenti. Per valutare la precisione del fold change bisogna quindi guardare soprattutto `lfcSE` e `padj`, non soltanto `baseMean`.

## `lfcSE`: che cosa considera

`lfcSE` è l'incertezza del `log2FoldChange` stimato dal modello. Diventa in generale più piccola quando la differenza tra condizioni è coerente tra i campioni e i dati contengono sufficiente informazione.

Tiene conto di:

- numero di campioni e replicati biologici;
- variabilità dei campioni dopo aver considerato il disegno sperimentale;
- livello di conteggio del gene, soprattutto per i geni poco espressi;
- dispersione dei conteggi RNA-seq, cioè la variabilità naturale oltre il semplice rumore di campionamento;
- covariate incluse nel disegno, come il paziente nell'analisi appaiata.

Un modello mentale utile è un intervallo approssimato:

```text
log2FoldChange ± 2 × lfcSE
```

Ad esempio, `log2FC = 2` e `lfcSE = 0.2` indicano una stima chiaramente positiva e precisa. Con `log2FC = 2` e `lfcSE = 1.5`, invece, la direzione e l'ampiezza dell'effetto sono molto più incerte.

Non esiste una soglia universale valida per `lfcSE`. Usalo per confrontare la precisione di effetti simili e per individuare fold change da interpretare con cautela. Per selezionare risultati statisticamente supportati, usa `padj`.

## Disegno appaiato: paziente + condizione

Il disegno dell'analisi è:

```text
~ patient + condition
```

In modo concettuale, per ciascun gene il modello separa:

```text
espressione osservata = effetto del paziente + effetto della condizione + variabilità residua
```

Questo significa che il confronto tumore–normale considera già le differenze di base tra pazienti. Non è un semplice confronto tra la media di tutti i tumori e la media di tutti i normali: confronta l'effetto della condizione all'interno di ciascun paziente e combina le differenze appaiate.

Se un paziente ha naturalmente espressione alta e un altro bassa, questa differenza può essere attribuita al termine `patient`, invece di essere confusa con l'effetto tumore/normalità. Di conseguenza, `lfcSE` rappresenta la precisione dell'effetto della condizione **dopo** aver considerato il paziente.

Questo non corregge automaticamente fattori che non sono nel disegno, come un batch tecnico confuso con la condizione, differenze sistematiche di qualità dell'RNA, contaminazione o etichette errate.

## Esperimenti con tre o più condizioni

Un `log2FoldChange` ha sempre bisogno di un confronto preciso. Con tre condizioni A, B e C non esiste un unico fold change universale. Occorrono contrasti distinti, per esempio:

```text
B vs A
C vs A
C vs B
```

Per ogni contrasto, il segno dipende dall'ordine: un valore positivo indica espressione maggiore nella prima condizione nominata rispetto alla seconda. Un test globale può dire che un gene differisce in almeno una condizione, ma non indica quale coppia di condizioni è diversa né la direzione; per questo servono i contrasti a coppie.

## Lettura pratica di un gene

Un risultato solido tende ad avere:

```text
padj basso       → evidenza statistica dopo correzione multipla
|log2FC| grande  → cambiamento potenzialmente rilevante biologicamente
lfcSE basso      → stima del cambiamento precisa
```

Questi elementi vanno interpretati insieme ai grafici PCA/heatmap, al controllo qualità e alla conoscenza biologica. La significatività statistica non dimostra da sola un meccanismo biologico o causalità.
