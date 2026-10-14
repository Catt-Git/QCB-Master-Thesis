# Interpretazione delle dimensioni latenti DRVI — `drvi_tum` e `drvi_epi_cnv`

Risponde alla domanda: **come interpreto le varie DR ottenute, cosa contengono?**, sui due
run principali della fase 05 (`05_drvi_tumoral_epi`):

- **`drvi_tum`** → run `drvi_tum_64_nomt` (64 dimensioni latenti, pannello HVG senza geni
  mitocondriali), sul sottoinsieme maligno (`CELL_SET=tum`, 19 coorti, **42 096 cellule**,
  35 cluster leiden `optscib_tum_leiden`).
- **`drvi_epi_cnv`** → run `drvi_epicnv_64_nomt`, sul controllo epiteliale completo dopo la
  chiamata CNV (`CELL_SET=epi`, 25 coorti, **60 233 cellule**, 44 cluster leiden), che include
  sia le cellule maligne sia quelle non maligne re-annotate con CellTypist.

Questo documento riporta **solo le DR "sicure"**: quelle in cui il programma genico è
confermato indipendentemente da due rotte di evidenza diverse (o, per l'identità di
cluster/tipo cellulare, da un match SMI pulito ≥ 0.7), **senza** segnali di confondimento
tecnico. Tutto il resto (assi confondati, singole rotte, target composti) è elencato a parte
come "cosa è stato scartato", perché il progetto stesso tratta quei casi come *lead*, non come
risultati.

---

## 1. Come è stata giudicata la "sicurezza" di una DR

Ogni dimensione DRVI produce due direzioni indipendenti (`DR n+`, `DR n-`, il lato positivo e
quello negativo del fattore), quindi ogni run di 64 dimensioni genera 128 righe da valutare.
Per ciascuna riga, la pipeline (`05_6_cell_first` → `05_7_factor_first` → `05_8_convergence`)
applica due rotte che possono sbagliare in modi diversi e quindi si controllano a vicenda:

| Rotta | Cosa chiede | Soglia |
|---|---|---|
| **A — cell-first** | le cellule con punteggio alto su questa DR hanno anche un punteggio alto sulla firma genica X? (correlazione di Spearman, per-coorte) | \|ρ\| ≥ 0.20 |
| **B — factor-first** | i geni con carico più alto su questa DR sono arricchiti per la firma X? (over-representation analysis sui geni top del decoder) | FDR globale < 0.05 |

Una riga è **`convergent`** solo se A e B superano la soglia **sulla stessa famiglia di
firme** (es. entrambe "stemness", non necessariamente la stessa lista). Questo è il criterio
usato in tutta la fase per chiamare una DR uno stato cellulare vero: né A né B da soli sono
sufficienti (A può convergere per correlazione tra punteggi molto simili tra loro, B per
sovrapposizione di geni senza controparte cellulare).

Oltre a questo, ogni riga convergente viene controllata contro due **confondenti tecnici**:
correlazione con la profondità di sequenziamento (`rho_n_genes_by_counts`) e con il ciclo
cellulare (`rho_S_score` / `rho_G2M_score`), entrambe con soglia 0.30. **Una DR convergente ma
con un flag non è considerata un risultato in questo report**, è un *lead* — coerente con come
la fase stessa la tratta (`05_8_convergence/README.md`).

Per l'attribuzione di **cluster / tipo cellulare** si usa una terza fonte, indipendente dalle
prime due: la **SMI** (*Specific Matching Information*, calcolata con `drvi.utils.metrics` in
`05_3_drvi_run`), che misura quanto una direzione di una DR corrisponde in modo pulito a UNA
sola categoria discreta (un cluster leiden o un tipo cellulare). Nella fase si usa soglia
**0.7** per i cluster leiden (partizione delle cellule del run stesso → match quasi 1:1 attesi)
e soglia **0.4** per le etichette "in prestito" (CellTypist, stato CNV). Solo le DR sopra
soglia sono riportate qui come identità di cluster/tipo cellulare.

**Una DR può comparire in entrambe le categorie** (programma genico confermato E identità di
cluster), in una sola, o in nessuna delle due (la maggioranza: su 128 righe per run, solo
39 per `drvi_tum` e 29 per `drvi_epi_cnv` superano il criterio di convergenza pulita).

---

## 2. Cosa NON è incluso (e perché)

- **Gli assi "target" delle collezioni SCIE ed EMT sono tutti confondati o assenti.**
  SCIE cerca un asse che sia *sia* stem-like *sia* immune-evasivo sullo stesso verso; su
  entrambi i run gli unici assi che soddisfano il criterio portano il flag `depth,cell_cycle`
  (`DR 60+`, `DR 61+`, `DR 59+`, `DR 16-` su `drvi_tum_64_nomt`; `DR 18-` su
  `drvi_epicnv_64_nomt`) — quindi **nessuno stato "stem + immune-evasivo" pulito emerge da
  questi run**, e non è incluso sotto. EMT non scrive nessun asse target in nessuno dei due run:
  nessuna dimensione soddisfa *insieme* i due criteri (epiteliale-alto e mesenchimale-alto sullo
  stesso lato) — la lettura corretta, spiegata nel README della fase, è che lo stato ibrido è
  un'**intersezione di due assi**, non una singola direzione dello spazio latente.
- **Le dimensioni del ciclo cellulare sono le correlazioni più forti del dataset, ma sono
  escluse per costruzione.** Su `drvi_tum_64_nomt`, `DR 3+` = `MP1_CELL_CYCLE_G2_M` (ρ = 0.845,
  la riga più forte dell'intero run) e `DR 5+` = `MP2_CELL_CYCLE_G1_S` (ρ = 0.775) sono
  entrambe `convergent` ma portano il flag `cell_cycle`; lo stesso vale su
  `drvi_epicnv_64_nomt` (`DR 10-` ρ = 0.772, `DR 7+` ρ = 0.733). Non sono "sbagliate": sono la
  domanda che la fase **non ancora scritta** `05_9_cycle_confound` esiste per rispondere
  (quanto di "staminalità"/proliferazione è ciclo cellulare puro).
- **Le firme EMT di Gavish (`MP12`/`MP13`/`MP14`) sono sistematicamente flaggate**
  `mesenchymal_may_be_ambient_or_doublet` — un flag specifico della collezione, non
  `depth`/`cell_cycle`. Sono lead per un possibile artefatto di ambient RNA/doublet coi
  fibroblasti, non stati mesenchimali confermati.
- Le due tabelle `target_axes` (dove esistono) e le righe convergenti-ma-flaggate sono comunque
  su disco (`tables/<collezione>/<run_id>/`) per chi vuole vederle: qui sono omesse perché la
  richiesta era "solo quelle sicure".

---

## 3. `drvi_tum` — run `drvi_tum_64_nomt` (42 096 cellule maligne, 35 cluster)

`cell_type` è costante (`malignant`) su questo oggetto, quindi — come richiesto — l'identità di
cluster usa **`optscib_tum_leiden`** al posto del tipo cellulare.

### 3.1 DR con programma genico confermato (Route A + Route B convergenti, nessun flag)

39 direzioni su 128 (30%). Raggruppate per tema biologico; "collezione" indica quale/i
collezione/i ha confermato la convergenza indipendentemente (SCIE = stemness/immunogenicità,
EMT = ibrido epiteliale-mesenchimale, Gavish-TNBC = metaprogrammi pan-cancro).

**Staminalità mammaria** (`LIM_STEM` = MSigDB *LIM_MAMMARY_STEM_CELL_UP*)

| DR | ρ (Route A) | FDR (Route B) | Collezione |
|---|---|---|---|
| DR 1+ | 0.348 | 1.3e-21 | SCIE |
| DR 19- | 0.317 | 1.5e-05 | SCIE |
| DR 21+ | 0.225 | 8.1e-03 | SCIE |
| DR 26- | 0.275 | 3.4e-14 | SCIE |
| DR 27- | 0.208 | 2.9e-04 | SCIE |
| DR 39- | 0.245 | 2.9e-04 | SCIE |
| DR 49+ | 0.221 | 3.9e-03 | SCIE |

Sette direzioni distinte portano la stessa firma di staminalità — non un'unica DR "stemness",
ma un piccolo fascio di assi che catturano la stessa biologia da lati diversi dello spazio
latente (nessuna condivide una firma B diversa dalla A, cioè `same_signature = True` per tutte).

**Interferone / immunogenicità** (`HALLMARK_IFNA`, `HALLMARK_IFNG`, `ISDS`,
`IMMUNOGENIC_CONSENSUS`, e i metaprogrammi `MP17`/`MP18_INTERFERON_MHC_II`)

| DR | Firma A | ρ | Firma B | FDR | Collezione |
|---|---|---|---|---|---|
| DR 18- | IMMUNOGENIC_CONSENSUS | 0.729 | HALLMARK_IFNG | 1.6e-14 | **SCIE + Gavish** (MP17, ρ 0.579, FDR 3.2e-25) |
| DR 23- | HALLMARK_IFNA | 0.472 | HALLMARK_IFNG | 1.3e-40 | SCIE |
| DR 30+ | IMMUNOGENIC_CONSENSUS | 0.271 | HALLMARK_IFNG | 5.3e-04 | SCIE |
| DR 31- | IMMUNOGENIC_CONSENSUS | 0.523 | HALLMARK_IFNG | 8.4e-10 | SCIE |
| DR 32+ | HALLMARK_IFNA | 0.287 | HALLMARK_IFNA | 3.2e-06 | SCIE |
| DR 34- | ISDS | 0.561 | HALLMARK_IFNA | 1.8e-36 | SCIE |
| DR 35- | HALLMARK_IFNA | 0.419 | HALLMARK_IFNG | 2.1e-34 | SCIE |
| DR 61- | IMMUNOGENIC_CONSENSUS | 0.369 | HALLMARK_IFNA | 1.4e-02 | **SCIE + Gavish** (MP18, ρ 0.407, FDR 1.4e-06) |
| DR 2- | MP17_INTERFERON_MHC_II_1 | 0.227 | MP18_INTERFERON_MHC_II_2 | 2.2e-04 | Gavish |
| DR 14- | MP18_INTERFERON_MHC_II_2 | 0.207 | MP18_INTERFERON_MHC_II_2 | 8.5e-06 | Gavish |
| DR 23- | MP17_INTERFERON_MHC_II_1 | 0.415 | MP17_INTERFERON_MHC_II_1 | 2.2e-28 | Gavish (stessa DR di sopra, doppia conferma) |
| DR 34- | MP17_INTERFERON_MHC_II_1 | 0.447 | MP17_INTERFERON_MHC_II_1 | 9.8e-24 | Gavish (stessa DR di sopra, doppia conferma) |
| DR 35- | MP17_INTERFERON_MHC_II_1 | 0.359 | MP17_INTERFERON_MHC_II_1 | 2.7e-22 | Gavish (stessa DR di sopra, doppia conferma) |
| DR 64- | MP17/MP18_INTERFERON_MHC_II | 0.236 | MP18_INTERFERON_MHC_II_2 | 2.2e-04 | Gavish (dimensione "vanished") |

`DR 18-`, `DR 23-`, `DR 34-`, `DR 35-` e `DR 61-` sono confermate **indipendentemente da due
collezioni curate in modo completamente diverso** (le firme IFN/immunogenicità del laboratorio
vs. i metaprogrammi pan-cancro di Gavish): è il gruppo di DR più solido del run. Il segno
negativo domina: il polo "immunogenico/IFN-alto" di queste dimensioni è quasi sempre il lato
`-`.

**Stress / ipossia / metabolismo**

| DR | Firma (Gavish) | ρ | FDR |
|---|---|---|---|
| DR 8+ | MP5_STRESS | 0.811 | 1.1e-38 |
| DR 44- | MP5_STRESS | 0.365 | 6.4e-17 |
| DR 53- | MP5_STRESS | 0.297 | 1.1e-15 |
| DR 13- | MP6_HYPOXIA | 0.575 | 6.1e-31 |
| DR 7- | MP7_STRESS_IN_VITRO | 0.268 | 1.5e-08 |
| DR 50+ | MP7_STRESS_IN_VITRO | 0.471 | 2.0e-16 |
| DR 4+ | MP21_RESPIRATION | 0.469 | 2.0e-03 |
| DR 29+ | MP21_RESPIRATION | 0.273 | 3.6e-02 |
| DR 9+ | MP8_PROTEASOMAL_DEGRADATION | 0.447 | 3.8e-05 |
| DR 43- | MP39_METAL_RESPONSE | 0.416 | 9.8e-06 |

`DR 8+` è la seconda correlazione più forte del run (0.811) dopo le dimensioni del ciclo
cellulare escluse sopra. `MP7_STRESS_IN_VITRO` è il metaprogramma che la fase tiene apposta per
distinguere uno stress biologico da un artefatto di dissociazione (vedi README fase 05): qui
compare come firma **secondaria** su `DR 7-`/`DR 50+`, non come firma dominante di nessuna DR —
lo stress catturato da `MP5`/`MP6` non sembra un semplice artefatto di dissociazione.

**Secrezione** (`MP22_SECRETED_1`, `MP23_SECRETED_2`)

| DR | Firma | ρ | FDR |
|---|---|---|---|
| DR 21+ | MP22_SECRETED_1 | 0.461 | 4.5e-14 |
| DR 29- | MP22_SECRETED_1 | 0.383 | 1.6e-15 |
| DR 9- | MP23_SECRETED_2 | 0.356 | 5.4e-07 |
| DR 27+ | MP23_SECRETED_2 | 0.242 | 9.1e-05 |
| DR 30+ | MP23_SECRETED_2 | 0.484 | 4.8e-14 |

**Senescenza epiteliale** (`MP19_EPITHELIAL_SENESCENCE`)

| DR | ρ | FDR |
|---|---|---|
| DR 36- | 0.337 | 4.3e-16 |
| DR 40- | 0.312 | 1.6e-13 |
| DR 41- | 0.281 | 4.4e-22 |
| DR 48- | 0.279 | 5.0e-09 |
| DR 51- | 0.343 | 8.5e-15 |

Cinque direzioni distinte per lo stesso metaprogramma — come la staminalità, non un singolo
asse ma un fascio.

**Stato ibrido epiteliale-mesenchimale** (`EMT_A_HYBRID` + `EMT_C_HYBRID`, le due liste
indipendenti della collezione EMT)

| DR | ρ (A_HYBRID) | FDR (C_HYBRID) |
|---|---|---|
| DR 21+ | 0.224 | 4.4e-02 |
| DR 32- | 0.248 | 8.1e-03 |

Solo due direzioni, deboli (ρ ≈ 0.22–0.25, vicino alla soglia): coerente con quanto scritto nel
README della fase — la EMT parziale nel carcinoma è un continuum con una coda mesenchimale
sparsa, non tre cluster netti. `DR 21+` è la stessa dimensione che converge anche su
`MP22_SECRETED_1` (tabella secrezione sopra): un fenotipo secretorio e uno ibrido-EMT deboli
coesistono sullo stesso lato di questa DR.

### 3.2 DR = identità di un cluster (match SMI ≥ 0.7 su `optscib_tum_leiden`)

Queste 4 direzioni non hanno (o non hanno bisogno di) un programma genico da Route A/B: SMI
dice che la direzione corrisponde quasi **1:1** a un singolo cluster leiden, quasi sempre
piccolo. I geni top e il termine Enrichr più significativo vengono da
`factor_processes_drvi_tum_64_nomt_top100.csv`.

| DR | Cluster | SMI | n cellule nel cluster | Geni top | Lettura |
|---|---|---|---|---|---|
| DR 55- | 27 | 0.993 | 24 | `IGLV3-21, IGLV3-9, IGHV3-48, IGHV3-21, IGHV3-15, FCRL5, IGHG1` | Immunoglobuline — molto probabilmente RNA ambientale da plasmacellule/B, non un programma tumorale |
| DR 17- | 24 | 0.855 | 598 | `PRR27, CSN3, SELENBP1, ELF3, LRRC26, CSN1S1, CLDN3, ZG16B` | Caseine (`CSN3`, `CSN1S1`) — differenziazione secretoria/lattazionale |
| DR 38+ | 31 | 0.848 | 979 | `HBB, HBA2, HBA1, ALAS2, HBM, SNCA, TRIM58, HEMGN` | Emoglobine — contaminazione da globuli rossi (RNA ambientale), non biologia tumorale |
| DR 39+ | 26 | 0.723 | 15 | `CALML5` (+ *Epidermis Development*) | Differenziazione squamosa/epidermica, cluster minuscolo |

Nessuna di queste quattro coincide con le DR "a programma" della sezione 3.1: sono un tipo di
finding diverso — non stati tumorali ricorrenti ma **identità di cluster molto specifiche**,
due delle quali (`DR 55-`, `DR 38+`) hanno la firma di contaminazione da RNA ambientale
piuttosto che di uno stato biologico del tumore.

---

## 4. `drvi_epi_cnv` — run `drvi_epicnv_64_nomt` (60 233 cellule, 44 cluster)

Qui `cell_type` (etichetta post-CNV) **non** è costante — include sia `malignant` sia le
etichette CellTypist reali delle cellule non maligne — quindi, dove disponibile, l'attribuzione
usa direttamente il tipo cellulare; `optscib_tum_leiden` resta il riferimento per le DR senza
match di tipo cellulare.

### 4.1 DR con programma genico confermato (Route A + Route B convergenti, nessun flag)

29 direzioni su 128 (23%).

**Staminalità / progenitori**

| DR | Firma A | ρ | Firma B | FDR |
|---|---|---|---|---|
| DR 2+ | LIM_STEM | 0.290 | LIM_STEM | 1.1e-12 |
| DR 5- | LIM_STEM | 0.361 | LIM_STEM | 9.4e-59 |
| DR 23- | LIM_STEM | 0.355 | LIM_STEM | 6.7e-09 |
| DR 29- | LIM_STEM | 0.230 | LIM_STEM | 9.0e-06 |
| DR 6- | FMASC | 0.360 | FMASC | 2.8e-02 |
| DR 10- | FMASC | 0.291 | ESC_WONG | 3.0e-26 |

**Interferone / immunogenicità**

| DR | Firma A | ρ | Firma B | FDR | Nota |
|---|---|---|---|---|---|
| DR 18+ | HALLMARK_IFNG | 0.256 | HALLMARK_IFNG | 3.1e-02 | |
| DR 24- | HALLMARK_IFNA | 0.487 | HALLMARK_IFNG | 1.0e-38 | **conferma anche in Gavish**: MP17, ρ 0.436, FDR 1.7e-23 |
| DR 26- | IMMUNOGENIC_CONSENSUS | 0.699 | HALLMARK_IFNG | 3.2e-17 | **conferma anche in Gavish**: MP17, ρ 0.578, FDR 7.9e-27 |
| DR 34- | ISDS | 0.551 | HALLMARK_IFNA | 1.4e-34 | **conferma anche in Gavish**: MP17, ρ 0.461, FDR 1.6e-20 |
| DR 48- | HALLMARK_IFNG | 0.272 | HALLMARK_IFNG | 2.2e-05 | |
| DR 54+ | HALLMARK_IFNA | 0.385 | HALLMARK_IFNA | 6.2e-13 | vedi anche DR54+ = MP18/17 sotto |

`DR 24-`, `DR 26-` e `DR 34-` sono confermate da **due collezioni indipendenti** (come nel run
`tum`): sono le DR più solide di questo run per l'asse interferone/immunogenicità.

**Stress / ipossia / secrezione / senescenza** (Gavish)

| DR | Firma | ρ | FDR |
|---|---|---|---|
| DR 15- | MP5_STRESS | 0.784 | 1.3e-33 |
| DR 1+ | MP5_STRESS | 0.319 | 1.2e-11 |
| DR 46+ | MP5_STRESS | 0.363 | 3.2e-15 |
| DR 17- | MP6_HYPOXIA | 0.536 | 1.9e-26 |
| DR 40+ | MP7_STRESS_IN_VITRO | 0.298 | 4.1e-10 |
| DR 3- | MP23_SECRETED_2 | 0.305 | 8.6e-09 |
| DR 18+ | MP23_SECRETED_2 | 0.398 | 2.5e-06 |
| DR 21- | MP23_SECRETED_2 | 0.385 | 8.6e-09 |
| DR 42+ | MP23_SECRETED_2 | 0.251 | 6.6e-08 |
| DR 48- | MP23_SECRETED_2 | 0.393 | 1.2e-11 |
| DR 51+ | MP23_SECRETED_2 | 0.265 | 6.6e-08 |
| DR 35- | MP22_SECRETED_1 | 0.395 | 1.8e-18 |
| DR 52- | MP22_SECRETED_1 | 0.288 | 1.1e-13 |
| DR 9- | MP19_EPITHELIAL_SENESCENCE | 0.397 | 3.1e-22 |
| DR 38+ | MP19_EPITHELIAL_SENESCENCE | 0.302 | 7.6e-18 |
| DR 43- | MP19_EPITHELIAL_SENESCENCE | 0.236 | 1.1e-07 |
| DR 44+ | MP39_METAL_RESPONSE | 0.405 | 1.1e-04 |
| DR 54+ | MP18_INTERFERON_MHC_II_2 | 0.336 | — (vs MP17, ρ 1.6e-16) | |

**Stato "epiteliale puro"** (`EMT_B_EPITHELIAL`, il polo epiteliale della collezione EMT — non
uno stato ibrido)

| DR | ρ | FDR |
|---|---|---|
| DR 15- | 0.265 | 2.8e-02 |
| DR 18+ | 0.221 | 2.8e-02 |

Da leggere con attenzione: qui la convergenza è sul **polo epiteliale** della collezione EMT
(differenziazione epiteliale forte), non sul polo ibrido/mesenchimale — a differenza del run
`tum`, in questo run non emerge nessuno stato ibrido pulito. `DR 15-` e `DR 18+` sono anche le
DR dominanti per `MP5_STRESS` e `MP23_SECRETED_2` rispettivamente (tabelle sopra): la
differenziazione epiteliale coesiste con quei programmi sugli stessi assi.

### 4.2 DR = identità di tipo cellulare (match SMI ≥ 0.4 su `cell_type`)

A differenza di `drvi_tum`, qui il tipo cellulare non è costante e tre direzioni hanno un match
diretto, tutte con geni coerenti con l'etichetta:

| DR | Tipo cellulare | SMI | Geni top | Coerenza |
|---|---|---|---|---|
| DR 1+ | LummHR-major (luminale ormono-responsivo) | 0.625 | `TFF1, STC2, PI15, AGR3, PHGR1, ESR1, AGR2, AZGP1` | **Alta** — `ESR1`, `TFF1`, `AGR2` sono marcatori canonici luminal HR+ |
| DR 3- | Lumsec-major (luminale secretorio) | 0.576 | `SCGB3A1, SNORC, CSN1S1, AQP5, PIGR, TMEM213, FAM3D, KIT` | **Alta** — `SCGB3A1`, `KIT`, `PIGR` sono marcatori canonici luminal-secretorio |
| DR 5- | basale/mioepiteliale | 0.513 | `C2orf40, MYLK, MYH11, ACTA2, TAGLN, CXCL14, GRP, TPM2` | **Alta** — `MYH11`, `ACTA2`, `MYLK`, `TAGLN` sono geni contrattili mioepiteliali classici |

Queste tre DR sono anche in tabella 4.1 (`DR 3-` → `MP23_SECRETED_2`, `DR 5-` → `MP6_HYPOXIA`):
in questo run, quindi, alcune DR portano **sia** un'identità di lineage pulita **sia** un
programma biologico convergente — la lettura più forte possibile per una DR in questo report.

Un quarto match, **`DR 1+` → stato CNV**, è più debole e ambiguo: SMI 0.468 con
`non_malignant` e 0.419 con `malignant` (soglia 0.4, quindi sopra ma vicina, e sui *due* lati)
— `DR 1+` separa parzialmente maligno da non maligno ma non in modo pulito 1:1; da leggere come
un accenno, non come un'identità di stato CNV.

### 4.3 DR = identità di cluster leiden (match SMI ≥ 0.7 su `optscib_tum_leiden`)

| DR | Cluster | SMI | n cellule | Label dominante nel cluster | Geni top | Lettura |
|---|---|---|---|---|---|---|
| DR 23+ | 32 | 0.861 | 1 028 | malignant (93%) | `PRR27, CSN3` | Lattazione/secrezione |
| DR 16+ | 35 | 0.859 | 1 665 | malignant (61%) | `HBB, HBA2, HBA1, ALAS2, CA1, HBM, HBD, SLC4A1` | Emoglobine — contaminazione RBC |
| DR 5- | 6 | 0.808 | 1 525 | basal (94%) | `MYLK, MYH11, ACTA2, TAGLN, CXCL14` | Stesso cluster del match cell-type `DR 5-` sopra: coerenza interna |
| DR 8- | 41 | 0.795 | 2 724 | malignant (75%) | `IFI27` (isolato) | Interferone-stimolato, poco supportato (1 solo gene sopra soglia) |
| DR 39+ | 26 | 0.778 | 369 | basal (99%) | `UCP2, IGF1, TMSB4X, CORO1A, BASP1, AIF1, LCP1, SPI1` | Geni mieloidi/immunitari (`AIF1`, `SPI1`, `LCP1`) — probabile infiltrato immunitario nel cluster basale, non biologia epiteliale |
| DR 37- | 38 | 0.743 | 24 | malignant (100%) | `IGKV3-15, IGHV1-2, IGLV3-1, IGKC, IGHA1` | Immunoglobuline — ambient RNA da plasmacellule, cluster minuscolo |
| DR 9- | 31 | 0.736 | 26 | malignant (100%) | `CST6, SBSN, SPRR1B, KLK10, KLK7, SPRR1A, KLK5, KLK6` | Cheratinizzazione/differenziazione squamosa, cluster minuscolo |
| DR 41- | 26 | 0.716 | 369 | basal (99%) | `HIST3H2A, ENO3, FSCN1, CDH13, PCOLCE, CKS1B, CPE, PHGDH` | Misto, nessun tema dominante chiaro |

Pattern che si ripete dal run `tum`: i match SMI più puliti su cluster piccoli sono quasi
sempre firme di contaminazione (emoglobine, immunoglobuline) piuttosto che stati tumorali.

### 4.4 Un'identità aggiuntiva: contaminazione/doppietto immunitario

| DR | Match | SMI (soglia 0.4, etichetta in prestito) | Geni top | Lettura |
|---|---|---|---|---|
| DR 12- | CD8-Tem (CellTypist, etichetta pre-CNV) | 0.815 | `CD3E, NKG7, CD69, GNLY, CD3D, GZMK, GZMA, CD7` | Firma T CD8 completa e pulita (tutti marcatori canonici) — cellule T infiltranti catturate nell'oggetto epiteliale, non un programma epiteliale/tumorale |

---

## 5. Temi biologici che ricorrono in entrambi i run

**Avvertenza:** `drvi_tum` e `drvi_epi_cnv` sono due modelli DRVI addestrati indipendentemente
(celle diverse, in parte HVG diversi), quindi **lo stesso numero di DR nei due run non è la
stessa dimensione biologica** — la corrispondenza va letta a livello di *programma*, non di
indice. Detto questo, quattro temi emergono come DR sicure, indipendenti, in **entrambi** i
compartimenti:

1. **Interferone/MHC-II (`HALLMARK_IFNA`/`IFNG`, `IMMUNOGENIC_CONSENSUS`, `MP17`/`MP18`)** — il
   tema più robusto in assoluto: 5 DR in `tum` e 3 in `epi_cnv`, quasi tutte confermate da due
   collezioni curate indipendentemente.
2. **Staminalità mammaria (`LIM_STEM`, `FMASC`)** — 7 DR in `tum`, 6 in `epi_cnv`.
3. **Stress/ipossia (`MP5`, `MP6`, `MP7`)** e **secrezione (`MP22`, `MP23`)** — presenti con
   forza comparabile in entrambi.
4. **Senescenza epiteliale (`MP19`)** — 5 DR in `tum`, 3 in `epi_cnv`.

Due differenze notabili:

- Solo `drvi_tum` produce DR convergenti sul polo **ibrido** EMT (`EMT_A/C_HYBRID`, seppur
  deboli); `drvi_epi_cnv` converge invece sul polo **epiteliale puro** (`EMT_B_EPITHELIAL`).
  Non contraddittorio — sono letture compatibili con un continuum E→M che il compartimento
  `tum` (solo cellule aneuploidi) intercetta più vicino al lato mesenchimale.
- Solo `drvi_epi_cnv` ha DR con un'attribuzione di **lineage pulita** (`LummHR-major`,
  `Lumsec-major`, basale): atteso, perché solo qui `cell_type` non è costante.

---

## 6. Limiti che si applicano a tutto quanto sopra

Ereditati dai CAVEAT scritti in ogni tabella della fase (`utils/signature_common.py`):

- Le cellule "maligne" sono definite da una **soglia** sul punteggio CNV per coorte
  (`05_1_infercnv`), non da un clone validato — niente normale abbinato, niente WGS bulk.
  19 delle 29 coorti epiteliali raggiungono la soglia dei 200 cellule richiesta.
- Le etichette CellTypist "pre-CNV" (`cell_type_01_4`) vengono da un atlante di mammella
  **normale**, senza classe maligna: sono landmark, non identità — per questo, dove possibile,
  questo report nomina le DR con i loro geni (Route A/B) o con `optscib_tum_leiden`/`cell_type`
  post-CNV, non con quell'etichetta.
- Nessun controllo di embedding indipendente è stato ancora eseguito su questa fase (l'analogo
  di `04_9_embedding_control` non esiste ancora per la 05): la stabilità di queste DR rispetto a
  un seed diverso o a un'altra architettura non è stata testata.
- I match SMI sono descrittivi (numero di cellule, non un test statistico con p-value): un
  match su un cluster da 15–26 cellule (`DR 39+`/tum, `DR 37-`/`DR 9-`/`DR 41-`/epi_cnv) va
  letto come un'osservazione puntuale, non come uno stato tumorale generalizzabile.

---

## 7. File sorgente (per riferimento/citazione)

```
05_drvi_tumoral_epi/tables/{scie,emt,gavish_tnbc}/drvi_tum_64_nomt/convergence_*.csv
05_drvi_tumoral_epi/tables/{scie,emt,gavish_tnbc}/drvi_epicnv_64_nomt/convergence_*.csv
05_drvi_tumoral_epi/tables/05_3_drvi_tum_64_nomt/{smi_leiden_matches,smi_celltypist_01_4_matches,factor_processes}_drvi_tum_64_nomt*.csv
05_drvi_tumoral_epi/tables/05_3_drvi_epicnv_64_nomt/{smi_leiden_matches,smi_cell_type_matches,smi_cnv_status_matches,smi_celltypist_01_4_matches,factor_processes}_drvi_epicnv_64_nomt*.csv
05_drvi_tumoral_epi/tables/05_2_visualization/clusters_vs_label_{tum,epicnv}.csv
05_drvi_tumoral_epi/tables/05_2_visualization/summary_{tum,epicnv}.csv
05_drvi_tumoral_epi/05_8_convergence/README.md
05_drvi_tumoral_epi/README.md
```
