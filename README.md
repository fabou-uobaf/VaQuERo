# VaQuERo
SARS-CoV-2 (Va)riant (Qu)antification in s(E)wage, designed for (Ro)bustness

## Synopsis

Takes allele frequencies in table formate and associated meta data and quantifies different virus variants as defined in provided marker mutation file. The workflow first calls 'detected' variants and subsequently quantifies the abundance of the detected variants based on the observed non-zero frequencies using a SIMPLEX regression.

## Input

## mutation file

| Col | HEADER               | EXAMPLE VALUE      |
| --- |                  --: |  :--               |
| $1  |                CHROM |  NC_045512.2       | 
| $2  |                  POS |  3002              | 
| $3  |                  REF |  G                 | 
| $4  |                  ALT |  T                 | 
| $5  |             ANN.GENE |  ORF1ab            | 
| $6  |        ANN.FEATUREID |  pp1ab_GU280_gp01  | 
| $7  |           ANN.EFFECT |  stop_gained       | 
| $8  |               ANN.AA |  p.E913stop        | 
| $9  |           ANN.GENE.1 |  ORF1ab            | 
| $10 |      ANN.FEATUREID.1 |  pp1a_ORF1a        | 
| $11 |             ANN.AA.1 |  p.E913stop        | 
| $12 |           ANN.GENE.2 |  ORF1ab            | 
| $13 |      ANN.FEATUREID.2 |  nsp3              | 
| $14 |             ANN.AA.2 |  p.E95stop         |
| $15 |           ANN.GENE.3 |  ""                | 
| $16 |      ANN.FEATUREID.3 |  ""                | 
| $17 |             ANN.AA.3 |  ""                | 
| $18 |                  EFF |  Gag/Tag           | 
| $19 |   CoV_2212_S74699:AF |  .                 | 


## Output
