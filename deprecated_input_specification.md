### allele frequency file

Previously an sparse table was used for allele frequency input. This file is still supported with the --data2 option. See here the specification of this deprecated file format. 

A TAB separated table with columns in four blocks. First block gives the annotatio of the associated mutation (column 1-18). The next block of N columns (with n samples) specify the allele frequencies. Thereby the column header is named SAMPLE:AF. The next block of N columns specify the read depth in each sample for taht mutation. Thereby the column header is named SAMPLE:DP.  The last block of n columns provides the quality value for the mutation calling. Thereby the column header is named SAMPLE:PQ. 

| Col | HEADER               | EXAMPLE VALUE      | DESCRIPTION                  |
| --- |                  --: |  :--               | :--                          |
| $1  |                CHROM |  NC_045512.2       | chromosome name              |
| $2  |                  POS |  3002              | position                     |
| $3  |                  REF |  G                 | reference base               |
| $4  |                  ALT |  T                 | alternative base             |
| $5  |             ANN.GENE |  ORF1ab            | annotation gene              |
| $6  |        ANN.FEATUREID |  pp1ab_GU280_gp01  | annotation ID                |
| $7  |           ANN.EFFECT |  stop_gained       | annotation effect            |
| $8  |               ANN.AA |  p.E913stop        | AA mutation                  |
| $9  |           ANN.GENE.1 |  ORF1ab            | annotation level 2           |
| $10 |      ANN.FEATUREID.1 |  pp1a_ORF1a        | annotation level 2           |
| $11 |             ANN.AA.1 |  p.E913stop        | annotation level 2           |
| $12 |           ANN.GENE.2 |  ORF1ab            | annotation level 3           |
| $13 |      ANN.FEATUREID.2 |  nsp3              | annotation level 3           |
| $14 |             ANN.AA.2 |  p.E95stop         | annotation level 3           |
| $15 |           ANN.GENE.3 |  ""                | annotation level 4           |
| $16 |      ANN.FEATUREID.3 |  ""                | annotation level 4           |
| $17 |             ANN.AA.3 |  ""                | annotation level 4           |
| $18 |                  EFF |  Gag/Tag           | phased mutation              |
| $19 |      SAMPLE_NAME:AF  |  SAMPLE1:AF        | sample name and allele freq  |
| +N  |                      |                    | N more sample and AF thereof |
| $25 |      SAMPLE_NAME:DP  |  SAMPLE1:DP        | sample name and read depth   |
| +N  |                      |                    | N more sample and DP thereof |
| $32 |      SAMPLE_NAME:PQ  |  SAMPLE1:PQ        | sample name and quality      |
| +N  |                      |                    | N more sample and PQ thereof |
