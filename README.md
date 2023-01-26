# VaQuERo
SARS-CoV-2 (Va)riant (Qu)antification in s(E)wage, designed for (Ro)bustness

## Synopsis

Rscript scripts/VaQuERo_v2.r > vaquero.log

## Description
Takes allele frequencies in table formate and associated meta data and quantifies different virus variants as defined in provided marker mutation file. The workflow first calls 'detected' variants and subsequently quantifies the abundance of the detected variants based on the observed non-zero frequencies using a SIMPLEX regression.

## Usage 


scripts/VaQuERo_v2.r [options]

## Options
	--dir=CHARACTER
		Directory to write results [default: ExampleOutput]

	--country=CHARACTER
		Name of country used to produce map [default: Austria]

	--bbsouth=CHARACTER
		Bounding box most south point [default: 46.38]

	--bbnorth=CHARACTER
		Bounding box most norther point [default: 49.01]

	--bbwest=CHARACTER
		Bounding box most western point [default: 9.53]

	--bbeast=CHARACTER
		Bounding box most easter point [default: 17.15]

	--metadata=CHARACTER
		Path to meta data input file [default: data/metaDataSub.tsv]

	--marker=CHARACTER
		Path to marker mutation input file [default: resources/mutations_list.csv]

	--smarker=CHARACTER
		Path to special mutation input file [default: resources/mutations_special.csv]

	--pmarker=CHARACTER
                Path to problematic mutations input file, which will be omitted throughout [default: resources/mutations_problematic_all.csv]

	--data=CHARACTER
		Path to data input file [default: data/mutationDataSub.tsv]
		
	--data2=CHARACTER
		Deprecated, for backwards compatibility only. Path to data input file in deprected old sparse table file format. Only effective together with --inputformat=sparse. [default: data/mutationDataSub.tsv]

       --inputformat=CHARACTER
                For backwards compatibility, set to 'sparse' to use input from --data2. Description of deprecated sparse table file format see file deprecated_input_specification.md. [default: tidy]

	--plotwidth=CHARACTER
		Base size of plot width [default: 8]

	--plotheight=CHARACTER
		Base size of plot height [default: 4.5]

	--ninconsens=CHARACTER
		Minimal fraction of genome covered by reads to be considered (0-1) [default: 0.4]

	--zero=DOUBLE
		Minimal allele frequency to be considered [default: 0.02]

	--depth=CHARACTER
		Minimal depth at mutation locus to be considered [default: 75]

	--recent=CHARACTER
		How old (in days) most recent sample might be to be still considered in overview maps [default: 99]

	--plottp=CHARACTER
		Produce timecourse plots only if more than this timepoints are available [default: 3]

	--minuniqmark=CHARACTER
		Minimal absolute number of uniq markers that variant is considered detected [default: 3]

	--minuniqmarkfrac=CHARACTER
		Minimal fraction of uniq markers that variant is considered detected [default: 0]

	--minqmark=CHARACTER
		Minimal absolute number of markers that variant is considered detected [default: 3]

	--minmarkfrac=CHARACTER
		Minimal fraction of markers that variant is considered detected [default: 0]

	--smoothingsamples=CHARACTER
		Number of previous timepoints use for smoothing [default: 1]

	--smoothingtime=CHARACTER
		Previous timepoints for smoothing are ignored if more days than this days apart [default: 8]

	--voi=CHARACTER
		List of variants which should be plotted in more detail. List separated by semicolon [default: B.1.1.7;B.1.617.2;P.1;B.1.351]

	--highlight=CHARACTER
		List of variants which should be plotted at the bottom axis. List separated by semicolon [default: B.1.1.7;B.1.617.2]

	-h, --help
		Show help message and exit.
		
	--debug=LOGIC
		Set to TRUE to use provided input example file [default: FALSE]. 
		


## Input

### allele frequency file

A TAB separated tidy table, specifying for each mutation and each sample its allele frequency and sequencing depth. The helper script vcf2tsv_long.py assists to create this file directly from a list of vcf files.

| Col | HEADER               | EXAMPLE VALUE      | DESCRIPTION                  |
| --- |                  --: |  :--               | :--                          |
| $1  |             SAMPLEID |  SAMPLE1           | sample name                  |
| $2  |                CHROM |  NC_045512.2       | chromosome name              |
| $3  |                  POS |  3002              | position                     |
| $4  |                  REF |  G                 | reference base               |
| $5  |                  ALT |  T                 | alternative base             |
| $6  |                 GENE |  ORF1ab            | annotation gene              |
| $7  |                   AA |  E913stop          | AA substitution              |
| $8  |                   AF |  0.997152          | Allele frequency             |
| $9  |                   DP |  16855             | Sequencing depth             |
| $10 |                   PQ |  49314             | positional quality           |

A previously used sparse table file format is still supported for backwards compatibility. See --data2 and --inputformat for details how to invoke it. See the file deprecated_input_specification.md for a specification of this file format.

#### vcf2tsv_long.py

`vcf2tsv_long.py` needs the package `pysam` to be installed (eg. via `pip` or `conda`).

It takes a vcf file, a directory with vcfs, or a list with paths to vcf files as an input and creates a tidy table fit as input for VaQuERo.

You can test it with the vcf files in `data/vcf2tsv` in the following ways on files with and without SNPEFF annotations, merged samples, all vcf files in a directory or a list of files:

```
python3 scripts/vcf2tsv_long.py -i data/vcf2tsv/CoV_29633_S117400.vcf.gz | less
python3 scripts/vcf2tsv_long.py -i data/vcf2tsv/CoV_29634_S117319_ann.vcf | less
python3 scripts/vcf2tsv_long.py -i data/vcf2tsv/merged_samples.vcf.gz | less
python3 scripts/vcf2tsv_long.py -i data/vcf2tsv/ | less
python3 scripts/vcf2tsv_long.py -i data/vcf2tsv/vcf_files.txt | less
```

To redirect output to a new file, append to an existing file and filter with a certain minimal allele frequency (default 0.01):

```
python3 scripts/vcf2tsv_long.py -i data/vcf2tsv/merged_samples.vcf.gz -m 0.1 -o test_input.tsv.gz
python3 scripts/vcf2tsv_long.py -i data/vcf2tsv/CoV_29634_S117319_ann.vcf --append -m 0.1 -o test_input.tsv.gz
```



### meta data file

The meta data file connects each sample as defined in the mutation file with the sampling location and time. Each samples included in allele frequency file must be specified here.

| Col | HEADER                  |  EXAMPLE VALUE   | DESCRIPTION                 |
| --- |                     --: |  :--             | :--                         |
| $1  |                 BSF_run |  BSF_0895        | Sequencing batch            |
| $2  |         BSF_sample_name |  SAMPLE1         | Seq. sample name            |
| $3  |          BSF_start_date |  2021-01-18      | Seq. date                   |
| $4  |              LocationID |  ATTP_10-Krezlin | Sample loction ID           |
| $5  |            LocationName |  Krezlin         | Sample location name        |
| $6  |          N_in_Consensus |  125             | Nr. Of N in consensus seq   |
| $7  |              RNA_ID_int |  SAMPLE1         | Sample name                 |
| $8  |  additional_information |  Ct = xx.2       | any add. info               |
| $9  |             adress_town |  Krezlin         | Sampling Town               |
| $10 |        connected_people |  1234567         | Nr. of conncected people    |
| $11 |             dcpLatitude |  48.50452605     | latitude                    |
| $12 |            dcpLongitude |  11.70669547     | lingitude                   |
| $13 |       include_in_report |  TRUE            | should be reported (F|T)    |
| $14 |         report_category |  fictional       | included in which report    |
| $15 |             sample_date |  2020-12-20      | sample date                 |
| $16 |                  status |  passed_qc       | status (passed|failed)      |

### marker mutation definition

The marker mutation file links the mutation expected to be seen in a variant (sensitivity) for all variants. Table is comma separated. If one mutation is sensitive for more than one variants, a list of variants (semicolon separated) can be provided in columns 1. Same format applies to special mutations which are used for plotting and problematic sites which are ignored during analysis.


| Col | HEADER         |  EXAMPLE VALUE   | DESCRIPTION                         |
| --- |            --: |  :--             | :--                                 |
| $1  |       Variants |  AV.1;B.1.1.318  | List of Variants (seperated by ";") | 
| $2  |      Chromosom |  NC_045512.2     | Chromosome name                     |
| $3  |        Postion |  21990           | Chromosome position                 |
| $4  |            REF |  TTTA            | Reference base(s)                   |
| $5  |            ALT |  T               | Alternative base(s)                 |
| $6  |           Gene |  S               | Gene                                |
| $7  |  Sensitivities |  1;0.87;0.96     | Sensitivity in each of the variants |
| $8  |             AA |  S:Y144del       | Mutation in AA nomenclature         |
| $9  |            NUC |  TTTA21990T      | Mutation in nucc nomenclature       |



## Output

### globalFittedData.csv

Table holding the deduced variant frequencies.

### globalFullData.csv

Table holding the observed allele frequencies and the deduced variant frequencies.

### globalSpecialmutData.csv

Table holding the observed allele frequencies and the allels specified in special mutation file.

### summary.csv

List of all plots produced.

### Figure directory

Timecourse plots and map represention of results, to be included in a report.

## Install

The following R packages must be pre-installed: 

- tidyr
- ggplot2
- reshape2
- dplyr
- data.table
- gamlss
- ggmap
- tmaptools
- ggrepel
- scales
- betareg
- ggspatial
- sf
- rnaturalearth
- rnaturalearthdata
- optparse
- stringr

For convenience a Rscupt for package installation is provided (Rscript R_package_dependency_install.r)
