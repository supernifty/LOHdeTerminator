# LOH caller

## Installation
```
python3 -m venv venv
. ./venv/bin/activate
pip install -r requirements.txt
```

## Usage
```
python loh_determinator/loh.py --germline NORMAL --tumour TUMOR --filtered_variants --min_dp_germline 10 --min_dp_tumour 20 --neutral --germline_vcf strelka.vcf < strelka.vcf | sort -k1,1 -k2,2n > loh.tsv
python loh_determinator/loh_merge.py --noheader --min_len 1000 --min_prop 0.1 < loh.tsv > loh.merged.bed
python loh_determinator/combine_loh.py TODO
```

* the sort is required if you use --germline_vcf

## Method

|tumour &downarrow; normal &rightarrow;| ref | het | hom |
|-|-|-|-|
|ref |      neutral | accept | support |
|het |      reject  | reject | reject |
|hom |      support | accept | neutral |

## Output

### loh output - individual variants
* chromosome
* position
* germline af
* tumour af
* germline dp
* tumour dp
* assessment of LOH potential

### loh_merged output - predicted areas of LOH
* chrom: chromosome
* start: start of LOH region
* end: end of LOH region
* accept_pct: proportion of variants in the region in the "accept" category compared to accept + support + neutral
* accepts: number of accept variants
* supports: number of supporting variants
* neutrals: number of neutral variants
* length: length of region

### combine_loh output
* sample
* gene
* accept
* length
* accept_per_mb

## Authors
* Romy Walker
* Peter Georgeson
