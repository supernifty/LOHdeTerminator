# LOH caller

## Installation
```
python3 -m venv venv
. ./venv/bin/activate
pip install -r requirements.txt
```

## Usage
```
python loh.py --germline NORMAL --tumour TUMOR --filtered_variants --min_dp_germline 10 --min_dp_tumour 20 --neutral < strelka.vcf > loh.tsv
python loh_merge.py --noheader --min_len 1000 --min_prop 0.1 < loh.tsv > loh.merged.bed
```

## Method

|tumour &downarrow; normal &rightarrow;| ref | het | hom |
|-|-|-|-|
|ref |      neutral | accept | support |
|het |      reject  | reject | reject |
|hom |      support | accept | neutral |

## Output

* chrom: chromosome
* start: start of LOH region
* end: end of LOH region
* accept_pct: proportion of variants in the region in the "accept" category compared to accept + support + neutral
* accepts: number of accept variants
* supports: number of supporting variants
* neutrals: number of neutral variants
* length: length of region

## Authors
* Romy Walker
* Peter Georgeson
