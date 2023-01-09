#!/usr/bin/env python
'''
  generate a list of loh evidence based on an input vcf
  output from this file is used as input to loh_merge
  - supported callers: strelka
'''

import argparse
import collections
import csv
import cyvcf2
import gzip
import logging
import math
import sys

import version

PASSED='+'
FILTERED=''

# het in the germline
MIN_AF_HET_GERMLINE=0.3
MAX_AF_HET_GERMLINE=0.7

MIN_AF_HET_TUMOUR=0.3
MAX_AF_HET_TUMOUR=0.7

MIN_HET_HOM_DIFF=0.3 # het -> hom diff required (accept)
MIN_HOM_REF_DIFF=0.6 # hom -> hom diff required (support)

def calculate_af(variant, sample_id):
  # mafs
  if hasattr(variant, 'REF_COUNT') and hasattr(variant, 'ALT_COUNT'):
    return variant.REF_COUNT, variant.ALT_COUNT

  try:
    ad_ref, ad_alt = variant.format("AD")[sample_id][:2]
  except KeyError:
    # strelka snvs
    try:
      tier1RefCounts = variant.format('{}U'.format(variant.REF))[sample_id] # e.g 12,12 tier1,tier2
      tier1AltCounts = variant.format('{}U'.format(variant.ALT[0]))[sample_id] # assume not multiallelic
      if len(variant.ALT) > 1:
        logging.warn('variant %i is multi-allelic', variant_count + 1)

      ad_ref = ad_alt = 0
      for idx, refCountList in enumerate(tier1RefCounts):
        ad_alt += int(tier1AltCounts[idx])
        ad_ref += int(refCountList)
    except KeyError: # strelka indels
      ad_alt = sum(variant.format('TIR')[sample_id])
      ad_ref = sum(variant.format('TAR')[sample_id])
 
  return ad_ref, ad_alt

def write(stats, variant, germline_af, tumour_af, germline_dp, tumour_dp, cat):
  if variant.FILTER is None:
    pass_marker = PASSED
  else:
    pass_marker = FILTERED
  sys.stdout.write('{}\t{}\t{:.2f}\t{:.2f}\t{}\t{}\t{}{}\n'.format(variant.CHROM, variant.POS, germline_af, tumour_af, germline_dp, tumour_dp, cat, pass_marker))
  stats[cat] += 1

def main(vcf_in, tumour, germline, neutral_variants, filtered_variants, min_dp_germline, min_dp_tumour, min_af, tumour_cellularity, min_het_hom, min_hom_hom, min_af_het_tumour_threshold, max_af_het_tumour_threshold, germline_vcf, min_germline_het, max_germline_het):
  # calculate tumour thresholds based on cellularity
  min_af_het_tumour= min_af_het_tumour_threshold * tumour_cellularity
  max_af_het_tumour= max_af_het_tumour_threshold * tumour_cellularity
  logging.info('tumour cutoffs are %.2f to %.2f with cellularity %.2f het diff %.2f hom diff %.2f min tumour het %.2f max tumour het %.2f', min_af_het_tumour, max_af_het_tumour, tumour_cellularity, min_het_hom, min_hom_hom, min_af_het_tumour_threshold, max_af_het_tumour_threshold)

  skipped_pass = count = 0
  sys.stdout.write('chrom\tpos\tg_af\tt_af\tg_dp\tt_af\tstatus\n')

  if tumour is not None: # none for mafs
    tumour_id = vcf_in.samples.index(tumour)
  else:
    tumour_id = None
  if germline is not None:
    germline_id = vcf_in.samples.index(germline)

  stats = {'af_g_min': 1.0, 'af_g_max': 0.0, 'af_t_min': 1.0, 'af_t_max': 0.0, 'accept': 0, 'support': 0, 'reject': 0, 'accept': 0, 'neutral': 0, 'min_dp': 0, 'min_af': 0}

  seen = set()
  for count, variant in enumerate(vcf_in):
    if (count + 1) % 10000 == 0:
      logging.debug('processed %i variants. skipped %i', count + 1, skipped_pass)

    if variant.FILTER is None:
      pass_marker = '+'
    else:
      if not filtered_variants:
        skipped_pass += 1
        continue
      else:
        pass_marker = ''

    seen.add((variant.CHROM, variant.POS))

    if germline is not None:
      germline_ad_ref, germline_ad_alt = calculate_af(variant, germline_id)
      if germline_ad_ref + germline_ad_alt < min_dp_germline:
        stats['min_dp'] += 1
        continue
      else:
        germline_af = germline_ad_alt / (germline_ad_ref + germline_ad_alt)
    else: # no germline, assume het
      germline_af = 0.5
      germline_ad_ref, germline_ad_alt = 100, 100
  
    stats['af_g_min'] = min(stats['af_g_min'], germline_af)
    stats['af_g_max'] = max(stats['af_g_max'], germline_af)

    tumour_ad_ref, tumour_ad_alt = calculate_af(variant, tumour_id)
    if tumour_ad_ref + tumour_ad_alt < min_dp_tumour:
      stats['min_dp'] += 1
      continue 
    else:
      tumour_af = tumour_ad_alt / (tumour_ad_ref + tumour_ad_alt)

    if germline_af < min_af and tumour_af < min_af:
      stats['min_af'] += 1
      continue

    stats['af_t_min'] = min(stats['af_t_min'], tumour_af)
    stats['af_t_max'] = max(stats['af_t_max'], tumour_af)

    # het -> hom
    if min_germline_het < germline_af < max_germline_het and abs(germline_af - tumour_af / tumour_cellularity) > min_het_hom:
      write(stats, variant, germline_af, tumour_af, germline_ad_ref + germline_ad_alt, tumour_ad_ref + tumour_ad_alt, 'accept')
      stats['accept'] += 1

    # hom -> ref or ref -> hom
    elif abs(germline_af - tumour_af / tumour_cellularity) > min_hom_hom:
      write(stats, variant, germline_af, tumour_af, germline_ad_ref + germline_ad_alt, tumour_ad_ref + tumour_ad_alt, 'support')

    # * -> het
    elif min_af_het_tumour < tumour_af < max_af_het_tumour:
      write(stats, variant, germline_af, tumour_af, germline_ad_ref + germline_ad_alt, tumour_ad_ref + tumour_ad_alt, 'reject')

    # het -> ambiguous, hom -> hom, ref -> ref
    elif neutral_variants:
      write(stats, variant, germline_af, tumour_af, germline_ad_ref + germline_ad_alt, tumour_ad_ref + tumour_ad_alt, 'neutral')

  logging.info('done with tumour. total %i skipped %i stats %s', count + 1, skipped_pass, ', '.join(['{}: {}'.format(k, stats[k]) for k in stats]))

  if germline_vcf is not None:
    logging.info('processing %s...', germline_vcf)
    vcf_in = cyvcf2.VCF(germline_vcf)
    for count, variant in enumerate(vcf_in):
      if (variant.CHROM, variant.POS) in seen:
        logging.debug('skipping already seen variant: %s:%s', variant.CHROM, variant.POS)
        continue
      
      if variant.FILTER is None:
        pass_marker = '+'
      else:
        if not filtered_variants:
          skipped_pass += 1
          continue
        else:
          pass_marker = ''

      germline_ad_ref, germline_ad_alt = calculate_af(variant, 0)
      if germline_ad_ref + germline_ad_alt < min_dp_germline:
        stats['min_dp'] += 1
        continue
      else:
        germline_af = germline_ad_alt / (germline_ad_ref + germline_ad_alt)

      if min_germline_het < germline_af < max_germline_het:
        write(stats, variant, germline_af, germline_af, germline_ad_ref + germline_ad_alt, -1, 'reject')
    logging.info('processing %s: done', germline_vcf)
 
def no_chr(chrom):
  if chrom == 'MT':
    return 'M'
  # deals with chrUn_gl000220.1
  return chrom.split('.')[0].split('_', maxsplit=1)[-1].replace('chr', '').upper()

def get_value(header, col, row):
  return row[header.index(col)]

def open_file(fn, is_gzipped):
  if is_gzipped:
    return gzip.open(fn, 'rt')
  else:
    return open(fn, 'rt')

def maf_to_vcf(maf, sample, sample_col, chrom_col, pos_col, ref_col, alt_col, is_not_zipped, germline_vcf):

  Variant = collections.namedtuple('Variant', 'CHROM POS REF ALT FILTER AF REF_COUNT ALT_COUNT')

  # enumeration a maf into a variant
  header = None
  for line, row in enumerate(csv.reader(open_file(maf, not is_not_zipped), delimiter='\t')):
    if line % 1000 == 0:
      logging.debug('processed %i lines of %s...', line, maf)

    if row[0].startswith('#'):
      continue
    if header is None:
      header = row
      continue

    #Hugo_Symbol     Entrez_Gene_Id  Center  NCBI_Build      Chromosome      Start_Position  End_Position    Strand  Variant_Classification  Variant_Type    Reference_Allele        Tumor_Seq_Allele1       Tumor_Seq_Allele2       dbSNP_RS        dbSNP_Val_Status        Tumor_Sample_Barcode    Matched_Norm_Sample_Barcode     Match_Norm_Seq_Allele1  Match_Norm_Seq_Allele2  Tumor_Validation_Allele1        Tumor_Validation_Allele2        Match_Norm_Validation_Allele1   Match_Norm_Validation_Allele2   Verification_Status     Validation_Status       Mutation_Status Sequencing_Phase        Sequence_Source Validation_Method       Score   BAM_File        Sequencer       Tumor_Sample_UUID       Matched_Norm_Sample_UUID        HGVSc   HGVSp   HGVSp_Short     Transcript_ID   Exon_Number     t_depth t_ref_count     t_alt_count     n_depth n_ref_count     n_alt_count     all_effects     Allele  Gene    Feature Feature_type    One_Consequence Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids Codons  Existing_variation      ALLELE_NUM      DISTANCE        TRANSCRIPT_STRAND       SYMBOL  SYMBOL_SOURCE   HGNC_ID BIOTYPE CANONICAL       CCDS    ENSP    SWISSPROT       TREMBL  UNIPARC RefSeq  SIFT    PolyPhen        EXON    INTRON  DOMAINS GMAF    AFR_MAF AMR_MAF ASN_MAF EAS_MAF EUR_MAF SAS_MAF AA_MAF  EA_MAF  CLIN_SIG        SOMATIC PUBMED  MOTIF_NAME      MOTIF_POS       HIGH_INF_POS    MOTIF_SCORE_CHANGE      IMPACT  PICK    VARIANT_CLASS   TSL     HGVS_OFFSET     PHENO   MINIMISED       ExAC_AF ExAC_AF_Adj     ExAC_AF_AFR     ExAC_AF_AMR     ExAC_AF_EAS     ExAC_AF_FIN     ExAC_AF_NFE     ExAC_AF_OTH     ExAC_AF_SAS     GENE_PHENO      FILTER  CONTEXT src_vcf_id      tumor_bam_uuid  normal_bam_uuid case_id GDC_FILTER      COSMIC  MC3_Overlap     GDC_Validation_Status

    row_sample = get_value(header, sample_col, row)
    if sample is not None and row_sample != sample:
      continue

    chrom = no_chr(get_value(header, chrom_col, row))
    pos = int(get_value(header, pos_col, row))
    ref = get_value(header, ref_col, row)
    if ref == '-':
      pos += 1 # fix for TCGA mafs
    ref = ref.replace('-', '')
    alt = get_value(header, alt_col, row).replace('-', '')
    filtr = get_value(header, 'FILTER', row)
    if filtr == 'PASS':
      filtr = None
    ref_count = float(get_value(header, 't_ref_count', row))
    alt_count = float(get_value(header, 't_alt_count', row))
    af = alt_count / (ref_count + alt_count)

    yield Variant(chrom, pos, ref, (alt,), filtr, af, ref_count, alt_count)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='LOH caller')
  parser.add_argument('--version', action='version', version=version.PROGRAM_VERSION)
  parser.add_argument('--tumour', required=True, help='tumour sample name')
  parser.add_argument('--tumour_cellularity', required=False, type=float, default=1.0, help='tumour cellularity')
  parser.add_argument('--germline', required=False, help='germline sample name')
  parser.add_argument('--neutral_variants', action="store_true", help='include neutral variants')
  parser.add_argument('--filtered_variants', action="store_true", help='include filtered variants')
  parser.add_argument('--min_dp_germline', type=int, default=1, help='min dp for germline')
  parser.add_argument('--min_dp_tumour', type=int, default=1, help='min dp for tumour')
  parser.add_argument('--min_af', type=float, default=0.0, help='min af for either germline or tumour')
  parser.add_argument('--min_het_hom', type=float, default=MIN_HET_HOM_DIFF, help='change to be considered het to hom default 0.3')
  parser.add_argument('--min_hom_hom', type=float, default=MIN_HOM_REF_DIFF, help='change to be considered hom to hom default 0.7')
  parser.add_argument('--min_tumour_het', type=float, default=MIN_AF_HET_TUMOUR, help='min tumour af to be het')
  parser.add_argument('--max_tumour_het', type=float, default=MAX_AF_HET_TUMOUR, help='max tumour af to be het')
  parser.add_argument('--min_germline_het', type=float, default=MIN_AF_HET_GERMLINE, help='min germline af to be het')
  parser.add_argument('--max_germline_het', type=float, default=MAX_AF_HET_GERMLINE, help='max germline af to be het')
  parser.add_argument('--germline_vcf', required=False, help='additional germline vcf for germline only hets')

  parser.add_argument('--maf_filename', required=False, help='vcf is actually a maf with this filename')
  parser.add_argument('--maf_sample', required=False, help='vcf is actually a maf with this sample of interest')
  parser.add_argument('--maf_sample_column', required=False, default='Tumor_Sample_Barcode', help='maf chrom column name')
  parser.add_argument('--maf_chrom_column', required=False, default='Chromosome', help='maf chrom column name')
  parser.add_argument('--maf_pos_column', required=False, default='Start_Position', help='maf pos column name')
  parser.add_argument('--maf_ref_column', required=False, default='Reference_Allele', help='maf ref column name')
  parser.add_argument('--maf_alt_column', required=False, default='Tumor_Seq_Allele2', help='maf alt column name')
  parser.add_argument('--vcf_not_zipped', action='store_true', help='do not try to unzip (only matters for maf)')

  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  if args.tumour_cellularity > 1.0:
    args.tumour_cellularity /= 100

  if args.maf_sample:
    vcf_in = maf_to_vcf(args.maf_filename, args.maf_sample, args.maf_sample_column, args.maf_chrom_column, args.maf_pos_column, args.maf_ref_column, args.maf_alt_column, args.vcf_not_zipped, args.germline_vcf)
    args.tumour = None
  else:
    logging.info('reading vcf from stdin...')
    vcf_in = cyvcf2.VCF('-')

  main(vcf_in, args.tumour, args.germline, args.neutral_variants, args.filtered_variants, args.min_dp_germline, args.min_dp_tumour, args.min_af, args.tumour_cellularity, args.min_het_hom, args.min_hom_hom, args.min_tumour_het, args.max_tumour_het, args.germline_vcf, min_germline_het=args.min_germline_het, max_germline_het=args.max_germline_het)
