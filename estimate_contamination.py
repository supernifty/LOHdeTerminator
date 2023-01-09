#!/usr/bin/env python
'''
  estimate % tumour from allele frequencies
'''

import argparse
import csv
import collections
import gzip
import logging
import sys

import cyvcf2
import numpy as np
import scipy.stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def estimate_percentile(values):
  '''
    ultra simple approach of just taking the 99th percentile
  '''
  distribution = [v for v in values if v > 0.1]
  if len(distribution) == 0:
    return {'tumour': 1.0}
  return {'tumour': np.percentile(np.array(distribution), 99)}

def generate_estimate(params):
  # now make our estimate based on our distributions
  l = np.linspace(0, 1, 100)
  #logging.debug('%s %s', l, params)
  cdf_signal = scipy.stats.norm.cdf(l, loc=params[0], scale=params[1]) * params[4]
  cdf_noise = scipy.stats.gamma.cdf(l, loc=params[2], a=params[3]) * (1-params[4])
  cdf = cdf_signal + cdf_noise # TODO can I just do this?
  predictions = [cdf[i+1] - cdf[i] for i in range(0, 99)]
  return predictions
 
def measure_error(values):
  def error_fn(params):
    #logging.debug('trying %s...', params)
    # make the values into a histogram
    targets, xs = np.histogram(values, bins=99, range=(0, 1), density=True)
    targets /= sum(targets)
    targets = np.array([max(target, 0.001) for target in targets])
    #logging.debug('targets %s', len(targets))
  
    # now make our estimate based on our distributions
    predictions = generate_estimate(params)
    
    # rmse
    error = np.sqrt(np.mean((predictions-targets)**2))
    # TODO log might be better
    #error = np.sqrt(np.mean((np.log(predictions)-np.log(targets))**2))

    logging.debug('trying %s: error %.2f', params, error)
    return error
  return error_fn

def estimate_model(values):
  '''
    model a sum of distributions
  '''
  # gamma + gaussian = four parameters
  # gamma_alpha, gamma_beta, gaussian_mean, gaussian_sd
  # start off with some reasonable estimates
  # todo some proportion of each as well
  #params = { 'gamma_alpha': 0.1, 'gamma_beta': 1.0, 'gaussian_mu': 0.4, 'gaussian_sd': 0.1, 'gamma_proportion': 0.5 }
  params = [0.4, 0.1, 0.1, 1.0, 0.5]
  bounds=[(0.0, np.inf), (0.0, np.inf), (0.01, 0.6), (0.001, np.inf), (0.0, 1.0)]
  minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds)
  result = scipy.optimize.basinhopping(measure_error(values), params, minimizer_kwargs=minimizer_kwargs)
  logging.info('%s has error %.2f', result.x, measure_error(values)(result.x))

  # measure the error and solve
  return {'tumour': result.x[2] * 2, 'distribution': generate_estimate(result.x)}

ESTIMATE = {
  'percentile': estimate_percentile,
  'model': estimate_model
}

def read_vcf(vcf_in, pass_only, dp_threshold, info_af):
  logging.info('reading vcf from stdin...')
  skipped_dp = skipped_pass = 0

  values = []

  for variant_count, variant in enumerate(vcf_in):
    # calculate vaf
    if len(variant.ALT) > 1:
      logging.warn('variant %i is multi-allelic', variant_count + 1)

    is_pass = variant.FILTER is None or variant.FILTER == 'alleleBias'
    if pass_only and not is_pass:
      skipped_pass += 1
      continue

    if hasattr(variant, "DP"):
      depth = variant.DP
    else:
      depth = variant.INFO["DP"]

    if depth < dp_threshold: # somatic + germline
      skipped_dp += 1
      continue

    if hasattr(variant, "AF"):
      value = variant.AF
    elif info_af:
      value = variant.INFO["AF"]
    else:
      ad = variant.format("AD")[sample_id]
      ref = ad[0]
      alt = ad[1]
      if ref + alt > 0:
        value = alt / (ref + alt)
      else:
        value = 0

    values.append(value)

  return values

def estimate(method, values, pass_only, dp_threshold, info_af):
  est = ESTIMATE[method](values)
  logging.info('done')
  return est

def main(vcf_in, pass_only, dp_threshold, info_af, plot, trials, prefix='Estimated tumour percentage (percentile):'):
  values = read_vcf(vcf_in, pass_only, dp_threshold, info_af)
  result = estimate('percentile', values, pass_only, dp_threshold, info_af)
  sys.stdout.write('{}\t{:.2f}\n'.format(prefix, result['tumour']))
  if trials > 0:
    results = []
    for _ in range(trials):
      result = estimate('model', values, pass_only, dp_threshold, info_af)
      results.append(result)
    sys.stdout.write('Estimated tumour percentage (model):\t{}\n'.format(' '.join(['{:.2f}'.format(result['tumour']) for result in results])))

  if plot:
    # draw histogram of normalised values and fitted distribution
    targets, xs = np.histogram(values, bins=99, range=(0, 1), density=True)
    targets /= sum(targets)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    l = np.linspace(0, 1, 99)
    ax.plot(l, targets, label='Observed')
    if trials > 0:
      ax.plot(l, results[0]['distribution'], label='Predicted')
    plt.tight_layout()
    plt.savefig(plot)

def open_file(fn, is_gzipped):
  if is_gzipped:
    return gzip.open(fn, 'rt')
  else:
    return open(fn, 'rt')

def no_chr(chrom):
  if chrom == 'MT':
    return 'M'
  # deals with chrUn_gl000220.1
  return chrom.split('.')[0].split('_', maxsplit=1)[-1].replace('chr', '').upper()

def get_value(header, col, row):
  return row[header.index(col)]

def maf_to_vcf(maf, sample, sample_col, chrom_col, pos_col, ref_col, alt_col, is_not_zipped):

  Variant = collections.namedtuple('Variant', 'CHROM POS REF ALT FILTER AF REF_COUNT ALT_COUNT DP')

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
    dp = ref_count + alt_count
    af = alt_count / dp

    yield Variant(chrom, pos, ref, (alt,), filtr, af, ref_count, alt_count, dp)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='estimate tumour %')
  parser.add_argument('--dp_threshold', required=False, default=50, help='use af in info field')
  parser.add_argument('--pass_only', action='store_true', help='only pass variants')
  parser.add_argument('--info_af', action='store_true', help='use af in info field')
  parser.add_argument('--trials', required=False, type=int, default=1, help='how many runs of model')
  parser.add_argument('--plot', required=False, help='use af in info field')
  parser.add_argument('--prefix', required=False, help='name')

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

  if args.maf_sample is not None:
    logging.info('opening maf file %s with sample %s...', args.maf_filename, args.maf_sample)
    vcf_in = maf_to_vcf(args.maf_filename, args.maf_sample, args.maf_sample_column, args.maf_chrom_column, args.maf_pos_column, args.maf_ref_column, args.maf_alt_column, args.vcf_not_zipped)
  else:
    logging.info('opening vcf from stdin...')
    vcf_in = cyvcf2.VCF('-')

  main(vcf_in, args.pass_only, args.dp_threshold, args.info_af, args.plot, args.trials, args.prefix)
