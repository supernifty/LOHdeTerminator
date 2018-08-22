#!/usr/bin/env python
'''
  given tumour and germline vcf pairs, explore msi status
'''

import argparse
import cyvcf2
import logging
import math
import sys

PASSED='+'
FILTERED=''

# het in the germline
MIN_AF_HET_GL=0.3
MAX_AF_HET_GL=0.7

MIN_HET_HOM_DIFF=0.3 # het -> hom diff required (accept)
MIN_HOM_REF_DIFF=0.6 # hom -> hom diff required (support)

# het in the tumour
MIN_AF_HET_TUMOUR=0.3
MAX_AF_HET_TUMOUR=0.7

def calculate_af(variant, sample_id):
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

def main(tumour, germline, neutral_variants, filtered_variants, min_dp_germline, min_dp_tumour, min_af):
  logging.info('reading vcf from stdin...')
  vcf_in = cyvcf2.VCF('-')
  skipped_pass = count = 0

  tumour_id = vcf_in.samples.index(tumour)
  germline_id = vcf_in.samples.index(germline)

  stats = {'af_g_min': 1.0, 'af_g_max': 0.0, 'af_t_min': 1.0, 'af_t_max': 0.0, 'accept': 0, 'support': 0, 'reject': 0, 'accept': 0, 'neutral': 0, 'min_dp': 0, 'min_af': 0}

  sys.stdout.write('chrom\tpos\tg_af\tt_af\tg_dp\tt_af\tstatus\n')

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

    germline_ad_ref, germline_ad_alt = calculate_af(variant, germline_id)
    if germline_ad_ref + germline_ad_alt < min_dp_germline:
      stats['min_dp'] += 1
      continue
    else:
      germline_af = germline_ad_alt / (germline_ad_ref + germline_ad_alt)

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
    if MIN_AF_HET_GL < germline_af < MAX_AF_HET_GL and abs(germline_af - tumour_af) > MIN_HET_HOM_DIFF:
      write(stats, variant, germline_af, tumour_af, germline_ad_ref + germline_ad_alt, tumour_ad_ref + tumour_ad_alt, 'accept')
      stats['accept'] += 1

    # hom -> ref or ref -> hom
    elif abs(germline_af - tumour_af) > MIN_HOM_REF_DIFF:
      write(stats, variant, germline_af, tumour_af, germline_ad_ref + germline_ad_alt, tumour_ad_ref + tumour_ad_alt, 'support')

    # * -> het
    elif MIN_AF_HET_TUMOUR < tumour_af < MAX_AF_HET_TUMOUR:
      write(stats, variant, germline_af, tumour_af, germline_ad_ref + germline_ad_alt, tumour_ad_ref + tumour_ad_alt, 'reject')

    # het -> ambiguous, hom -> hom, ref -> ref
    elif neutral_variants:
      write(stats, variant, germline_af, tumour_af, germline_ad_ref + germline_ad_alt, tumour_ad_ref + tumour_ad_alt, 'neutral')

  logging.info('done. total %i skipped %i stats %s', count + 1, skipped_pass, ', '.join(['{}: {}'.format(k, stats[k]) for k in stats]))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='LOH caller')
  parser.add_argument('--tumour', required=True, help='tumour sample name')
  parser.add_argument('--germline', required=True, help='germline sample name')
  parser.add_argument('--neutral_variants', action="store_true", help='include neutral variants')
  parser.add_argument('--filtered_variants', action="store_true", help='include filtered variants')
  parser.add_argument('--min_dp_germline', type=int, default=1, help='min dp for germline')
  parser.add_argument('--min_dp_tumour', type=int, default=1, help='min dp for tumour')
  parser.add_argument('--min_af', type=float, default=0.0, help='min af for either germline or tumour')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.tumour, args.germline, args.neutral_variants, args.filtered_variants, args.min_dp_germline, args.min_dp_tumour, args.min_af)
