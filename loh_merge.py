#!/usr/bin/env python
'''
  loh
'''

import argparse
import collections
import cyvcf2
import logging
import math
import sys

def write_merged(start, last, loh_status, stats, min_len, min_prop):
  if loh_status['loh'] and loh_status['accepts'] > 0:
    length = last[1] - start[1]
    prop = loh_status['accepts'] / (loh_status['accepts'] + loh_status['supports'] + loh_status['neutrals'])
    if length > min_len and prop > min_prop:
      sys.stdout.write("{chrom}\t{start}\t{end}\t{accept_pct:.1f}\t{accepts}\t{supports}\t{neutrals}\t{length}\n".format(chrom=start[0], start=start[1], end=last[1], accept_pct=100. * prop, accepts=loh_status['accepts'], supports=loh_status['supports'], neutrals=loh_status['neutrals'], length=length))
      stats['count'] += 1
      stats['length'] += length
      stats['count-{}'.format(start[0])] += 1
      stats['length-{}'.format(start[0])] += length
      return True
  return False

def calculate_segments(min_len, min_prop, noheader):
  start = None # last potential end
  last = None # last potential end
  last_chrom = None
  loh_status = {'loh': False, 'accepts': 0, 'supports': 0, 'neutrals': 0}
  stats = collections.defaultdict(int)

  if not noheader:
    sys.stdout.write("Chr\tStart\tEnd\tAccept_Pct\tAccepts\tSupports\tNeutrals\tLength\n")

  for line in sys.stdin:
    if line.startswith('chrom'): # header
      continue
    fields = line.strip('\n').split('\t')
    chrom = fields[0]
    pos = int(fields[1])
    status = fields[6]

    if status.startswith('reject') or last_chrom != chrom: # end of loh
      write_merged(start, last, loh_status, stats, min_len, min_prop)
      loh_status = {'loh': False, 'accepts': 0, 'supports': 0, 'neutrals': 0}
      start = None
      last = None
    else: # potential loh region
      if start is None:
        start = (chrom, pos)
      loh_status['loh'] = True
      if status.startswith('accept'):
        loh_status['accepts'] += 1
      elif status.startswith('support'):
        loh_status['supports'] += 1
      elif status.startswith('neutral'):
        loh_status['neutrals'] += 1

      last = (chrom, pos)

    last_chrom = chrom

  write_merged(start, last, loh_status, stats, min_len, min_prop)
  logging.info('done: %s', ', '.join(['{}: {}'.format(k, stats[k]) for k in stats]))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='LOH caller')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--min_len', default=1000, type=int, help='minimum segment size')
  parser.add_argument('--min_prop', type=float, default=0.1,  help='minimum accept proportion')
  parser.add_argument('--noheader', action='store_true',  help='no header for bed compatibility')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  calculate_segments(args.min_len, args.min_prop, args.noheader)
