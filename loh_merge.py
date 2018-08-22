#!/usr/bin/env python
'''
  loh
'''

import argparse
import collections
import logging
import math
import sys

import cyvcf2
import intervaltree

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def write_merged(start, last, loh_status, stats, min_len, min_prop, bafs):
  if loh_status['loh'] and loh_status['accepts'] > 0:
    length = last[1] - start[1]
    prop = loh_status['accepts'] / (loh_status['accepts'] + loh_status['supports'] + loh_status['neutrals'])
    if length > min_len and prop > min_prop:
      sys.stdout.write("{chrom}\t{start}\t{end}\t{accept_pct:.1f}\t{accepts}\t{supports}\t{neutrals}\t{length}\n".format(chrom=start[0], start=start[1], end=last[1], accept_pct=100. * prop, accepts=loh_status['accepts'], supports=loh_status['supports'], neutrals=loh_status['neutrals'], length=length))
      stats['count'] += 1
      stats['length'] += length
      stats['count-{}'.format(start[0])] += 1
      stats['length-{}'.format(start[0])] += length
      bafs['calls'].append((start[1], last[1]))
      return True
  return False

def plot_bafs(filename, bafs):
  if len(bafs['x']) > 0:
    plt.figure(figsize=(24, 8))
    plt.scatter(bafs['x'], bafs['y'], c=bafs['c'], alpha=0.5)
    # add calls
    for call in bafs['calls']:
      plt.axvspan(call[0], call[1], facecolor='#00cc00', alpha=0.2)

    import matplotlib.patches as mpatches
    l1 = mpatches.Patch(color='green', label='LOH')
    l2 = mpatches.Patch(color='red', label='Not LOH')
    l3 = mpatches.Patch(color='blue', label='Favours LOH')
    l4 = mpatches.Patch(color='#c0c0c0', label='Neutral')
    plt.legend( [l1, l2, l3, l4],['LOH', 'Not LOH', 'Favours LOH', 'Neutral'], loc = 'lower center', ncol=5, labelspacing=0. )

    plt.savefig(filename)
    plt.close()

    bafs = {'x': [], 'y': [], 'c': []}
    logging.info('wrote %s', filename)

def calculate_segments(min_len, min_prop, noheader, plot):
  start = None # last potential end
  last = None # last potential end
  last_chrom = None
  loh_status = {'loh': False, 'accepts': 0, 'supports': 0, 'neutrals': 0}
  stats = collections.defaultdict(int)
  bafs = {'x': [], 'y': [], 'c': [], 'calls': []}

  if not noheader:
    sys.stdout.write("Chr\tStart\tEnd\tAccept_Pct\tAccepts\tSupports\tNeutrals\tLength\n")

  for line in sys.stdin:
    if line.startswith('chrom'): # header
      continue
    fields = line.strip('\n').split('\t')
    chrom = fields[0]
    pos = int(fields[1])
    baf = float(fields[3])
    status = fields[6]

    if plot is not None:
      bafs['x'].append(pos)
      bafs['y'].append(baf)
      if status.startswith('accept'):
        bafs['c'].append('green')
      elif status.startswith('reject'):
        bafs['c'].append('red')
      elif status.startswith('support'):
        bafs['c'].append('blue')
      else: # neutral
        bafs['c'].append('#c0c0c0')

    if status.startswith('reject') or last_chrom != chrom: # end of loh
      write_merged(start, last, loh_status, stats, min_len, min_prop, bafs)
      loh_status = {'loh': False, 'accepts': 0, 'supports': 0, 'neutrals': 0}
      start = None
      last = None

      if plot is not None and last_chrom != chrom:
        plot_bafs('{}.{}.png'.format(plot, last_chrom), bafs)
        bafs = {'x': [], 'y': [], 'c': [], 'calls': []}

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

  write_merged(start, last, loh_status, stats, min_len, min_prop, bafs)
  if plot is not None:
    plot_bafs('{}.{}.png'.format(plot, last_chrom), bafs)

  logging.info('done: %s', ', '.join(['{}: {}'.format(k, stats[k]) for k in stats]))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='LOH caller')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--min_len', default=1000, type=int, help='minimum segment size')
  parser.add_argument('--min_prop', type=float, default=0.1,  help='minimum accept proportion')
  parser.add_argument('--noheader', action='store_true',  help='no header for bed compatibility')
  parser.add_argument('--plot', required=False, help='plot results with image prefix')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  calculate_segments(args.min_len, args.min_prop, args.noheader, args.plot)
