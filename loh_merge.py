#!/usr/bin/env python
'''
  plots areas that are proposed to contain loh
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
matplotlib.rcParams.update({'font.size': 18})

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

def plot_bafs(filename, bafs, title, start, finish, gene_starts, gene_finishes, gene_names, annotate):
  if len(bafs['x']) > 0:
    x_mb = [x / 1000000 for x in bafs['x']]
    plt.figure(figsize=(24, 8))
    if start is None:
      plt.scatter(x_mb, bafs['y'], c=bafs['c'], alpha=0.5)
      start = 0
      finish = max(bafs['x'])
      if annotate:
        for x, y, l in zip(x_mb, bafs['y'], bafs['l']):
          plt.annotate(l, x, y)
    else:
      xycs = [xyc for xyc in zip(x_mb, bafs['y'], bafs['c']) if start / 1000000 <= xyc[0] <= finish / 1000000]
      if len(xycs) > 0:
        plt.scatter([xyc[0] for xyc in xycs], [xyc[1] for xyc in xycs], c=[xyc[2] for xyc in xycs], alpha=0.5)
        if annotate:
          for x, y, l in zip(x_mb, bafs['y'], bafs['l']):
            plt.annotate(s=l, xy=(x, y))

    plt.title('Genomic regions containing evidence for LOH across {}'.format(title))
    plt.ylabel('Tumour Allele Frequency')
    plt.xlabel('Genomic Position (MB)')
    plt.gca().set_ylim([0.0, 1.0])

    # add calls
    for call in bafs['calls']:
      if start is None:
        plt.axvspan(call[0] / 1000000, call[1] / 1000000, facecolor='#00cc00', alpha=0.2)
      else:
        if call[0] <= finish and call[1] >= start:
          call_s = max(start, call[0])
          call_f = min(finish, call[1])
          plt.axvspan(call_s / 1000000, call_f / 1000000, facecolor='#00cc00', alpha=0.2)

    # highlight gene
    if gene_starts is not None:
      for gene_start, gene_finish, gene_name in zip(gene_starts, gene_finishes, gene_names):
        logging.debug('annotating %s at %i %i in %i %i', gene_name, gene_start, gene_finish, start, finish)
        if (gene_finish - gene_start) < (finish - start) * 0.005:
          gene_finish = gene_start + (finish - start) * 0.005
        logging.debug('annotating %s at %i %i', gene_name, gene_start, gene_finish)
        plt.axvspan(gene_start / 1000000, gene_finish / 1000000, ymin=0.16, ymax=0.21, facecolor='#000000', alpha=0.8)
        plt.text(gene_start / 1000000, 0.11, gene_name, horizontalalignment='left', verticalalignment='bottom', bbox={'facecolor': 'yellow', 'alpha': 0.5})

    import matplotlib.patches as mpatches
    l1 = mpatches.Patch(color='green', label='LOH')
    l2 = mpatches.Patch(color='red', label='Not LOH')
    l3 = mpatches.Patch(color='blue', label='Favours LOH')
    l4 = mpatches.Patch(color='#c0c0c0', label='Neutral')
    plt.legend( [l1, l2, l3, l4],['LOH', 'Not LOH', 'Favours LOH', 'Neutral'], loc = 'lower right', ncol=5, labelspacing=0. )

    plt.savefig(filename)
    plt.close()

    bafs = {'x': [], 'y': [], 'c': []}
    logging.info('wrote %s', filename)

def check_plot(plot, last_chrom, chrom, regions, region_names, region_padding, bafs, plot_chromosomes, annotate, plot_custom, sample):
  # additional regions
  if regions is not None:
    gene_starts = []
    gene_finishes = []
    gene_names = []
    for name, region, padding in zip(region_names, regions, region_padding):
      region_chrom, start_finish = region.split(':')
      if region_chrom == last_chrom:
        start, finish = [int(x) for x in start_finish.split('-')]
        if ',' in padding:
          before, after = [int(x) for x in padding.split(',')]
          plot_start = max(0, start - before)
          plot_finish = finish + after
        else:
          plot_start = max(0, start - int(padding))
          plot_finish = finish + int(padding)
        if sample is not None:
          plot_bafs('{}.{}.png'.format(plot, name), bafs, '{} for sample {}'.format(name, sample), plot_start, plot_finish, [start], [finish], [name], annotate)
        else:
          plot_bafs('{}.{}.png'.format(plot, name), bafs, name, plot_start, plot_finish, [start], [finish], [name], annotate)
        gene_starts.append(start)
        gene_finishes.append(finish)
        gene_names.append(name)
    if plot_chromosomes:
      logging.debug('plotting %s with regions: %s', last_chrom, ', '.join(gene_names))
      if sample is not None:
        plot_bafs('{}.{}.png'.format(plot, last_chrom), bafs, 'Chromosome {} for sample {}'.format(last_chrom, sample), None, None, gene_starts, gene_finishes, gene_names, False) # don't annotate chroms
      else:
        plot_bafs('{}.{}.png'.format(plot, last_chrom), bafs, 'Chromosome {}'.format(last_chrom), None, None, gene_starts, gene_finishes, gene_names, False) # don't annotate chroms
    if plot_custom is not None:
      chromosome, coords = plot_custom.split(':')
      if chromosome == last_chrom:
        logging.debug('plotting %s with regions: %s to %s', last_chrom, ', '.join(gene_names), plot_custom)
        start, finish = coords.split('-')
        if sample is not None:
          plot_bafs('{}.custom.png'.format(plot), bafs, 'Chromosome {} from {} to {:.0f}M for sample {}'.format(chromosome, start, int(finish) / 1000000, sample), int(start), int(finish), gene_starts, gene_finishes, gene_names, False) # don't annotate custom
        else:
          plot_bafs('{}.custom.png'.format(plot), bafs, 'Chromosome {} from {} to {:.0f}M'.format(chromosome, start, int(finish) / 1000000), int(start), int(finish), gene_starts, gene_finishes, gene_names, False) # don't annotate custom
  elif plot_chromosomes:
    if sample is not None:
      plot_bafs('{}.{}.png'.format(plot, last_chrom), bafs, 'Chromosome {} for sample {}'.format(last_chrom, sample), False) # don't annotate chroms
    else:
      plot_bafs('{}.{}.png'.format(plot, last_chrom), bafs, 'Chromosome {}'.format(last_chrom), False) # don't annotate chroms

def calculate_segments(min_len, min_prop, noheader, plot, regions, region_names, region_padding, plot_chromosomes, annotate, plot_custom, sample):
  start = None # last potential end
  last = None # last potential end
  chrom = None
  last_chrom = None
  loh_status = {'loh': False, 'accepts': 0, 'supports': 0, 'neutrals': 0}
  stats = collections.defaultdict(int)
  bafs = {'x': [], 'y': [], 'c': [], 'calls': [], 'l': []}

  logging.info('regions: %s', ', '.join(regions))
  logging.info('region names: %s', ', '.join(region_names))

  if not noheader:
    sys.stdout.write("Chr\tStart\tEnd\tAccept_Pct\tAccepts\tSupports\tNeutrals\tLength\n")

  logging.info('reading from stdin...')
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

      if annotate:
        bafs['l'].append(str(pos))

    if status.startswith('reject') or last_chrom != chrom: # end of loh
      write_merged(start, last, loh_status, stats, min_len, min_prop, bafs)
      loh_status = {'loh': False, 'accepts': 0, 'supports': 0, 'neutrals': 0}
      start = None
      last = None

      if plot is not None and last_chrom != chrom:
        check_plot(plot, last_chrom, chrom, regions, region_names, region_padding, bafs, plot_chromosomes, annotate, plot_custom, sample)
        bafs = {'x': [], 'y': [], 'c': [], 'calls': [], 'l': []}

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
  if plot is not None and chrom is not None:
    check_plot(plot, last_chrom, chrom, regions, region_names, region_padding, bafs, plot_chromosomes, annotate, plot_custom, sample)

  logging.info('done: %s', ', '.join(['{}: {}'.format(k, stats[k]) for k in stats]))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='LOH caller')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--min_len', default=1000, type=int, help='minimum segment size')
  parser.add_argument('--min_prop', type=float, default=0.1,  help='minimum accept proportion')
  parser.add_argument('--noheader', action='store_true',  help='no header for bed compatibility')
  parser.add_argument('--regions', nargs='+',  help='regions of interest to check for LOH and plot e.g. chr2:47628206-47712367')
  parser.add_argument('--region_names', nargs='+',  help='names of regions of interest')
  parser.add_argument('--region_padding', nargs='+', help='names of regions of interest')
  parser.add_argument('--plot_custom', required=False, help='custom region to plot like a chromosome')
  parser.add_argument('--sample', required=False, help='sample name to put in title')
  parser.add_argument('--plot', required=False, help='plot results with image prefix')
  parser.add_argument('--plot_chromosomes', action='store_true', help='plot results with image prefix')
  parser.add_argument('--annotate', action='store_true', help='include variant locations on region plots')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  if args.regions is not None and args.region_padding is None:
    region_padding = [10000] * len(args.regions)
  else:
    region_padding = args.region_padding

  if args.regions is not None and len(args.regions) != len(region_padding):
    logging.error('Region lengths do not match')
    sys.exit(1)
  calculate_segments(args.min_len, args.min_prop, args.noheader, args.plot, args.regions, args.region_names, region_padding, args.plot_chromosomes, args.annotate, args.plot_custom, args.sample)
