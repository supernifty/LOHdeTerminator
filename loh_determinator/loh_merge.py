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

#import version
PROGRAM_VERSION="0.1"

#LOH_COLOR='#00852b'
#GL_COLOR='#303030'
#SUPPORT_COLOR='#1e5aa8'
#NOT_COLOR='#ca0020'
#NEUTRAL_COLOR='#bababa'

NOT_COLOR='#00852b' # greenish
GL_COLOR='#303030' # grey
SUPPORT_COLOR='#590696' # maroon
LOH_COLOR='#ca0020' # red
NEUTRAL_COLOR='#bababa' #light grey

ALPHA = {
  LOH_COLOR: 0.7,
  SUPPORT_COLOR: 0.1
}

CONNECT=(LOH_COLOR, SUPPORT_COLOR)

def write_merged(start, last, loh_status, stats, min_len, min_prop, min_count, bafs, max_loh_without_evidence, min_accept):
  if loh_status['loh'] and loh_status['accepts'] >= min_accept:
    if 'last_accept_or_support' not in loh_status:
      new_last = last[1]
    else:
      new_last = min(last[1], loh_status['last_accept_or_support'][1] + max_loh_without_evidence)
    length = new_last - start[1]
    prop = loh_status['accepts'] / (loh_status['accepts'] + loh_status['supports'] + loh_status['neutrals'])
    if length > min_len and prop > min_prop and loh_status['accepts'] + loh_status['supports'] >= min_count:
      sys.stdout.write("{chrom}\t{start}\t{end}\t{accept_pct:.1f}\t{accepts}\t{supports}\t{neutrals}\t{length}\n".format(chrom=start[0], start=start[1], end=new_last, accept_pct=100. * prop, accepts=loh_status['accepts'], supports=loh_status['supports'], neutrals=loh_status['neutrals'], length=length))
      stats['count'] += 1
      stats['length'] += length
      stats['count-{}'.format(start[0])] += 1
      stats['length-{}'.format(start[0])] += length
      bafs['calls'].append((start[1], new_last))
      return True
  return False

def plot_bafs(filename, bafs, sample, start, finish, gene_starts, gene_finishes, gene_names, annotate, width, height, no_legend, annotation_style, title, legend_location, markersize, dpi, tumour_cellularity):
  somatic_patch = germline_patch = None
  if len(bafs['x']) > 0:
    x_mb = [x / 1000000 for x in bafs['x']]
    plt.figure(figsize=(width, height))
    if start is None:
      somatic_patch = plt.scatter(x_mb, bafs['y'], c=bafs['c'], alpha=0.5, label='Somatic')
      start = 0
      finish = max(bafs['x'])
      if annotate:
        for x, y, l in zip(x_mb, bafs['y'], bafs['l']):
          plt.annotate(l, x, y)
      # germline
      germline_patch = plt.scatter(x_mb, bafs['g'], c=GL_COLOR, alpha=0.7, marker='+', label='Germline')

      for xyc in zip(x_mb, bafs['y'], bafs['c'], bafs['g']):
        if xyc[2] in CONNECT:
          plt.vlines(xyc[0], xyc[1], xyc[3], colors=xyc[2], alpha=ALPHA[xyc[2]])

    else:
      xycs = [xyc for xyc in zip(x_mb, bafs['y'], bafs['c'], bafs['g']) if start / 1000000 <= xyc[0] <= finish / 1000000]
      if len(xycs) > 0:
        somatic_patch = plt.scatter([xyc[0] for xyc in xycs], [xyc[1] for xyc in xycs], c=[xyc[2] for xyc in xycs], alpha=0.5, label='Somatic')
        if annotate:
          for x, y, l in zip(x_mb, bafs['y'], bafs['l']):
            plt.annotate(s=l, xy=(x, y))
        # germline
        germline_patch = plt.scatter([xyc[0] for xyc in xycs], [xyc[3] for xyc in xycs], c=GL_COLOR, alpha=0.7, marker='+', label='Germline')

        # connect germline/somatic calls of interest
        for xyc in xycs:
          if xyc[2] in CONNECT:
            plt.vlines(xyc[0], xyc[1], xyc[3], colors=xyc[2], alpha=ALPHA[xyc[2]])

    if title is None:
      plt.title('Genomic regions containing evidence for LOH across {}'.format(sample))
    else:
      plt.title(title)
    #plt.title('{}'.format(title))
    if tumour_cellularity is not None:
      plt.ylabel('Adjusted Allele Fraction')
    else:
      plt.ylabel('Allele Fraction')
    plt.xlabel('Genomic Position (MB)')
    plt.gca().set_ylim([0.0, 1.0])
    plt.gca().set_xlim([start / 1000000, finish / 1000000])

    plt.grid(which='major', axis='y', linewidth=1)

    # add calls
    for call in bafs['calls']:
      if start is None:
        plt.axvspan(call[0] / 1000000, call[1] / 1000000, facecolor=LOH_COLOR, alpha=0.2)
      else:
        if call[0] <= finish and call[1] >= start:
          call_s = max(start, call[0])
          call_f = min(finish, call[1])
          plt.axvspan(call_s / 1000000, call_f / 1000000, facecolor=LOH_COLOR, alpha=0.2)

    # highlight gene
    if gene_starts is not None:
      for gene_start, gene_finish, gene_name in zip(gene_starts, gene_finishes, gene_names):
        logging.debug('annotating %s at %i %i in %i %i', gene_name, gene_start, gene_finish, start, finish)
        if (gene_finish - gene_start) < (finish - start) * 0.005:
          gene_finish = gene_start + (finish - start) * 0.005
        logging.debug('annotating %s at %i %i', gene_name, gene_start, gene_finish)
        if annotation_style == 'flag':
          # like a flag
          plt.axvspan(gene_start / 1000000, gene_finish / 1000000, ymin=0.0, ymax=0.15, facecolor='#000000', alpha=0.8)
          plt.text(gene_start / 1000000, 0.16, gene_name, horizontalalignment='left', verticalalignment='bottom', bbox={'facecolor': 'yellow', 'alpha': 0.5})
        elif annotation_style == 'float':
          # not like a flag
          plt.axvspan(gene_start / 1000000, gene_finish / 1000000, ymin=0.13, ymax=0.16, facecolor='#000000', alpha=0.8)
          plt.text(gene_start / 1000000, 0.18, gene_name, horizontalalignment='left', verticalalignment='bottom', bbox={'facecolor': 'yellow', 'alpha': 0.5})
        else:
          logging.warn('unsupported annotation style %s', annotation_style)

    import matplotlib.patches as mpatches
    l1 = mpatches.Patch(color=LOH_COLOR, label='LOH')
    l2 = mpatches.Patch(color=NOT_COLOR, label='Not LOH')
    l3 = mpatches.Patch(color=SUPPORT_COLOR, label='Favours LOH')
    l4 = mpatches.Patch(color=NEUTRAL_COLOR, label='Neutral')
    
    if not no_legend:
      if somatic_patch is not None and germline_patch is not None:
        plt.legend([l1, l2, l3, l4, somatic_patch, germline_patch],['LOH', 'Not LOH', 'Favours LOH', 'Neutral', 'Somatic', 'Germline'], loc=legend_location, ncol=7, labelspacing=0., framealpha=0.7)
      else:
        plt.legend([l1, l2, l3, l4],['LOH', 'Not LOH'], loc=legend_location, ncol=5, labelspacing=0., framealpha=0.7)

    plt.tight_layout()
    plt.savefig(filename, dpi=dpi)
    plt.close()

    bafs = {'x': [], 'y': [], 'c': []}
    logging.info('wrote %s', filename)

def check_plot(plot, last_chrom, chrom, regions, region_names, region_padding, bafs, plot_chromosomes, annotate, plot_custom, sample, width, height, no_legend, annotation_style, title, legend_location, markersize, dpi, plot_type, tumour_cellularity):
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
          #plot_bafs('{}.{}.png'.format(plot, name), bafs, '{} for sample {}'.format(name, sample), plot_start, plot_finish, [start], [finish], [name], annotate, width, height)
          plot_bafs('{}.{}.{}'.format(plot, name, plot_type), bafs, 'Sample {}'.format(sample), plot_start, plot_finish, [start], [finish], [name], annotate, width, height, no_legend, annotation_style, title, legend_location, markersize, dpi, tumour_cellularity)
        else:
          plot_bafs('{}.{}.{}'.format(plot, name, plot_type, ), bafs, name, plot_start, plot_finish, [start], [finish], [name], annotate, width, height, no_legend, annotation_style, title, legend_location, markersize, dpi, tumour_cellularity)
        gene_starts.append(start)
        gene_finishes.append(finish)
        gene_names.append(name)
    if plot_chromosomes:
      logging.debug('plotting %s with regions: %s', last_chrom, ', '.join(gene_names))
      if sample is not None:
        plot_bafs('{}.{}.{}'.format(plot, last_chrom, plot_type), bafs, 'Sample {}'.format(sample), None, None, gene_starts, gene_finishes, gene_names, False, width, height, no_legend, annotation_style, title, legend_location, markersize, dpi, tumour_cellularity)
      else:
        plot_bafs('{}.{}.{}'.format(plot, last_chrom, plot_type), bafs, 'Chromosome {}'.format(last_chrom), None, None, gene_starts, gene_finishes, gene_names, False, width, height, no_legend, annotation_style, title, legend_location, markersize, dpi, tumour_cellularity)
    if plot_custom is not None:
      chromosome, coords = plot_custom.split(':')
      if chromosome == last_chrom:
        logging.debug('plotting %s with regions: %s to %s', last_chrom, ', '.join(gene_names), plot_custom)
        start, finish = coords.split('-')
        if sample is not None:
          plot_bafs('{}.custom.{}'.format(plot, plot_type), bafs, 'Sample {}'.format(sample), int(start), int(finish), gene_starts, gene_finishes, gene_names, False, width, height, no_legend, annotation_style, title, legend_location, markersize, dpi, tumour_cellularity)
        else:
          plot_bafs('{}.custom.{}'.format(plot, plot_type), bafs, 'Chromosome {} from {:.0f}M to {:.0f}M'.format(chromosome, int(start) / 1000000, int(finish) / 1000000), int(start), int(finish), gene_starts, gene_finishes, gene_names, False, width, height, no_legend, annotation_style, title, legend_location, markersize, dpi, tumour_cellularity)
  elif plot_chromosomes:
    if sample is not None:
      plot_bafs('{}.{}.{}'.format(plot, last_chrom, plot_type), bafs, 'Chromosome {} for sample {}'.format(last_chrom, sample), None, None, [], [], [], False, width, height, no_legend, annotation_style, title, legend_location, markersize, dpi, tumour_cellularity)
    else:
      plot_bafs('{}.{}.{}'.format(plot, last_chrom, plot_type), bafs, 'Chromosome {}'.format(last_chrom), None, None, [], [], [], False, width, height, no_legend, annotation_style, title, legend_location, markersize, dpi, tumour_cellularity)

def calculate_segments(min_len, min_prop, min_count, noheader, plot, regions, region_names, region_padding, plot_chromosomes, annotate, plot_custom, sample, width, height, no_legend, annotation_style, title, legend_location, max_loh_without_evidence, fontsize, markersize, dpi, plot_type, tumour_cellularity, min_accept):
  matplotlib.rcParams.update({'font.size': fontsize})
  start = None # last potential end
  last = None # last potential end
  chrom = None
  last_chrom = None
  loh_status = {'loh': False, 'accepts': 0, 'supports': 0, 'neutrals': 0}
  stats = collections.defaultdict(int)
  bafs = {'x': [], 'y': [], 'g': [], 'c': [], 'calls': [], 'l': []}

  if regions is not None:
    logging.info('regions: %s', ', '.join(regions))

  if region_names is not None:
    logging.info('region names: %s', ', '.join(region_names))

  if not noheader:
    sys.stdout.write("Chr\tStart\tEnd\tAccept_Pct\tAccepts\tSupports\tNeutrals\tLength\n")

  if tumour_cellularity is not None and tumour_cellularity > 1.0:
    tumour_cellularity /= 100
    logging.info('adjusting bafs by %.2f', tumour_cellularity)

  logging.info('reading from stdin...')
  for line in sys.stdin: # each point of evidence
    if line.startswith('chrom'): # header
      continue
    fields = line.strip('\n').split('\t')
    chrom = fields[0]
    pos = int(fields[1])
    gaf = float(fields[2]) # germline af
    baf = float(fields[3]) # tumour af

    if tumour_cellularity is not None:
      baf = min(baf / tumour_cellularity, 1.0)

    status = fields[6]

    if plot is not None:
      bafs['x'].append(pos)
      bafs['y'].append(baf)
      bafs['g'].append(gaf)
      if status.startswith('accept'):
        bafs['c'].append(LOH_COLOR)
      elif status.startswith('reject'):
        bafs['c'].append(NOT_COLOR)
      elif status.startswith('support'):
        bafs['c'].append(SUPPORT_COLOR)
      else: # neutral
        bafs['c'].append(NEUTRAL_COLOR)

      if annotate:
        bafs['l'].append(str(pos))

    if status.startswith('reject') or last_chrom != chrom or pos - last_pos > max_loh_without_evidence: # end of loh
      write_merged(start, last, loh_status, stats, min_len, min_prop, min_count, bafs, max_loh_without_evidence, min_accept)
      loh_status = {'loh': False, 'accepts': 0, 'supports': 0, 'neutrals': 0}
      start = None
      last = None

      if plot is not None and last_chrom != chrom:
        check_plot(plot, last_chrom, chrom, regions, region_names, region_padding, bafs, plot_chromosomes, annotate, plot_custom, sample, width, height, no_legend, annotation_style, title, legend_location, markersize, dpi, plot_type, tumour_cellularity)
        bafs = {'x': [], 'y': [], 'g': [], 'c': [], 'calls': [], 'l': []}

    if not status.startswith('reject'): # potential loh region
      if start is None:
        start = (chrom, pos)
      loh_status['loh'] = True
      if status.startswith('accept'):
        loh_status['accepts'] += 1
        loh_status['last_accept_or_support'] = (chrom, pos)
      elif status.startswith('support'):
        loh_status['supports'] += 1
        loh_status['last_accept_or_support'] = (chrom, pos)
      elif status.startswith('neutral'):
        loh_status['neutrals'] += 1

      last = (chrom, pos) # last seen of any non-reject

    last_chrom = chrom
    last_pos = pos

  write_merged(start, last, loh_status, stats, min_len, min_prop, min_count, bafs, max_loh_without_evidence, min_accept)
  if plot is not None and chrom is not None:
    check_plot(plot, last_chrom, chrom, regions, region_names, region_padding, bafs, plot_chromosomes, annotate, plot_custom, sample, width, height, no_legend, annotation_style, title, legend_location, markersize, dpi, plot_type, tumour_cellularity)

  logging.info('done: %s', ', '.join(['{}: {}'.format(k, stats[k]) for k in stats]))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='LOH caller')
  parser.add_argument('--version', action='version', version=PROGRAM_VERSION)
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--min_len', default=1000, type=int, help='minimum segment size')
  parser.add_argument('--min_prop', type=float, default=0.1,  help='minimum accept proportion')
  parser.add_argument('--min_count', type=int, default=1,  help='minimum number of non-negative variants required')
  parser.add_argument('--min_accept', type=int, default=0,  help='minimum number of positive variants required')
  parser.add_argument('--max_loh_without_evidence', type=int, default=1000000, help='maximum size of loh without evidence')
  parser.add_argument('--nolegend', action='store_true',  help='no legend on plot')
  parser.add_argument('--legend_location', default='lower left', help='location')
  parser.add_argument('--noheader', action='store_true',  help='no header for bed compatibility')
  parser.add_argument('--regions', nargs='+',  help='regions of interest to check for LOH and plot e.g. chr2:47628206-47712367')
  parser.add_argument('--region_names', nargs='+',  help='names of regions of interest')
  parser.add_argument('--region_padding', nargs='+', help='names of regions of interest')
  parser.add_argument('--plot_custom', required=False, help='custom region to plot like a chromosome')
  parser.add_argument('--sample', required=False, help='sample name to put in title')
  parser.add_argument('--title', required=False, help='title for plot')
  parser.add_argument('--plot', required=False, help='plot results with image prefix')
  parser.add_argument('--plot_type', required=False, default='png', help='plot format png pdf svg')
  parser.add_argument('--plot_chromosomes', action='store_true', help='plot results with image prefix')
  parser.add_argument('--width', default=24, type=int, help='width of image')
  parser.add_argument('--height', default=8, type=int, help='height of image')
  parser.add_argument('--fontsize', default=8, type=int, help='font size')
  parser.add_argument('--markersize', default=12, type=int, help='font size')
  parser.add_argument('--annotate', action='store_true', help='include variant locations on region plots')
  parser.add_argument('--tumour_cellularity', required=False, type=float, help='if specified, tumour AF will be adjusted')
  parser.add_argument('--dpi', default=300, type=int, help='dpi')
  parser.add_argument('--annotation_style', choices=['float', 'flag'], default='float', required=False, help='how to annotate regions of interest')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  if args.regions is not None and args.region_padding is None:
    region_padding = [10000] * len(args.regions)
  elif args.regions is not None and len(args.region_padding) == 1 and len(args.regions) != 1:
    region_padding = args.region_padding * len(args.regions)
  else:
    region_padding = args.region_padding

  if args.regions is not None and len(args.regions) != len(region_padding):
    logging.error('Region lengths do not match')
    sys.exit(1)
  calculate_segments(args.min_len, args.min_prop, args.min_count, args.noheader, args.plot, args.regions, args.region_names, region_padding, args.plot_chromosomes, args.annotate, args.plot_custom, args.sample, args.width, args.height, args.nolegend, args.annotation_style, args.title, args.legend_location, args.max_loh_without_evidence, args.fontsize, args.markersize, args.dpi, args.plot_type, args.tumour_cellularity, args.min_accept)
