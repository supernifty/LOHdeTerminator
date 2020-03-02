#!/usr/bin/env python
'''
  estimate % tumour from allele frequencies
'''

import argparse
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
  return {'tumour': np.percentile(np.array([v for v in values if v > 0.1]), 99)}

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

def read_vcf(fn, pass_only, dp_threshold, info_af):
  logging.info('reading vcf from stdin...')
  skipped_dp = skipped_pass = 0

  vcf_in = cyvcf2.VCF(fn)
  values = []

  for variant_count, variant in enumerate(vcf_in):
    # calculate vaf
    if len(variant.ALT) > 1:
      logging.warn('variant %i is multi-allelic', variant_count + 1)

    is_pass = variant.FILTER is None or variant.FILTER == 'alleleBias'
    if pass_only and not is_pass:
      skipped_pass += 1
      continue

    if variant.INFO["DP"] < dp_threshold: # somatic + germline
      skipped_dp += 1
      continue

    if info_af:
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

def main(pass_only, dp_threshold, info_af, plot, trials):
  values = read_vcf('-', pass_only, dp_threshold, info_af)
  result = estimate('percentile', values, pass_only, dp_threshold, info_af)
  sys.stdout.write('Estimated tumour percentage (percentile):\t{:.2f}\n'.format(result['tumour']))
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

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='estimate tumour %')
  parser.add_argument('--dp_threshold', required=False, default=50, help='use af in info field')
  parser.add_argument('--pass_only', action='store_true', help='only pass variants')
  parser.add_argument('--info_af', action='store_true', help='use af in info field')
  parser.add_argument('--trials', required=False, type=int, default=1, help='how many runs of model')
  parser.add_argument('--plot', required=False, help='use af in info field')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.pass_only, args.dp_threshold, args.info_af, args.plot, args.trials)
