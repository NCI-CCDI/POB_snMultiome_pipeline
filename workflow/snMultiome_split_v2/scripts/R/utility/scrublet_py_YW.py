
import matplotlib
matplotlib.use('agg')
import numpy as np
import scrublet as scr
import scipy.io
import os
from scipy.sparse import csc_matrix
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from tabulate import tabulate
import h5py
import pandas as pd
from scipy.sparse import csc_matrix
import math
from sklearn.mixture import GaussianMixture
from scipy.optimize import root_scalar
from scipy.stats import norm
from collections import Counter
import random

random.seed(42) 

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica', 'sans-serif']
plt.rc('font', size=8)

#np.random.seed(42) 
# Define the PDFs of the two normal distributions
def pdf1(x, mu1, sigma1):
    return np.exp(-0.5 * ((x - mu1) / sigma1) ** 2) / (np.sqrt(2 * np.pi) * sigma1)

def pdf2(x, mu2, sigma2):
    return np.exp(-0.5 * ((x - mu2) / sigma2) ** 2) / (np.sqrt(2 * np.pi) * sigma2)

# Define the function to find the difference between the two PDFs
def difference(x, mu1, sigma1, mu2, sigma2):
    return pdf1(x, mu1, sigma1) - pdf2(x, mu2, sigma2)
  
def fit_Gaussian(data, mu1, sigma1, mu2, sigma2,doublet_rate,doublet_rate_adj,doublet_scores,cutoff_auto,sample):
    # Generate x values
    #x = np.linspace(min(mu1 - 3 * sigma1, mu2 - 3 * sigma2), max(mu1 + 3 * sigma1, mu2 + 3 * sigma2), 1000)
    min_val = min(mu1 - 3 * sigma1, mu2 - 3 * sigma2)
    if min_val <0:
        min_val = 0
    max_val = max(mu1 + 3 * sigma1, mu2 + 3 * sigma2)
    if max_val >1:
        min_val = 1
    x = np.linspace(min_val, max_val, 1000) 
    # # Compute PDFs for the two normal distributions
    pdf_v1 = norm.pdf(x, mu1, sigma1)
    pdf_v2 = norm.pdf(x, mu2, sigma2)

    #print(pdf_v1)
    #print(pdf_v2)
    # Find the intersection point numerically
    max_attempts = 100
    attempt = 1
    while attempt <= max_attempts:
        try:
            intersection = root_scalar(difference, args=(mu1, sigma1, mu2, sigma2), bracket=[0, 1])
            # If the operation succeeds, break out of the loop
            break
        except ValueError:
            attempt += 1    # Increment the attempt counter
    else:
        # This block executes if the loop completes without a successful attempt
        raise ValueError("Maximum attempts reached. Unable to find intersection.")
    intersection_val = round(intersection.root,3)
    cutoff_auto_val =  round(cutoff_auto,3)
    # Print the intersection point
    print("Intersection point:", intersection.root)
    print ("Find intersection attempt:", attempt)
    # Define the percentile value
    percentile =  round(100 - doublet_rate * 100, 2)  # Example: 75th percentile
    percentile_adj = round(100 - doublet_rate_adj * 100, 2)
    # Compute the percentile score
    percentile_score = round(np.percentile(doublet_scores, percentile),3)
    print(f"{percentile}th percentile score:", percentile_score)
    percentile_score_adj = round(np.percentile(doublet_scores, percentile_adj),3)
    print(f"Adj {percentile_adj}th percentile score:", percentile_score_adj)
    
    filter1 = doublet_scores[doublet_scores > intersection.root]
    print("filtered out by gmm cutoff:",len(filter1))
    filter2 = doublet_scores[doublet_scores > percentile_score]
    print("filtered out by percentile cutoff:",len(filter2))
    filter3 = doublet_scores[doublet_scores > cutoff_auto]
    print("filtered out by auto cutoff:",len(filter3))
    filter4 = doublet_scores[doublet_scores > percentile_score_adj]
    print("filtered out by ajust percentile cutoff:",len(filter4))
    
    # Plot the two normal distributions
    plt.figure(figsize=(8, 6)) 
    plt.hist(data, np.linspace(0, 1, 50), linewidth=0, alpha=0.8,density=True,color='skyblue',  label = 'Histogram')
    #plt.hist(data, bins=50, density=True, alpha=0.5, label='Histogram')
    plt.plot(x, pdf_v1,color='orange',  label='First Normal Distribution\n(μ={}, σ={})'.format(mu1, sigma1))
    plt.plot(x, pdf_v2, color='green', label='Second Normal Distribution\n(μ={}, σ={})'.format(mu2, sigma2))
    plt.xlabel('x')
    plt.ylabel('Probability Density')
    plt.title(f'{sample}\nTwo Normal Distributions')

    plt.grid(False)
    plt.axvline(intersection.root, color="red",label=f'intersection_{intersection_val}\n(Doublet:{len(filter1)})')
    plt.axvline(percentile_score_adj, color="red",linestyle='--', label=f'{percentile_adj}th Percentile_{percentile_score_adj}\n(Doublet:{len(filter4)})')
    plt.axvline(percentile_score, color="blue", label=f'{percentile}th Percentile_{percentile_score}\n(Doublet:{len(filter2)})')
    plt.axvline(cutoff_auto, color="blue",linestyle='--', label=f'scrublet auto threshold_{cutoff_auto_val}\n(Doublet:{len(filter3)})')
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=2)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=2, fontsize=8)
    plt.subplots_adjust(bottom=0.3)
    #plt.legend(loc='upper right')
    #plt.show()
    intersection_val =  round(intersection.root,3)
    percentile_cutoff = percentile_score
    percentile_cutoff_adj = percentile_score_adj
    print("round to 3 decimal")
    print(intersection.root, intersection_val, percentile_score, percentile_cutoff, percentile_score_adj, percentile_cutoff_adj)
    return intersection_val, percentile_cutoff,percentile, percentile_cutoff_adj, percentile_adj, attempt, plt

def combine_plot(fig1,fig2,sample,criteria,path):
    file_name = f'{path}/{sample}_scrublet_hist_umap_{criteria}.pdf'
    print(file_name)
    with PdfPages(file_name) as pdf:
        pdf.savefig(fig1)
        pdf.savefig(fig2)
    plt.close()

def scrublet_py_YW(i, j, val, dim, expected_doublet_rate,doublet_rate_adj, min_counts, min_cells, min_gene_variability_pctl, n_prin_comps, sim_doublet_ratio, n_neighbors,sample,path):
    data = csc_matrix((val, (i, j)), shape = dim)

    doublet_rate = expected_doublet_rate
    doublet_rate_adj = doublet_rate_adj
    scrub = scr.Scrublet(data, expected_doublet_rate=doublet_rate, sim_doublet_ratio=int(sim_doublet_ratio), n_neighbors=int(n_neighbors))
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=int(min_counts), min_cells=int(min_cells), min_gene_variability_pctl=min_gene_variability_pctl, n_prin_comps=int(n_prin_comps))
    cutoff_auto = round(scrub.threshold_,3)
    n_cells=data.shape[0]
    print(n_cells)
    print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    print('Done.')
    score_data = scrub.doublet_scores_sim_
    #print(score_data)
    # Fit a Gaussian Mixture Model with 2 components
    #max_attempts = 100
    attempt = 1

    gmm = GaussianMixture(n_components=2)
    gmm.fit(score_data.reshape(-1, 1))
    means = gmm.means_.flatten()
    variances = gmm.covariances_.flatten()
    # Known means and standard deviations of the two normal distributions
    mean1 = round(means[0], 3)
    sigma1 = round(np.sqrt(variances[0])*1.5 , 3)
    mean2 = round(means[1], 3)
    sigma2 = round(np.sqrt(variances[1]), 3)
    if (mean1 > mean2):
        mean1 = round(means[1], 3)
        sigma1 = round(np.sqrt(variances[1])*1.5 , 3)
        mean2 = round(means[0], 3)
        sigma2 = round(np.sqrt(variances[0]), 3)
        
#    while True:
#        try:
#            gmm = GaussianMixture(n_components=2)
#            gmm.fit(score_data.reshape(-1, 1))
#            # Get means and variances of the two components
#            means = gmm.means_.flatten()
#            variances = gmm.covariances_.flatten()
#            # Known means and standard deviations of the two normal distributions
#            mean1 = round(means[0], 3)
#            sigma1 = round(np.sqrt(variances[0])*1.5 , 3)
#            mean2 = round(means[1], 3)
#            sigma2 = round(np.sqrt(variances[1]), 3)
#            print(mean1, sigma1, mean2, sigma2)
#        except Exception as e:
#            # Handle any exceptions raised during fitting
#            print(f"Error occurred: {e}")
#            # Continue to the next iteration of the loop to retry
#        # Check if the condition is met (mean2 < mean1)
#        if mean2 > mean1:
#            # If the condition is met, break out of the loop
#            break
#        else:
#            if attempt == max_attempts:
#                print("Attempt" + str(attempt) + "; Condition not met. STOPPED!!!")
#                sys.exit()
#            else:
#                # If the condition is not met, print a message and continue trying
#                print("Attempt" + str(attempt) + "; Condition not met. Retrying...")
#                attempt += 1
    print(attempt)
    print(mean1,sigma1, mean2, sigma2)
    cutoff_gmm1, cutoff_percentile1,percentile,cutoff_percentile1_adj,percentile_adj,attempt2,plot1 = fit_Gaussian (score_data,mean1,sigma1,mean2,sigma2,doublet_rate,doublet_rate_adj,doublet_scores,cutoff_auto,sample)
    filename = f'{path}/{sample}_scrublet_1st_round_fit.pdf'
    #print(filename)
    plot1.savefig(filename, format='pdf')
    plot1.close()
    criteria = "gmm"
    scrub.call_doublets(threshold=cutoff_gmm1)
    # Generate the histogram plot
    fig1_tuple = scrub.plot_histogram()
    fig1 = fig1_tuple[0]
    # Generate the embedding plot
    fig2_tuple = scrub.plot_embedding('UMAP', order_points=True)
    fig2 = fig2_tuple[0]
    combine_plot(fig1,fig2,sample,criteria,path)
    
    criteria = "percentile"
    scrub.call_doublets(threshold=cutoff_percentile1)
    fig3_tuple = scrub.plot_histogram()
    fig3 = fig3_tuple[0]
    fig4_tuple = scrub.plot_embedding('UMAP', order_points=True)
    fig4 = fig4_tuple[0]
    combine_plot(fig3,fig4,sample,criteria,path)

    criteria = "percentile_adj"
    scrub.call_doublets(threshold=cutoff_percentile1_adj)
    fig7_tuple = scrub.plot_histogram()
    fig7 = fig7_tuple[0]
    fig8_tuple = scrub.plot_embedding('UMAP', order_points=True)
    fig8 = fig8_tuple[0]
    combine_plot(fig7,fig8,sample,criteria,path)
    
    criteria = "auto"
    scrub.call_doublets(threshold=cutoff_auto)
    fig5_tuple = scrub.plot_histogram()
    fig5 = fig5_tuple[0]
    fig6_tuple = scrub.plot_embedding('UMAP', order_points=True)
    fig6 = fig6_tuple[0]
    combine_plot(fig5,fig6,sample,criteria,path)
    cutoff_final = cutoff_gmm1
    scr_cri = "gmm"
#    if (attempt > 1 or attempt2 > 1):
    if (attempt2 > 1) :
      print("attempt fit:", attempt, "\n", "attempt intersection:",attempt2,"\n")
      cutoff_final = cutoff_percentile1_adj
      scr_cri = "percentile_adj"
    print ("auto,gmm,percentile,percentile_adj")
    print (cutoff_auto,cutoff_gmm1, cutoff_percentile1,cutoff_percentile1_adj)
    print (attempt,cutoff_final,scr_cri)
    return(doublet_scores, predicted_doublets,cutoff_auto,cutoff_gmm1, cutoff_percentile1,cutoff_final,scr_cri,percentile,cutoff_percentile1_adj,percentile_adj)
    
