import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.mixture import BayesianGaussianMixture
from scipy.stats import norm
import pandas as pd

def find_weighted_intersection(mu1, sigma1, w1, mu2, sigma2, w2):
    # Coefficients for the quadratic equation a*x^2 + b*x + c = 0
    a = 1 / (2 * sigma1**2) - 1 / (2 * sigma2**2)
    b = mu2 / (sigma2**2) - mu1 / (sigma1**2)
    c = mu1**2 / (2 * sigma1**2) - mu2**2 / (2 * sigma2**2) - np.log((sigma2 / sigma1) * (w1 / w2))
    roots = np.roots([a, b, c])
    # Return the root that lies between the two means
    return roots[np.logical_and(roots > mu1, roots < mu2)]
  
  # Step 5: Classify data into three levels based on thresholds
def classify_expression(value, gauss_n):
    if gauss_n == 4:
      if value <= threshold1:
          return 0  # 0 expression
      elif value <= threshold2:
          return 1  # Low expression
      elif value <= threshold3:
          return 2  # Medium expression
      else:
          return 3  # High expression
    elif gauss_n == 3:
      if value <= threshold1:
          return 0  # 0 expression
      elif value <= threshold2:
          return 1  # Low expression
      else:
          return 2  # Medium expression

# Generate or load gene expression data
gene_expression = np.genfromtxt('/Users/vickymoya/Desktop/Local code/HCR Annotation/Classification with Gaussian mixture/Npnt.csv', delimiter=',')  
gene_expression_nz = gene_expression[gene_expression > -5]

n_components_range = range(1, 7)

# Step 2: Fit GMMs and compute BIC and AIC for each
bic_scores = []
aic_scores = []
fit_scores = []

for n_components in n_components_range:
    gmm = GaussianMixture(n_components=n_components, random_state=42, max_iter = 500, init_params = 'k-means++', n_init = 5)
    gmm.fit(gene_expression_nz.reshape(-1, 1))
    bic_scores.append(gmm.bic(gene_expression_nz.reshape(-1, 1)))
    aic_scores.append(gmm.aic(gene_expression_nz.reshape(-1, 1)))
    fit_scores.append(gmm.score(gene_expression_nz.reshape(-1, 1)))

percent_bic_gain = np.abs(np.array(bic_scores[1:6])-np.array(bic_scores[0:5]))/np.array(bic_scores[0:5])
percent_aic_gain = np.abs(np.array(aic_scores[1:6])-np.array(aic_scores[0:5]))/np.array(aic_scores[0:5])
print(fit_scores)

# Step 3: Plot BIC and AIC scores
plt.figure(figsize=(10, 6))
plt.plot(n_components_range[1:], percent_bic_gain, label='BIC', marker='o')
plt.plot(n_components_range[1:], percent_aic_gain, label='AIC', marker='o')
plt.hlines(0.01, 2, 6)
plt.xlabel('Number of Components')
plt.ylabel('Score')
plt.title('BIC and AIC for Gaussian Mixture Models')
plt.legend()
plt.show()


# Step 1: Fit a two-component Gaussian Mixture Model
gmm = GaussianMixture(n_components=2, random_state=42, max_iter = 500, init_params = 'k-means++', n_init = 5)
#gmm = BayesianGaussianMixture(n_components=3, random_state=0, max_iter=500)
gmm.fit(gene_expression_nz.reshape(-1, 1))

# Step 2: Retrieve means and standard deviations
means = gmm.means_.flatten()
std_devs = np.sqrt(gmm.covariances_).flatten()
weights = gmm.weights_

# Ensure the lower mean is first
sorted_indices = np.argsort(means)
means = means[sorted_indices]
std_devs = std_devs[sorted_indices]
weights = weights[sorted_indices]

# Step 3: Find the intersection point
threshold = find_weighted_intersection(means[0], std_devs[0], weights[0], means[1], std_devs[1], weights[1])[0]

# Step 4: Binarize the data
binarized_expression = (gene_expression > threshold).astype(int)

# Output threshold and binary data
print("Threshold for binarization:", threshold)
print("Means:", means)
print("Stdevs:", std_devs)
print("Weight:", gmm.weights_)
#print("Binarized gene expression:", binarized_expression)

DF = pd.DataFrame(binarized_expression)
DF.to_csv("~/Desktop/Local code/HCR Annotation/Classification with Gaussian mixture/Ccdc80_bin.csv", header = False)

###### Trinarization ######

# Step 1: Fit a three-component Gaussian Mixture Model
gmm = GaussianMixture(n_components=3, random_state=42, max_iter = 500, init_params = 'k-means++', n_init = 5, means_init = np.array([[0.05],[2.2],[3.2]])) #
gmm.fit(gene_expression_nz.reshape(-1, 1))

# Step 2: Retrieve and sort means and standard deviations
means = gmm.means_.flatten()
std_devs = np.sqrt(gmm.covariances_).flatten()
weights = gmm.weights_

# Sort means and std_devs for consistent ordering
sorted_indices = np.argsort(means)
means = means[sorted_indices]
std_devs = std_devs[sorted_indices]
weights = weights[sorted_indices]

# Step 4: Find intersections between adjacent Gaussian components
threshold1 = find_weighted_intersection(means[0], std_devs[0], weights[0], means[1], std_devs[1], weights[1])[0]
threshold2 = find_weighted_intersection(means[1], std_devs[1], weights[1], means[2], std_devs[2], weights[2])[0]
threshold3 = find_weighted_intersection(means[2], std_devs[2], weights[2], means[3], std_devs[3], weights[3])[0]
threshold4 = find_weighted_intersection(means[3], std_devs[3], weights[3], means[4], std_devs[4], weights[4])[0]
gauss_n = 3

trinarized_expression = np.array([classify_expression(x, gauss_n) for x in gene_expression])

# Output thresholds and trinarized data
print("Thresholds for trinarization:", threshold1, threshold2)#, threshold3, threshold4)
print("Means:", means)
print("Stdevs:", std_devs)
print("Weight:", gmm.weights_)
#print("Trinarized gene expression levels:", trinarized_expression)

DF = pd.DataFrame(trinarized_expression)
DF.to_csv("~/Desktop/Local code/HCR Annotation/Classification with Gaussian mixture/Npnt_bin.csv", header = False)
