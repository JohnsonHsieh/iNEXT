\name{iNEXT-package}
\alias{iNEXT-package}
\docType{package}
\title{
Interpolation and extrapolation for species diversity
}
\description{
iNEXT (iNterpolation and EXTrapolation) provides functions to compute and plot two types (sample-size- and coverage-based) interpolation and extrapolation sampling curves of Hill numbers for three most widely used members of Hill numbers (species richness, Shannon diversity and Simpson diversity) for individual-based abundance data or sampling-unit-based incidence data. iNEXT also computes bootstrap confidence intervals around the diversity for rarefied/extrapolated samples, facilitating the comparisons of diversities across multiple assemblages/sites. The estimated asymptote along with a confidence interval for each of the three diversity measures is also provided. An auxiliary function is included to compute/compare diversities across multiple assemblages for a particular user-specified sample size or sample coverage. The sample-size-based rarefaction and extrapolation for species richness were developed by Colwell et al. (2012) and the corresponding coverage-based methodologies were developed by Chao and Jost (2012). Chao et al. (2014) extended the previous work for species richness to Hill numbers. The statistical methods and tools provided in iNEXT efficiently use all data to make more robust and detailed inferences about the sampled assemblages, and also to make objective comparisons of multiple assemblages. A short review of the theoretical background and a brief description of methods are included in an application paper by Hsieh, Ma & Chao (2016). An online version (https://chao.shinyapps.io/iNEXTOnline/) is also available for users without an R background. 
}
\author{
T. C. Hsieh \cr
K. H. Ma \cr
Anne Chao \cr
Maintainer: T. C. Hsieh <euler96@gmail.com>
}
\references{
Chao, A., Gotelli, N.J., Hsieh, T.C., Sander, E.L., Ma, K.H., Colwell, R.K. & Ellison, A.M. (2014) Rarefaction and extrapolation with Hill numbers: a framework for sampling and estimation in species diversity studies. Ecological Monographs, 84, 45-67. \cr\cr
Chao, A. & Jost, L. (2012) Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size. Ecology, 93, 2533-2547.\cr\cr
Colwell, R.K., Chao, A., Gotelli, N.J., Lin, S.-Y., Mao, C.X., Chazdon, R.L. & Longino, J.T. (2012) Models and estimators linking individual-based and sample-based rarefaction, extrapolation and comparison of assemblages. Journal of Plant Ecology, 5, 3-21. \cr\cr
Hsieh, T.C., Ma, K.H. & Chao, A. (2016) iNEXT: An R package for rarefaction and extrapolation of species diversity (Hill numbers). Methods in Ecology and Evolution, 7, 1451-1456. 
}
