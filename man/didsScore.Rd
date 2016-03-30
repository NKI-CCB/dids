\name{didsScore}
\alias{didsScore}
\title{
Detection of Imbalanced Differential Signal (DIDS) - Scoring
}
\description{
Performs a DIDS analysis on an ExpressionSet to identify genes that are differentially expressed in (a subset of) the samples the set.
}
\usage{
didsScore(eset, modF = 2, modType = "tanh", alt = "two.sided") 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
ExpressionSet object on which the analysis should be performed. The phenoData of the object should contain a 'condition' column that marks each sample as being either 'sensitive' or 'resistant'.
}
  \item{modF}{
    Scoring function parameter. Corresponds to the exponent if the quadratic scoring function is used (x^modF) or a scaling factor if the tangent method is used (tanh((x-1)*(1+modF))+1).
}
  \item{modType}{
  Type of scoring function used. Options are 'tanh' and 'quad'. The tangent-based method (tanh) weighs small excess expression values (< 1) slightly lower than a linear function, differences between one and two slightly more than the linear case and all differences larger than two are capped at two. It is our experience that this scoring function detects most aberrant patterns. The quadratic scoring function (quad) penalizes small differences (< 1), but exaggerates differences larger than one. This scoring function is suited to detect patterns of infrequent aberrant samples that show a large difference compared to the control group.
}
  \item{alt}{
Specifies whether a one-sided or two-sided test should be performed (i.e. 'greater', 'less' or 'two.sided'). A two-sided test will perform 'less' and 'greater' separately and return the highest score between the two resulting scores for each feature/gene.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a dataframe with a score and p-value for each gene that reflect how likely the corresponding gene is actually differentially expressed. 
}
\references{
Identifying subgroup markers in heterogeneous populations. (J de Ronde, et al. 2013.)
}
\author{
  Julian de Ruiter wrote this wrapper function of the original DIDS code. DIDS itself was written by Jorma de Ronde.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
res <- didsScore(eset)
didsPlot(eset, gene='Fgfr2')
}

