\name{didsPlot}
\alias{didsPlot}
\title{
Detection of Imbalanced Differential Signal (DIDS) - Plot 
}
\description{
Create a DIDS plot for the specified gene.
}
\usage{
didsPlot(eset, gene='Fgfr2')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
ExpressionSet object for which the plot should be made. The featureData of the object should either contain a column `symbol' or another specified column (identified) by the geneSymbolCol parameter) that contains the gene symbol of each gene. The phenoData of the object should contain a `condition' column that marks each sample as being either `sensitive' or `resistant'. Note that this column may be ommitted if alternative grouping is provided by the groups parameter.
}
\item{gene}{
    Specifies which gene should be plotted (given as a gene symbol).
}
\item{geneSymbolCol}{
Specifies an alternative column in which gene symbol information can be found.
}
\item{groups}{
Specifies an alternative grouping of samples into two groups (overrides the condition column if present). Should be provided as a logical vector.
}
\item{groupNames}{
    Specifies the names of the two groups. Defaults to `sensitive' and `resistant'.
}
\item{densPlot}{
Specifies if the density plot should be created.
}
\item{binclass}{
    
}
\item{sorted}{
    
}
\item{main}{
    Title of the plot.
}
\item{ylim}{
    Y-limits of the plot.
}
\item{ylab}{
    Y-label of the plot.
}
\item{cols}{
    
}
\item{colCode}{
    
}
\item{minkw}{
}
\item{selection}{
    
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Creates a DIDS plot of the expression values of the specified gene.
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

