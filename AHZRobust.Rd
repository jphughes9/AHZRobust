\name{AHZ.glmerMod}
\alias{AHZ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Approximate Hotelling (Zhang) degrees of freedom calculation for F-tests using a cluster-robust variance.
}
\description{
\code{AHZ} returns an F-statistic based on a cluster-robust variance, denominator degrees of freedom and p-value for glmer objects.
}
\usage{
AHZ(obj,L,cluster,type="classic",kadjust=FALSE,Fadjust=FALSE)
}

\arguments{
 \item{obj}{
Fitted glmer model 
}
\item{L}{
A q x p matrix of contrasts, where q is the number of contrasts and p is the number of fixed effects
}
\item{cluster}{
Optional expression or vector indicating which observations belong to the same cluster. If not specified, will be set to clubSandwich:::get_outer_group(obj)
}
\item{type}{
Character string specifying which small-sample adjustment should be used for the robust variance calculation, with available options "classic", "DF", "KC", "MD", "FG(\eqn{r_{fg}})", or "MBN(DF,d,\eqn{r_{mbn}})". See documentation for \code{vcovCR.glmerMod}. Default is "classic".
}
\item{kadjust}{
Should the constant \eqn{k} be set to 1 (kadjust=FALSE) or \eqn{\frac{\eta_Z - q + 1}{\eta_Z}} (kadjust=TRUE) when computing \eqn{\cal F} (see details)
}
\item{Fadjust}{
Should \eqn{F_j} be set equal to the identity matrix (Fadjust=FALSE) or the value used to compute \eqn{V^r} (Fadjust=TRUE) in computing \eqn{B_j} (see details)
}
}

\details{
Write the generalized linear mixed model as 
\deqn{E(Y | \gamma) = g^{-1}(X \beta + Z \gamma) = g^{-1}(\eta)}
\deqn{Var(Y | \gamma) = \Sigma}
where \eqn{Y} is \eqn{n \times 1}, \eqn{X} is \eqn{n \times p}, \eqn{Z} is \eqn{n \times q}, and \eqn{\gamma \sim N(0,R(\theta))}. \eqn{g} is a link function; only identity, log and logit link are accepted at this time.

Define
\deqn{\Delta = \left( \frac{\partial g^{-1}(\eta)}{\partial \eta} \right)_{\beta, \gamma} \mbox{a diagonal matrix}}
\deqn{P = \Delta^{-1} (Y - g^{-1}(\eta)) + X \beta + Z \gamma}
\deqn{V =  Z R(\theta) Z^{\prime} + \Delta^{-1} \Sigma \Delta^{-1}}
\deqn{e_i = P_i - X_i \beta }
Note that \eqn{M = (X^{\prime} V^{-1} X)^{-1}} is the usual model-based variance.

The robust variance of \eqn{\hat{\beta}} is given by the expression
\deqn{
V^r = c \times M \left( \sum_g A_g X_g^{\prime} V_g^{-1} F_g^{\prime} e_g e_g^{\prime} F_g V_g^{-1} X_g A_g \right) M + \delta_m \phi M}
where the sum is over the \eqn{m} independent groups (clusters) and all terms are evaluated at the estimates of \eqn{\hat{\beta}, \hat{\theta}, \mbox{and } \hat{\Sigma}}. 
Different values of c, F, A, \eqn{\delta_m} and \eqn{\phi} lead to different versions of the robust variance - see the documentation for \code{vcovCR.glmerMod}.

Consider the null hypothesis
\deqn{
Ho: L \beta = 0
}
where \eqn{L} is a \eqn{q \times p} matrix defining \eqn{q} independent contrasts. 
Let
\deqn{
{\cal F} = k*\hat{\beta}^\prime L^\prime (L V^r L^\prime)^{-1} L \hat{\beta}/q .
}
We assume \eqn{\cal F} is distributed as according to an F-distribution with \eqn{q} 
and \eqn{ddf} degrees of freedom. We use the 
 "Approximate Hotelling (Zhang)" (AHZ) (Tipton and Pustejovsky, 2015) method
for calculating the denominator degrees of freedom.

Define
\deqn{
Var(d_{st}) = \sum_{i=1}^m \sum_{j=1}^m u_s^{\prime} B_i V B_j^{\prime} u_t u_t^{\prime} B_i V B_j^{\prime} u_s + u_s^{\prime} B_i V B_j^{\prime} u_s u_t^{\prime} B_i V B_j^{\prime} u_t
}
where the q-vector \eqn{u_s} has the \eqn{s}'th entry equal to 1 and all other entries 0, and
\deqn{
B_j = \Gamma^{-1/2} L M X_j W_j F_j (I_K - XMX^{\prime} W)_j
}
where \eqn{\Gamma = LML^{\prime}} and \eqn{W = V^{-1}}.
Then compute 
\deqn{
\eta_Z = \frac{q(q+1)}{\sum_{s=1}^q \sum_{t=1}^q Var(d_{st})}
}
The denominator degrees of freedom are then calculated as \eqn{ddf = \eta_Z - q + 1} and
p-values are computed under the assumption that \eqn{{\cal F} \sim F(q, ddf)}.
}

\value{
For 1-degree of freedom tests, a list giving
\itemize{
\item{estimate - the estimated value of the contrast}
\item{se - the standard error of the estimated contrast}
\item{t stat - the t-statistic for the contrast}
\item{df - the degrees of freedom for the t-statistic}
\item{pvalue - the p-value}
\item{robustVar - the robust variance matrix for the fitted model}}

For multi-degree of freedom tests, a list giving
\itemize{
\item{F value - the F statistic}
\item{ndf - the numerator degrees of freedom for the F-statistic}
\item{ddf - the denominator degrees of freedom for the F-statistic}
\item{pvalue - the p-value}
\item{robustVar - the robust variance matrix for the fitted model}}
}

\references{
Tipton E and Pustejovsky JE. Small sample adjustments for tests of 
moderators and model fit using robust variance estimation in meta-regression.
 J of Educational and Behavioral Statistics 40:604-634, 2015.

}
\author{
James P Hughes
}
%%\note{
%%  ~~further notes~~
%%}

\examples{
L=matrix(c(0,0,0,0,0,1,0,0,-1,
           0,0,0,0,0,1,0,-1,0,
           0,0,0,0,0,1,-1,0,0),3,9,byrow=TRUE)
rslt = glmer(ct ~ ftime + ftot + (1 | fcluster), family=binomial, data=ctdata)
AHZ(rslt,L,type="MBN")
}
\keyword{ robust variance }
\keyword{ cluster randomized trial }
\keyword{ F test }
\keyword{ degrees of freedom }% __ONLY ONE__ keyword per line
