\name{whop.eg.selectOrganism}
\alias{whop.eg.selectOrganism}
\title{
Select the organism to query with subsequent whop.eg calls and load the appropiate database(s).
}
\description{
Select the organism to query with subsequent whop.eg calls and load the appropiate database(s).
}
\usage{
whop.eg.selectOrganism(organismname, dontload = FALSE, install.if.missing = F)
}
\arguments{
  \item{organismname}{Organism to query}
  \item{dontload}{Whether to load the database}
  \item{install.if.missing}{Whether to install the database, if it does not exist locally}
}
\value{
Success status
}
\author{
Ulrich Wittelsbuerger
}
