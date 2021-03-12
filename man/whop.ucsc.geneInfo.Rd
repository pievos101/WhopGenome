\name{whop.ucsc.geneInfo}
\alias{whop.ucsc.geneInfo}
\title{
Return information from UCSC about a gene named precisely as specified
}
\description{
Information about a gene (and optionally required to be located on a certain chromosome) is returned.
}
\usage{
whop.ucsc.geneInfo(gen, chr = NA)
}
\arguments{
  \item{gen}{Gene name to query information about}
  \item{chr}{If specified, the identifier of the chromosome, on which this gene is located}
}
\details{
Gene name, chromosome of origin, strand, and start and end positions of transcription site, coding sequence and exons are returned.
}
\value{
  \item{geneName}{Gene name}
  \item{name}{Gene identifier}
  \item{chrom}{Chromosome, on which the gene is located}
  \item{strand}{Whether this gene is located on the + or - strand}
  \item{txStart}{Transcription start site}
  \item{txEnd}{Transcription end}
  \item{cdsStart}{Coding sequence start}
  \item{cdsEnd}{coding sequence end}
  \item{exonCount}{Number of exons of this gene}
  \item{exonStarts}{comma-separated list of exon start position}
  \item{exonEnds}{comma-separated list of exon end positions}
}
\author{
Ulrich Wittelsbuerger
}
\seealso{
whop.ucsc.geneInfoSimilar
}
