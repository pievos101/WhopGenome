\name{whop.ucsc.geneInfoSimilar}
\alias{whop.ucsc.geneInfoSimilar}
\title{
Return information UCSC has about any genes with similar names
}
\description{
Information about any genes named similarly as specified in 'gen' (and optionally required to be located on chromosome 'chr') is returned.
}
\usage{
whop.ucsc.geneInfoSimilar(gen, chr = NA)
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
whop.ucsc.geneInfo
}
