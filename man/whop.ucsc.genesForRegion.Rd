\name{whop.ucsc.genesForRegion}
\alias{whop.ucsc.genesForRegion}
\title{
Return a list of genes located in a certain region on a certain chromosome
}
\description{
Details on all genes falling into the positions between 'beg' and 'end' on chromosome 'chrom' are returned.
}
\usage{
whop.ucsc.genesForRegion(chrom, beg, end)
}
\arguments{
  \item{chrom}{Chromosome on which to look in "chr1" notation}
  \item{beg}{First position of the region a gene may fall into}
  \item{end}{Last position of the region a gene may fall into}
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
