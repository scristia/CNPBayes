\name{granges}
\alias{granges}
\docType{data}
\title{
Genomic Ranges
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(granges)}
\format{
  The format is:
Formal class 'GRanges' [package "GenomicRanges"] with 6 slots
  ..@ seqnames       :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. ..@ values         : Factor w/ 22 levels "chr1","chr2",..: 1 2 3 4 5 6 7 8 9 10 ...
  .. .. ..@ lengths        : int [1:22] 28 31 32 31 30 22 27 22 20 20 ...
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
  .. .. ..@ start          : int [1:429] 51598 1617778 12764515 13247932 16024214 16741950 16878321 22190770 25455928 34875241 ...
  .. .. ..@ width          : int [1:429] 36868 44686 173568 9288 3384 117489 256514 19943 80004 11615 ...
  .. .. ..@ NAMES          : NULL
  .. .. ..@ elementType    : chr "integer"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ strand         :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 3
  .. .. ..@ lengths        : int 429
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ elementMetadata:Formal class 'DataFrame' [package "IRanges"] with 6 slots
  .. .. ..@ rownames       : NULL
  .. .. ..@ nrows          : int 429
  .. .. ..@ listData       :List of 1
  .. .. .. ..$ numberProbes: int [1:429] 18 19 47 15 11 50 129 11 35 21 ...
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ seqinfo        :Formal class 'Seqinfo' [package "GenomicRanges"] with 4 slots
  .. .. ..@ seqnames   : chr [1:22] "chr1" "chr2" "chr3" "chr4" ...
  .. .. ..@ seqlengths : int [1:22] 247249719 242951149 199501827 191273063 180857866 170899992 158821424 146274826 140273252 135374737 ...
  .. .. ..@ is_circular: logi [1:22] NA NA NA NA NA NA ...
  .. .. ..@ genome     : chr [1:22] "hg18" "hg18" "hg18" "hg18" ...
  ..@ metadata       :List of 1
  .. ..$ genome: chr "hg18"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(granges)
## maybe str(granges) ; plot(granges) ...
}
\keyword{datasets}
