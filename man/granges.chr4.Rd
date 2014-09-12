\name{granges.chr4}
\alias{granges.chr4}
\docType{data}
\title{
Regions of interest
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(granges.chr4)}
\format{
  The format is:
Formal class 'GRanges' [package "GenomicRanges"] with 6 slots
  ..@ seqnames       :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. ..@ values         : Factor w/ 22 levels "chr1","chr2",..: 4
  .. .. ..@ lengths        : int 2
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
  .. .. ..@ start          : int [1:2] 9832514 10002252
  .. .. ..@ width          : int [1:2] 11853 7515
  .. .. ..@ NAMES          : NULL
  .. .. ..@ elementType    : chr "integer"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ strand         :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 3
  .. .. ..@ lengths        : int 2
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ elementMetadata:Formal class 'DataFrame' [package "IRanges"] with 6 slots
  .. .. ..@ rownames       : NULL
  .. .. ..@ nrows          : int 2
  .. .. ..@ listData       : Named list()
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
data(granges.chr4)
## maybe str(granges.chr4) ; plot(granges.chr4) ...
}
\keyword{datasets}
