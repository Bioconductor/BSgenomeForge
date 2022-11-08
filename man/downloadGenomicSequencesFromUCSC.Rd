\name{downloadGenomicSequencesFromUCSC}

\alias{downloadGenomicSequencesFromUCSC}

\title{The downloadGenomicSequencesFromUCSC function}

\description{
    A function for downloading UCSC genome .2bit sequence files
 }

\usage{
    downloadGenomicSequencesFromUCSC(genome, goldenPath.url=getOption("UCSC.goldenPath.url"), destdir=".", method, quiet=FALSE)
 }

\arguments{
  \item{genome}{
    This is the name of the UCSC genome sequence to be downloaded.
    It is used to form the download URL.
  }
  
  \item{goldenPath.url}{
    goldenPath.url: a ‘character’ string set to getOption("UCSC.goldenPath.url)
    by default.getOption("UCSC.goldenPath.url) return the goldenPath URL, 
    http://hgdownload.cse.ucsc.edu/goldenPath. This argument can be changed as needed. 
  }
  
  \item{destdir}{
    destdir: It is a 'character' string indicating the path where 
    the downloaded file is to be stored
    destdir is set to "." by default. The file is saved as ".(genome).2bit" 
    This file is hidden but can be seen from command line by typing 
    "ls -a" in the parent directory. 
    "cp '.(genome).2bit' 'genome.2bit'" can then be used to change 
    the name to genome.2bit which can then be accessed easily.    
  }
  
  \item{method}{
   Method to be used for downloading files.  Current download
   methods are ‘"internal"’, ‘"wininet"’ (Windows only)
   ‘"libcurl"’, ‘"wget"’ and ‘"curl"’, and there is a value
   ‘"auto"’: see ‘Details’ and ‘Note’.
   } 
    
  \item{quiet}{
   If ‘TRUE’, suppress status messages (if any), and the
   progress bar.
   }  
 }

\details{
  This function are intended for Bioconductor users who want 
  to download the genomic sequence a UCSC genome in .2bit file format.
  The function ‘downloadGenomicSequencesFromUCSC’ can be used to download 
  a single file as described by ‘goldenPath.url’ internet and store 
  it in ‘destdir’.
  
  genome is used to make the URL used in the download.file function.
  
  The ‘goldenPath.url’ must start with a scheme such as ‘http://’, 
  ‘https://’, ‘ftp://’ or ‘file://’.  Which methods support which schemes 
  varies by R version, but ‘method = "auto"’ will try to find a method
  which supports the scheme.  
 }
  
  \seealso{
  The \code{\link[utils]{download.file}} function in the \pkg{utils} package, that
  \code{downloadGenomicSequencesFromUCSC} uses internally to download the 2bit file.
 }
  
\examples{
   
   ## Download a genome.2bit file 
   ## Downloading a genome.2bit file with only the genome passed
   
   library(BSgenomeForge)
   downloadGenomicSequencesFromUCSC("sacCer1")
   
   ## Downloading a genome.2bit file with genome and destdir arguments passed
   
   library(BSgenomeForge)
   downloadGenomicSequencesFromUCSC("sacCer1", destdir="/home/emmanuel/")
   
   ## How to import the file in R as a DNAStringSet object.
   ## genome_sequence <- import("/genomePath/genome.2bit") 
   library(BSgenomeForge)  
   genome_sequence <- import("/home/emmanuel/sacCer1.2bit") 
   genome_sequence
 }
   
   
   
   
   
