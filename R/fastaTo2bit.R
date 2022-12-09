### The function helps to convert a file from FASTA to 2bit
fastaTo2bit <- function(origfile, destfile)
  
{
    myfile = readDNAStringSet(origfile)
    export.2bit(myfile, destfile)
}

