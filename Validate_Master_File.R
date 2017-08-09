#!/usr/bin/env Rscript

#load the different tests in validations.R file
source("validations.R")

## Get User Input:
args=(commandArgs(TRUE))

cat("Command line arguments: ",args, "\n\n")

if(length(args)==0) {
    usage()
    quit()
}

for(i in 1:length(args)) {
    eval(parse(text=args[i]))
}

verifyFileExists(design.file)

# Read file to validate
x <- read.table(design.file, header=TRUE, comment.char="", sep="\t", stringsAsFactors=FALSE)

###########
## generalChecking
###########
checkColumns(x)                #test_files: TEST_file__numCols_check.txt
checkHeader(x)                #test_files:TEST_file__column_messup.txt

errorVal <- 0
errorValLater <- 0
checkInf(x)                #test_files: TEST_file__inf.txt
checkNA(x)                #test_files: TEST_file__NA.txt

if(errorVal != 0){
    cat("    Cannot continue validation until the above errors are fixed.\n")
    q()
}

checkColClasses(x)                #test_files: TEST_file__TilingABCD.txt, 


# ERROR variable. If there was an error and it was printed, change the value of errorVal=1 in function, so I can exit after most of the validations
errorVal <- 0

## Go column by column, and make sure each column has the right formatting.
checkGeneID(x)                #test_files: TEST_file__geneID.txt
checkTargetID(x)                #test_files: TEST_file__targetID.txt
checkIntervals(x)                #test_files: TEST_file__interval_format.txt
checkForNumber(x, "Length")                #test_files: TEST_file__LengthABCD.txt
checkIntLen_matching(x)                #test_files: TEST_file__intervals.txt, TEST_file__lengthMatch.txt, TEST_file__wrongCHR.txt
checkForNumber(x, "Tiling")                #test_files: TEST_file__TilingABCD.txt

## Make the Categories LVLs, then check
x$Category <- factor(x$Category)
checkCategory(x$Category)                #test_files: TEST_file__category_change.txt


## Specific Checks

## Category Fingerprint has GeneID "FP"
#checkGeneIDbyType(x,"Fingerprint", "rs")                #test_files: TEST_file__FP_geneid.txt

## Category Tiling has GeneID "Tiling"
#checkGeneIDbyType(x,"Tiling", "Tiling")                #test_files: TEST_file__tiling_geneid.txt

## TargetIDs with "_rs"### have a length of 1, AND gene name has to have rsid in it                #test_files: TEST_file__rsID_len.txt
checkSNPlen(x)
checkTargetUnique(x)


if (errorVal == 0 && errorValLater == 0)
{
    cat("If you've made it this far, it means that the file seems to pass all our validation criteria!\n")
}

