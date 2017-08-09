border="################################################"

usage <- function(){
    usage.str= "\nUsage: Validate_Master_File.R  \"design.file='FileToTest.txt'\"
    OPTION:
    \"design.file='[REQUIRED: the design file to validate]'\"

    DESCRIPTION OF HOW TO MAKE A MASTER FILE FORMAT FILE:

    File has 6 columns in this order: GeneID, TargetID, Interval, Length, Category, Tiling. 
    GeneID can have letters, numbers, underscore, dashes, periods, semi-colon and '/', no spaces.
    TargetID can have letters, numbers, underscores, dashes, periods, semi-colon, and '/'. No spaces
    Interval must be in one base, and in this format: chr#:###-### The end must be a larger number than the start
    Length must be the interval's stop - start + 1.
    Category can only be: Exon, Intron, PseudoExon, Fingerprint, Promoter, MSI, Tiling, and Bait. Tiling and Fingerprint
        categories are interchangable.
    Tiling column can be a number representing the desired weights of the target in the pool. Otherwise, it is set to 1.

    Additional rules:
	- Anything with an rsID in the TargetID is considered a SNP and must be 1bp in length
        - Tiling or FP with rsIDs must have the rsIDs in the GeneID column
	- There must be at least one record with Category 'Bait' only if you got a 'baits' file  
        - IF this is a custom design ( not just spike in genes ) you must have a least one Fingerprint Record
    \n\n"
    
    cat(usage.str)
}

# Initialize Variables:
design.file=NULL
header=c("GeneID", "TargetID", "Interval", "Length", "Category", "Tiling")
classesCheck=c("character", "character", "character", "integer", "character", "integer")
categoryList=c("Exon", "Intron", "MSI", "PseudoExon", "Fingerprint", "Promoter", "Tiling", "Bait")


# File to be validated should exist
verifyFileExists <- function(this.file) {
    if(!(file.exists(this.file))) {
        usage()
        cat("ERROR: This file either doesn't exist or hasn't been put into the program correctly:", this.file, "Try again.\n")
        q()
    }
}

# Make sure there are the right number of columns
checkColumns <- function(x) {
    if(length(names(x)) != 6) {
        cat("ERROR: File should have 6 columns, this file does not!\n\n")
        cat("    If you need more information, start this script with no command line arguments.\n")
        cat("    Cannot continue validation script. Fix these first, then rerun.\n")
        q()
    }
}

# Verify Header is correct
checkHeader <- function(x) {
    if(!(all(names(x) == header))) {
        cat("ERROR: The file header is wrong, it should be 'GeneID', 'TargetID', 'Interval', 'Length', 'Category', 'Tiling' in that order.\n")
        cat("    Cannot continue validation script. Fix these first, then rerun.\n")
        q()
    }
}

# Check for infinite
checkInf <- function(x) {
    if(any(sapply(x, is.infinite))) {
        cat("ERROR: There is a value that is infinite in this file\n")
        errorVal <<- 1
    }
}

# Check for NaN
checkNA <- function(x) {
	for (i in 1:(length(x)-1)) {
		if(any(is.na(x[i]))){
			cat(paste0("ERROR: Column ", names(x)[i], " has an NA (is Empty) in row number ", which(is.na(x[i])), ".\n" ))
			errorVal <<- 1
		}
	}
}

#Verify the Column Classes are correct
checkColClasses <- function(x) {
	classesX <- sapply(x, class)
	if(!(all(classesX == classesCheck))) {
		badCols <- which(!(classesX == classesCheck))
		cat("WARNING/ERROR: Column(s): \n", paste(names(classesX[badCols]), collapse="\n"), " are not formatted right.\nHopefully the rest of the script will tell you why.\n", sep="")
	}
}

###########
###########
## Column Checks
###########
###########

# Makes sure that the length of the interval is the corresponding length 
intervals_length_verification <- function(i){
    result_string=""
    interval=i[1]
    regionLength= as.numeric(i[2])
    pieces = unlist(strsplit(interval, c(":|-"), perl=TRUE))
    chr=gsub("chr","", pieces[1])
    if((chr != "X") & (chr != "Y") & (chr != "M")) {
        chrN=as.numeric(chr)
        if(chrN >= 23){
            result_string=paste0(result_string, paste("ERROR: This chromosome number looks wrong: ", chr, "\n", sep=" "), collapse="\n")
        }
    }
    
    start=as.numeric(pieces[2])
    end=as.numeric(pieces[3])
    if(!(end >= start)){
        result_string=paste0(result_string, paste("ERROR: Chromosome start position is bigger than end position! Start position must be smaller than end position: ", interval,"\n", sep=" "), collapse="\n")
    }
    
    verifyLen=(end-start) + 1
    if(verifyLen != regionLength) {
        result_string=paste0(result_string, paste("ERROR: The length in the Length column is different than the actual length of the intervals.", interval,"should be", verifyLen, "not", regionLength,"\n", sep=" "), collapse="\n")
    }
    
    result_string
}

# GeneID format checking
checkGeneID <- function(x) {
	if(any(!(grepl("^[a-z|A-Z|0-9|/|\\.|\\-|_|;]+$", x$GeneID, perl=TRUE)))) {
		badVals = x[which(!(grepl("^[a-z|A-Z|0-9|\\.|\\-|/|_|;]+$", x$GeneID, perl=TRUE))),]
		cat(border, "\nERROR: GeneIDs are formatted incorrectly:\n")
                cat(paste(header, collapse="\t"), "\n", paste(apply(badVals,1,paste,collapse="\t"), collapse="\n"), sep="")
                cat("\nGeneIDs must be [a-z|A-Z|0-9|.|-|_|;] only.\n")
		errorVal <<- 1
	}
}

# TargetID can have letters, numbers, ., ;, _ and -. 
checkTargetID <- function(x) {
	if(any(!(grepl("^[a-z|A-Z|0-9|\\.|\\-|/|_|;]+$", x$TargetID, perl=TRUE)))) {
		badVals = x[which(!(grepl("^[a-z|A-Z|0-9|\\.|\\-|/|_|;]+$", x$TargetID, perl=TRUE))),]
		cat(border, "\nERROR: TargetIDs are formatted incorrectly:\n")
                cat(paste(header, collapse="\t"), "\n",   paste(apply(badVals,1,paste, collapse="\t"), collapse="\n"), sep="")
                cat("\nTargetIDs must be [a-z|A-Z|0-9|.|-|_|;]  only.\n")
		errorVal <<- 1
	}
}

# First do general check that the interval is in the right format
checkIntervals <- function(x) {
	if(any(!(grepl("^chr([0-9]{1,2}|X|M|Y):[0-9]+-[0-9]+$", x$Interval, perl=TRUE)))) {
		badVals = x[which(!(grepl("^chr([0-9]{1,2}|X|Y):[0-9]+-[0-9]+$", x$Interval, perl=TRUE))),]
		cat("ERROR: The following records have intervals that are incorrectly formatted: \n")
                cat(paste(header, collapse="\t"), "\n",  paste(apply(badVals,1,paste, collapse="\t"),collapse="\n"), "\n\nIntervals must be 'chr#:#-#'.\n", sep="")
		cat("    Cannot continue validation script. Fix these first, then rerun.\n")
                q()
	}
}

# Then Make sure that Length is a number
checkForNumber <- function(x, colName) {
	if(any(!(grepl("[0-9]+", x[,match(colName, names(x))], perl=TRUE)))) {
		badVals = x[which(!(grepl("[0-9]+", x[,match(colName, names(x))] , perl=TRUE))),]
		cat(border, "\nERROR: In the", colName,"column, there aren't just numbers ?!?!?!?!?!\n")
                cat(paste(header, collapse="\t"), "\n", paste(apply(badVals,1,paste, collapse="\t"), collapse="\n"), "\n\n", sep="")
                cat("    Cannot continue validation script. Fix these first, then rerun.\n")
                q()
	}
}

# Then make sure that the interval makes sense and the length is the (end - start +1) of the interval.
checkIntLen_matching <- function(x) {
    z=cbind(Interval=x$Interval, Length=x$Length)
    res=(apply(z,1,intervals_length_verification))
    res=paste0(res, sep="", collapse="")
    if(nchar(res) != 0)
    {
        cat(border, "\n" , res, sep="")
        errorVal <<- 1
    }
}


checkCategory <- function(category) {
    lvls = levels(category)
    badVals = unique(subset(category, !(category %in% categoryList)))
    if (length(badVals) > 0) {
        cat(border, "\nERROR: The categories:", paste(badVals, collapse="\n"), "are not recognized as real categories. If they need to be added, contact BIC.\n")
        errorVal <<- 1
    }

    if (!("Bait" %in% category)) {
        cat(border, "\nERROR: There are no baits in this file! You must provide BAITS.\n")
        errorValLater <<- 1
    }

    if (!("Fingerprint" %in% category) && !("Tiling" %in% category)) {
        cat(border, "\nWARNING: If this is a new custom design you MUST have at least one Fingerprint record, or the pipeline will not work correctly!\nIf these are spike in genes, ignore this warning.\n\n")
    }
}



###########
###########
## Specific Checks
###########
###########

## Fingerprint category has GeneID "FP"
## Tiling category has GeneID "Tiling"
checkGeneIDbyType <- function(x, categ, geneName) {
    if(categ %in% x$Category) {
        FP.regions <- x[which(x$Category == categ),]
        if(any(FP.regions$GeneID != geneName)) {
            badVals = FP.regions[which(FP.regions$GeneID != geneName),c("TargetID","GeneID")]
            badVals2 = apply(badVals, 1, paste, collapse=" Gene Name: ")
            cat(border, "\nERROR: These ",categ," regions do not have GeneID '",geneName,"'\n", paste0(badVals2, collapse="\n"), "\nThis is necessary. Change the Category of these regions, or rename the GeneID to '",geneName,"'.\n", sep="")
            errorVal <<- 1
        }
    }
}

## TargetIDs with "_rs"### have a length of 1
checkSNPlen <- function(x) {
    if(any(grepl("_rs[0-9]+$", x$TargetID, perl=TRUE)) || any(x$Category == "Fingerprint")  || any(x$Category == "Tiling")) {
        snps = x[which(grepl("_rs[0-9]+$", x$TargetID, perl=TRUE) | x$Category == "Fingerprint" | x$Category == "Tiling") ,]
        if(any(snps$Length != 1)){
            badVals = snps[which(snps$Length!=1),]$TargetID
            cat(border, "\nERROR: These targetIDs are either SNPs or fingerprint regions that are not length 1:\n")
            cat(paste0(badVals, collapse=", "))
            cat("\nIf they have rsIDs in the TargedID or if the category is Fingerprint, it must have length 1.\n")
            errorVal <<- 1
        }
        if(any(!grepl("rs[0-9]+$", snps$GeneID))){
            badVals = snps[which(!grepl("rs[0-9]+", snps$GeneID)),]
            cat(border, "\nWARNING: These records are SNPs that do not have rsIDs in the geneID column. If these are to be used for genotyping, it MUST have a real rsID in the GeneID column:\n")
            cat(paste(apply(badVals,1,paste,collapse="\t"), collapse="\n"))
            cat("\n\nWARNING: We need the rsID to be in geneID for Genotyping. If you are using this for copy number and NOT Genotyping, ignore this warning.\n\n")
            errorVal <<- 1
        }
    }
}

## TargetIDs must be UNIQUE. 
checkTargetUnique <- function(x) {
    if(anyDuplicated(x$TargetID)){
        dup <- x[which(duplicated(x$TargetID)),]$TargetID
        cat(border, "\nERROR: Same TargetdIDs used in multiple records:\n", paste0(dup, collapse="\n"), "\n", sep="")
        cat("ERROR: TargetIDs must be unique\n\n")
        errorVal <<- 1
    }
}

