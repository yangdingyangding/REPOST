library("reutils")
library("data.table")


## protein sequences
{

    c("XP_004758943.1") -> SPECIAL.IDS
    
    "./SARS-CoV-2_CR22.csv" -> .;
    fread(.) -> .;
    .[Abbreviation != "Brown Trout"] -> .;
    esearch(term=.[(`Protein sequence ID` %in% SPECIAL.IDS) == FALSE, paste(collapse=" OR ", `Protein sequence ID`)], db="protein") -> .;
    efetch(., db="protein", rettype="fasta", retmode="text") -> .;
    write(content(.), file="efetched.mammals.protein.sequences.fasta")

}
##
writeLines(c(
    readLines("./efetched.mammals.protein.sequences.fasta"),
    readLines("./XP_004758943.1.fasta")
    ), "./all.17.mammals.protein.sequences.fasta")

