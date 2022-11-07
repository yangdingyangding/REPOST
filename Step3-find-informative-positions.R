library("Biostrings")
library("data.table")
library("magrittr")
library("UpSetR")

alignment.AAStringSet <- readAAStringSet("./results/Step3/SARS-CoV-2_first_round_sequences_alignment.fasta")
names(alignment.AAStringSet) <- sub(pattern="^([^ ]+) .*", replacement="\\1", x=names(alignment.AAStringSet))

## must be of the same length, because it is an alignment
alignment.length <- width(alignment.AAStringSet)[1]

info.dt <- fread("./results/Step2/SARS-CoV-2_first_round_info.csv")

alignment.matrix <- as.matrix(alignment.AAStringSet)
alignment.permissive.only.matrix <- alignment.matrix[info.dt[label=="permissive", accession], ]
alignment.nonpermissive.only.matrix <- alignment.matrix[info.dt[label=="non-permissive", accession], ]

all.positions.vector <- 1:alignment.length
all.positions.on.human.ACE2.vector <- alignment.permissive.only.matrix["NP_068576.1", ] %>% {temp <- cumsum(. != "-"); temp[. == "-"] <- NA; temp}

comparison.dt <- rbindlist(mapply(all.positions.vector, all.positions.on.human.ACE2.vector, FUN=function(temp.position, temp.position.on.human.ACE2){
    return(data.table(
        nonpermissive.sequence.protein_accession=rownames(alignment.nonpermissive.only.matrix),
        position=temp.position,
        position.on.human.ACE2=temp.position.on.human.ACE2,
        nonpermissive.sequence.residue=alignment.nonpermissive.only.matrix[, temp.position],
        all.possible.positive.sequence.residues=alignment.permissive.only.matrix[, temp.position] %>% unique %>% sort %>% paste(collapse=","),
        is.in.positive.sequences=alignment.nonpermissive.only.matrix[, temp.position] %in% alignment.permissive.only.matrix[, temp.position]
    ))
}, SIMPLIFY=FALSE))

nonpermissive.informative.dt <- merge(x=comparison.dt[is.in.positive.sequences==FALSE & is.na(position.on.human.ACE2) == FALSE], y=info.dt[, list(accession, species, `common name`)], by.x="nonpermissive.sequence.protein_accession", by.y="accession", all.x=TRUE, all.y=FALSE)
## 193 residues

fwrite(nonpermissive.informative.dt, "./results/Step3/SARS-CoV-2_first_round_non-permissive_informative.dt.csv")

png("./results/Step3/first_round_nonpermissive_informative_position_overlap.png", width=13, height=6, units="cm", res=600)
print(upset(dcast(data.table(nonpermissive.informative.dt, temp=1), position~species, value.var="temp", fill=0), nsets=4, mb.ratio=c(0.7, 0.3), mainbar.y.label="# intersected\ninformative positions", sets.x.label="# informative positions"))
dev.off()
