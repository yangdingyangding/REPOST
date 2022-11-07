library("data.table")
library("magrittr")
library("Biostrings")
library("readxl")
library("ggpubr")



## 1. extract all mammalian sequences (from both Step1 and Step2) and align them again using MUSCLE

#### 1.1. get mammalian accessions, including those from NCBI datasets and others
parent.path <- "results/Step1/NCBI_datasets_download_ortholog_ACE2_human/"
protein.identity.info.with.species.category.dt <- fread(paste(sep="", parent.path, "/protein.identity.info.with.species.category.dt.csv"))
writeLines(protein.identity.info.with.species.category.dt[is.mammals.info == "Mammals", protein_accession] %>% sort, "./results/Step4/ACE2_mammals_only_protein_accessions.txt")

liu2020.tested.dt <- {
    read_excel("./external/papers/DOI/10.1073/pnas.2025373118/pnas.2025373118.sd01.xlsx", skip=3) -> .;
    data.table(.) ->.;
    setnames(., 1, "ID")
    ## remove untested
    .[ID %in% c(NA, "") == FALSE] -> .;
    ## label permissives and non-permissives
    .[, label:="permissive"]
    .[`Common Name` %in% c("white-tufted-ear marmoset", "Tufted capuchin", "Bolivian squirrel monkey", "Koala"), label:="non-permissive"]
    ## return
    .
}
## the followings are missing, so we manually downloaded them
#### > setdiff(liu2020.tested.dt[, `NCBI Accession`], protein.identity.info.with.species.category.dt[, protein_accession])
#### [1] "XP_010334925.1" "XP_019781177.1" "XP_015974412.1" "XP_002930657.1" "MT560518"       "XP_017505746.1"



#### 1.2. use these accessions to extract mammalian fasta and align them using MUSCLE
#### NOTE: MT506518 is a CDS sequence with N; therefore, its protein sequence cannot be determined and we excluded this accession.  

{
    readAAStringSet("results/Step1/NCBI_datasets_download_ortholog_ACE2_human/ncbi_dataset/data/protein.faa") -> temp.all.AAStringSet;
    readLines("results/Step4/ACE2_mammals_only_protein_accessions.txt") -> temp.mammal.accessions
    grep(pattern=paste(sep="", collapse="|", temp.mammal.accessions), x=names(temp.all.AAStringSet)) -> temp.hit.names;
    temp.all.AAStringSet[temp.hit.names] -> temp.hit.AAStringSet
    writeXStringSet(x=temp.hit.AAStringSet, filepath="results/Step4/ACE2_mammals_only_protein_sequences.fasta")
}

system("cat ./results/Step4/ACE2_mammals_only_protein_sequences.fasta ./results/Step2/XP_007989304.1.fasta ./results/Step2/AAX63775.1.fasta ./results/Step2/ACT66275.1.fasta ./results/Step4/XP_002930657.1.fasta ./results/Step4/XP_010334925.1.fasta ./results/Step4/XP_015974412.1.fasta ./results/Step4/XP_017505746.1.fasta ./results/Step4/XP_019781177.1.fasta  > ./results/Step4/ACE2_mammals_only_plus_experimentally_tested_protein_sequences.fasta")
system("./tools/muscle3.8.31/src/muscle -in ./results/Step4/ACE2_mammals_only_plus_experimentally_tested_protein_sequences.fasta -out ./results/Step4/ACE2_mammals_only_plus_experimentally_tested_protein_sequence_alignment.fasta")

#### 1.3. extract the human ACE2-based alignment (i.e., leaving no gaps in human ACE2)
protein.mammals.only.plus.experimentally.tested.muscle.alignment.AAStringSet <- readAAStringSet("./results/Step4/ACE2_mammals_only_plus_experimentally_tested_protein_sequence_alignment.fasta")
names(protein.mammals.only.plus.experimentally.tested.muscle.alignment.AAStringSet) <- sub(pattern="^([^ ]+) .*", replacement="\\1", x=names(protein.mammals.only.plus.experimentally.tested.muscle.alignment.AAStringSet))
protein.mammals.only.plus.experimentally.tested.muscle.alignment.matrix <- as.matrix(protein.mammals.only.plus.experimentally.tested.muscle.alignment.AAStringSet)
protein.mammals.only.plus.experimentally.tested.muscle.alignment.human.ACE2.based.matrix <- protein.mammals.only.plus.experimentally.tested.muscle.alignment.matrix %>% {.[, .["NP_068576.1", ] != "-" ]}
write.table(protein.mammals.only.plus.experimentally.tested.muscle.alignment.human.ACE2.based.matrix, "./results/Step4/ACE2_mammals_only_plus_experimentally_tested_protein_sequence_alignment_human_ACE2_based_matrix.txt", sep=",",quote=FALSE)

## 2. compute the weighted distance

#### 2.1. get the positive (permissive) and negative (non-permissive) accessions from the first round
first.round.dt <- fread("./results/Step2/SARS-CoV-2_first_round_info.csv")
positive.protein.accessions.vector <- first.round.dt[label=='permissive', accession]
negative.protein.accessions.vector <- first.round.dt[label=='non-permissive', accession]

#### 2.2. build species info table for later annotation
protein.identity.info.with.species.category.dt <- fread(paste(sep="", parent.path, "/protein.identity.info.with.species.category.dt.csv"))
protein.with.species.info.dt <- rbindlist(list(
    protein.identity.info.with.species.category.dt[, list(protein_accession, scientific_name, common_name, category.corrected, is.mammals.info) ],
    data.table(protein_accession="AAX63775.1", scientific_name="Paguma larvata", common_name="masked palm civet",  category.corrected="carnivores", is.mammals.info="Mammals"),
    data.table(protein_accession="XP_007989304.1", scientific_name="Chlorocebus sabaeus", common_name="green monkey",  category.corrected="primates", is.mammals.info="Mammals"),
    data.table(protein_accession="ACT66275.1", scientific_name="Rhinolophus sinicus", common_name="Chinese rufous horseshoe bat",  category.corrected="bats", is.mammals.info="Mammals"),
    data.table(protein_accession="XP_002930657.1", scientific_name="Ailuropoda melanoleuca", common_name="Giant panda",  category.corrected="carnivores", is.mammals.info="Mammals"),
    data.table(protein_accession="XP_010334925.1", scientific_name="Saimiri boliviensis boliviensis", common_name="Bolivian squirrel monkey",  category.corrected="primates", is.mammals.info="Mammals"),
    data.table(protein_accession="XP_015974412.1", scientific_name="Rousettus aegyptiacus", common_name="Egyptian rousette",  category.corrected="bats", is.mammals.info="Mammals"),
    data.table(protein_accession="XP_017505746.1", scientific_name="Manis javanica", common_name="Malayan pangolin",  category.corrected="other placental mammals", is.mammals.info="Mammals"),
    data.table(protein_accession="XP_019781177.1", scientific_name="Tursiops truncatus", common_name="Common bottlenose dolphin",  category.corrected="whales and dolphins", is.mammals.info="Mammals")
    ))
fwrite(protein.with.species.info.dt, "./results/Step4/protein.with.species.info.dt.csv")

#### 2.3. enumerate all species pairs and add their annotation
#### first.protein: the species to predict
#### second.protein: the neighbor of `first.protein`; will keep first round-tested species only later
protein.mammals.only.plus.experimentally.tested.muscle.alignment.human.ACE2.based.protein.accessions.vector <- rownames(protein.mammals.only.plus.experimentally.tested.muscle.alignment.human.ACE2.based.matrix)
protein.mammals.only.plus.experimentally.tested.muscle.alignment.human.ACE2.based.protein.accession.pairs.dt <- protein.mammals.only.plus.experimentally.tested.muscle.alignment.human.ACE2.based.protein.accessions.vector %>%
    ## enumerate all pairs
    {expand.grid(first.protein.accession=., second.protein.accession=., stringsAsFactors=FALSE)} %>%
    data.table %>%
    ## add species info for `first.protein`
    {merge(x=., y=protein.with.species.info.dt[, list(first.protein.accession=protein_accession, first.protein.species.scientific.name=scientific_name, first.protein.species.common.name=common_name, first.protein.species.category.corrected=category.corrected, first.protein.species.is.mammals.info=is.mammals.info)], by="first.protein.accession", all.x=TRUE, all.y=FALSE)} %>%
    ## add species info for `second.protein`
    {merge(x=., y=protein.with.species.info.dt[, list(second.protein.accession=protein_accession, second.protein.species.scientific.name=scientific_name, second.protein.species.common.name=common_name, second.protein.species.category.corrected=category.corrected, second.protein.species.is.mammals.info=is.mammals.info)], by="second.protein.accession", all.x=TRUE, all.y=FALSE)}

#### 2.4. compute the weighted distances, using only those positions that are negative(non-permissive)-informative
negative.informative.residue.positions.vector <- fread("./results/Step3/SARS-CoV-2_first_round_non-permissive_informative.dt.csv")[, position.on.human.ACE2] %>% sort %>% unique
weighted.distance.dt <- copy(protein.mammals.only.plus.experimentally.tested.muscle.alignment.human.ACE2.based.protein.accession.pairs.dt)[
  , weighted.distance := sum((
      protein.mammals.only.plus.experimentally.tested.muscle.alignment.human.ACE2.based.matrix[first.protein.accession, negative.informative.residue.positions.vector] !=
      protein.mammals.only.plus.experimentally.tested.muscle.alignment.human.ACE2.based.matrix[second.protein.accession, negative.informative.residue.positions.vector]
  )), list(first.protein.accession, second.protein.accession)]
fwrite(weighted.distance.dt, paste(sep="", "./results/Step4/weighted_distance.dt.csv.gz"))

## 3. compute R for each species to test
#### 3.1. keep only first round-tested species for the `second.protein`
weighted.distance.second.tested.only.dt <- weighted.distance.dt[second.protein.accession %in% c(positive.protein.accessions.vector, negative.protein.accessions.vector)][second.protein.accession %in% positive.protein.accessions.vector, second.exp.result:="positive"][second.protein.accession %in% negative.protein.accessions.vector, second.exp.result:="negative"]
#### 3.2. compute for each `first.protein` its R by finding its closest positive and negative species from the first round-tested species
#### NOTE that here we do NOT consider those first round-species matching the `first.protein`
prioritization.within.different.species.dt <- weighted.distance.second.tested.only.dt[first.protein.species.scientific.name != second.protein.species.scientific.name, list(R=log10(min(.SD[second.exp.result=='negative', weighted.distance])/min(.SD[second.exp.result=='positive', weighted.distance]))), list(first.protein.accession, first.protein.species.scientific.name, first.protein.species.common.name, first.protein.species.category.corrected, first.protein.species.is.mammals.info)]
#### 3.3. annotate each `first.protein` whether it has been tested in the first round
prioritization.within.different.species.dt[, is.in.first.round:=first.protein.accession %in% c(positive.protein.accessions.vector, negative.protein.accessions.vector)]
#### save
fwrite(prioritization.within.different.species.dt, "./results/Step4/prioritization.within.different.species.SC2.dt.csv")

## 4. test on validation dataset

#### 4.1. Liu2020 test
{

    liu2020.tested.dt[`NCBI Accession` %in% c("MT560518") == FALSE, list(protein.accession=`NCBI Accession`, species.scientific.name=Species, exp.result=label)] -> .;
    merge(x=prioritization.within.different.species.dt, y=., by.x="first.protein.accession", by.y="protein.accession", all.x=FALSE, all.y=TRUE) -> .;
    .[, predicted.result:=c("non-permissive", "permissive")[(R>0) + 1]]
    .[, prediction.matched.exp:=c("No", "Yes")[(predicted.result == exp.result) + 1]]
    melt(., id.vars=c("first.protein.accession", "first.protein.species.scientific.name", "first.protein.species.common.name", "first.protein.species.category.corrected", "is.in.first.round"), measure.vars=c("exp.result", "predicted.result", "prediction.matched.exp"), variable.name="y.name", value.name="value") -> .;
    .[, value.ordered:=factor(value, levels=c("permissive", "non-permissive", "Yes", "No"))]
    .[, value.text:=c("permissive"="P", "non-permissive"="N", "Yes"="Y", "No"="N")[value]]
    .[, y.name.prettified:=factor(c("exp.result"="exp.", "predicted.result"="pred.", "prediction.matched.exp"="matched?")[y.name], levels=c("exp.", "pred.", "matched?"))]
    .[, first.protein.species.common.name.prettified:=first.protein.species.common.name][first.protein.species.common.name %in% c(NA, ""), first.protein.species.common.name.prettified:=paste(sep="", "(", first.protein.species.scientific.name, ")")][first.protein.species.common.name == "human" & first.protein.accession == "NP_001358344.1", first.protein.species.common.name.prettified:="human (NP_001358344.1)"]
    .[, first.protein.species.scientific.name.prettified:=first.protein.species.scientific.name][first.protein.species.scientific.name %in% c(NA, ""), first.protein.species.scientific.name.prettified:=paste(sep="", "(", first.protein.species.scientific.name, ")")][first.protein.species.common.name == "human" & first.protein.accession == "NP_001358344.1", first.protein.species.scientific.name.prettified:="human (NP_001358344.1)"]
    .[, first.protein.species.category.corrected.prettified:=first.protein.species.category.corrected][first.protein.species.category.corrected %in% c("even-toed ungulates", "odd-toed ungulates") | first.protein.species.scientific.name == "Ceratotherium simum simum", first.protein.species.category.corrected.prettified:="ungulates"][first.protein.species.category.corrected == "marsupials", first.protein.species.category.corrected.prettified:="marsup."][first.protein.species.category.corrected %in% c("other placental mammals", "other placentals") & first.protein.species.scientific.name != "Ceratotherium simum simum", first.protein.species.category.corrected.prettified:="o."]
    . -> to.plot.dt
    ##
    ggplot(to.plot.dt[is.in.first.round==FALSE], aes(x=first.protein.species.scientific.name.prettified, y=y.name.prettified,  fill=value.ordered)) +
        geom_tile() +
        geom_text(aes(label=value.text, color=value.ordered)) +
        facet_grid(~first.protein.species.category.corrected.prettified, scales="free_x", space="free_x") +
        theme_pubr(base_size=10, legend="none") + theme(axis.text.x=element_text(angle=45, hjust=1)) +
        labs(x="", y="") +
        scale_color_manual(values=c("white", "black", "black", "white")) +
        scale_fill_manual(values=c("black", "grey70", "lightgreen", "red"))-> .;
    ggsave(filename="./results/Step4/Liu2020.test.validation.tile.png", plot=., device="png", width=22, height=8, units="cm")
        
}

#### 4.2. non-permissives validated by ourselves: manually examine those with R<0 in prioritization.within.different.species.dt

{
    fread("results/Step4/SARS-CoV-2.second.round.test.protein.info.csv") -> second.round.dt
    merge(x=prioritization.within.different.species.dt, y=second.round.dt[, list(first.protein.accession=protein.accession, exp.result)], by="first.protein.accession", all.x=FALSE, all.y=TRUE) -> .;
    fwrite(., "./results/Step4/SARS-CoV-2.second.round.test.prediction.and.exp.comparison.dt.csv")
}

#### Figure S1

{
    
    weighted.distance.second.tested.only.dt -> .;
    ##.[first.protein.species.scientific.name != second.protein.species.scientific.name] -> .;
    .[first.protein.accession %in% negative.protein.accessions.vector] -> .;
    .[, second.is.negative:=c("permissive", "non-permissive")[as.integer(second.protein.accession %in% negative.protein.accessions.vector) + 1] ] -> .;
    . -> to.plot.dt
    ##
    to.plot.dt -> .;
    ggplot(., aes(x=second.protein.species.scientific.name, y=weighted.distance, fill=second.is.negative)) -> .;
    . + geom_bar(stat="identity") -> .;
    . + coord_flip() -> .;
    . + facet_grid(~first.protein.species.scientific.name) -> .;
    . + theme_pubr() -> .;
    . + labs(x="Pilot orthologs", y="Weighted distance to the pilot ortholog", fill="") -> .;
    ggsave("./results/Step4/SF1.png", ., device="png", width=21, height=12, units="cm")
    
}
