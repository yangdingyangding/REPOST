library("Biostrings")
library("magrittr")
library("data.table")
library("XML")
library("reutils")
library("ggpubr")
library("ggsci")

parent.path <- "results/Step1/NCBI_datasets_download_ortholog_ACE2_human/"

protein.muscle.alignment.AAStringSet <- readAAStringSet(paste(sep="", parent.path, "/protein.muscle.alignment.fasta"))
names(protein.muscle.alignment.AAStringSet) <- sub(pattern="^([^ ]+) .*", replacement="\\1", x=names(protein.muscle.alignment.AAStringSet))


protein.muscle.alignment.matrix <- as.matrix(protein.muscle.alignment.AAStringSet)

protein.muscle.alignment.human.ACE2.vector <- protein.muscle.alignment.matrix["NP_068576.1", ]

protein.identity.info.dt <- apply(protein.muscle.alignment.matrix, MARGIN=1, FUN=function(temp.vector){
    sum((temp.vector == protein.muscle.alignment.human.ACE2.vector) & (temp.vector != "-"))
}) %>%
    {data.table(protein_accession=names(.), identical_residue_counts_with_human_ACE2=.)} %>%
    {.[, identity_with_human_ACE2:=identical_residue_counts_with_human_ACE2/sum(protein.muscle.alignment.human.ACE2.vector != '-')]}


tax.core.info.dt <- fread(paste(sep="", parent.path, "/tax.core.info.dt.csv"))


tax.core.info.dcast.dt <- dcast(tax.core.info.dt, tax_id + common.name ~category, value.var="temp")

tax.core.info.dcast.dt[primates == 1, category.corrected:="primates"]
tax.core.info.dcast.dt[rodents == 1, category.corrected:="rodents"]
tax.core.info.dcast.dt[carnivores == 1, category.corrected:="carnivores"]
tax.core.info.dcast.dt[`whales and dolphins` == 1, category.corrected:="whales and dolphins"]
tax.core.info.dcast.dt[bats == 1, category.corrected:="bats"]
tax.core.info.dcast.dt[`even-toed ungulates` == 1 & is.na(category.corrected), category.corrected:="even-toed ungulates"]
tax.core.info.dcast.dt[marsupials == 1, category.corrected:="marsupials"]
tax.core.info.dcast.dt[mammals == 1 & is.na(category.corrected), category.corrected:="other placentals"]
tax.core.info.dcast.dt[birds == 1, category.corrected:="birds"]
tax.core.info.dcast.dt[amphibians == 1, category.corrected:="amphibians"]
tax.core.info.dcast.dt[amniotes == 1 & is.na(category.corrected), category.corrected:="reptiles"]
tax.core.info.dcast.dt[`bony fishes and tetrapods` == 1 & is.na(category.corrected), category.corrected:="bony fishes"]
tax.core.info.dcast.dt[chordates == 1 & is.na(category.corrected), category.corrected:="other chordates"]

tax.core.info.dcast.supplemented.dt <- rbindlist(list(
    tax.core.info.dcast.dt,
    data.table(tax_id=c(8673, 260995, 1415580, 1825980), common.name=c(as.character(NA), as.character(NA), as.character(NA), "Goodes thornscrub tortoise"), category.corrected=c("reptiles", "amphibians", "amphibians", "reptiles"))
), use.names=TRUE, fill=TRUE)

fwrite(tax.core.info.dcast.supplemented.dt, paste(sep="", parent.path, "/tax.core.info.dcast.supplemented.dt.csv"))
## Plot Figure 2A
{
    
    tax.core.info.dcast.supplemented.dt -> A;
    A[, list(tax_id, category.corrected, mammals)] -> A;
    A[mammals==1, rough.category.corrected:="mammals"] -> A;
    A[is.na(mammals)==TRUE, rough.category.corrected:=category.corrected] -> A;
    A[, list(count=.N), list(rough.category.corrected.prettified=factor(rough.category.corrected, levels=c("mammals", "birds", "reptiles", "amphibians", "bony fishes", "other chordates")))] -> Fig.2A.left.to.plot.dt;
    A[rough.category.corrected == "mammals", list(count=.N), list(category.corrected.prettified=factor(category.corrected, levels=c("primates", "rodents", "carnivores", "whales and dolphins", "bats", "even-toed ungulates", "marsupials", "other placentals")))] -> Fig.2A.right.to.plot.dt;
    ##
    Fig.2A.left.to.plot.dt -> A;
    ggpie(A, "count", label=A[, count], lab.pos="out", lab.font="black", fill="rough.category.corrected.prettified", color="white") -> A;
    A + theme(legend.position="right", legend.justification=c(0,0.4)) + labs(fill="") -> A;
    A -> Fig.2A.left.ggplot
    ##
    Fig.2A.right.to.plot.dt -> A;
    ggpie(A, "count", label=A[, count], lab.pos="out", lab.font="black", fill="category.corrected.prettified", color="white") -> A;
    A + theme(legend.position = "right", legend.justification=c(0,0.4)) + labs(fill="") -> A;
    A -> Fig.2A.right.ggplot
    ##
    ggarrange(plotlist=list(Fig.2A.left.ggplot, Fig.2A.right.ggplot), ncol=1, align="v") -> A;
    ggsave(paste(sep="", parent.path, "/Fig2A.png"), plot=A, device="png", width=20/2, height=29/2.5*0.9, units="cm")
    
}
    


## merge the species + category info

ortholog.info.dt <- fread(paste(sep="", parent.path, "/ncbi_dataset/data/data_table.tsv"))

protein.identity.info.with.species.category.dt <- merge(
    x=protein.identity.info.dt,
    y=merge(x=ortholog.info.dt[, tax_id.character:=as.character(tax_id)], y=tax.core.info.dcast.supplemented.dt[, list(tax_id.character=as.character(tax_id), category.corrected)], by.x="tax_id.character", by.y="tax_id.character", all.x=TRUE, all.y=FALSE),
    by="protein_accession", all.x=TRUE, all.y=FALSE)
protein.identity.info.with.species.category.dt[category.corrected %in% c("primates", "rodents", "whales and dolphins", "carnivores", "bats", "other placentals","even-toed ungulates", "marsupials"), is.mammals.info:="Mammals"]
protein.identity.info.with.species.category.dt[category.corrected %in% c("birds", "reptiles", "amphibians", "bony fishes", "other chordates"), is.mammals.info:="Others"]

fwrite(protein.identity.info.with.species.category.dt, paste(sep="", parent.path, "/protein.identity.info.with.species.category.dt.csv"))

## Plot Figure 2B

protein.identity.info.with.species.category.to.plot.dt <- copy(protein.identity.info.with.species.category.dt)
protein.identity.info.with.species.category.to.plot.dt[, category.corrected.ordered:=factor(category.corrected, levels=protein.identity.info.with.species.category.to.plot.dt[, list(median.identity=median(identity_with_human_ACE2)), category.corrected][order(median.identity), category.corrected])]
protein.identity.info.with.species.category.to.plot.dt[, identity_with_human_ACE2.in.pct:=identity_with_human_ACE2*100]

Fig2B.ggplot <- ggplot(protein.identity.info.with.species.category.to.plot.dt, aes(x=category.corrected.ordered, y=identity_with_human_ACE2.in.pct)) + geom_boxplot() + coord_flip() + theme_pubr(base_size=10) + scale_y_continuous(breaks=seq(0, 100, 20), limits=c(0, 100)) + labs(x="", y="Identity with human ACE2 (%)", fill="") + facet_grid(is.mammals.info~., scales="free_y", space="free_y")

saveRDS(Fig2B.ggplot, paste(sep="", parent.path, "/Fig2B.ggplot.RDS"))
ggsave(paste(sep="", parent.path, "/Fig2B.png"), plot=Fig2B.ggplot, device="png", units="cm", width=21*0.5, height=29.7*0.25)
