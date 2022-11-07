library("Biostrings")
library("data.table")
library("magrittr")
library("foreach")
library("glue")
library("ggpubr")



readAAStringSet("./all.17.mammals.protein.sequences.alignment.fasta") -> alignment.AAStringSet
sub(pattern="^([^ ]+) .*", replacement="\\1", x=names(alignment.AAStringSet)) -> names(alignment.AAStringSet)
{
    fread("./SARS-CoV-2_CR22.csv") -> virus.info.dt
    virus.info.dt[, Abbreviation] -> .;
    virus.info.dt[, `Protein sequence ID`] -> names(.)
    . -> renamer.vector
}
renamer.vector[names(alignment.AAStringSet)] -> names(alignment.AAStringSet)

alignment.matrix <- as.matrix(alignment.AAStringSet)

## must be of the same length, because it is an alignment
alignment.length <- width(alignment.AAStringSet)[1]
all.positions.vector <- 1:alignment.length
all.positions.on.human.ACE2.vector <- alignment.matrix["Human", ] %>% {temp <- cumsum(. != "-"); temp[. == "-"] <- NA; temp}
alignment.human.ACE2.based.matrix <- alignment.matrix[, which(is.na(all.positions.on.human.ACE2.vector) == FALSE)]

all.label.fixed.dt <- fread("./all.label.fixed.dt.csv")
setkey(all.label.fixed.dt, "variant.name")



foreach(temp.variant.name=sort(unique(all.label.fixed.dt[, variant.name]))) %do% {
    ##
    print(glue("{date()} : processing variant {temp.variant.name}"))

    all.label.fixed.dt[temp.variant.name] -> temp.label.fixed.dt
    ##
    alignment.matrix[temp.label.fixed.dt[label.fixed=="permissive", species.name.ordered], ] -> alignment.permissive.matrix
    ##
    alignment.matrix[temp.label.fixed.dt[label.fixed=="non-permissive", species.name.ordered], ] -> alignment.nonpermissive.matrix
    ##
    if (nrow(alignment.permissive.matrix) <= 6 | nrow(alignment.nonpermissive.matrix) <= 6) {
        return(NULL)
    }
    ##
    set.seed(123)
    alignment.permissive.matrix -> .;
    sample(x=setdiff(rownames(.), "Human"), size=3, replace=FALSE) -> permissive.test.species.vector
    .[rownames(.) %in% permissive.test.species.vector == FALSE, ] -> alignment.permissive.train.matrix
    .[rownames(.) %in% permissive.test.species.vector == TRUE, ] -> alignment.permissive.test.matrix
    alignment.nonpermissive.matrix -> .;
    sample(x=setdiff(rownames(.), "Human"), size=3, replace=FALSE) -> nonpermissive.test.species.vector
    .[rownames(.) %in% nonpermissive.test.species.vector == FALSE, ] -> alignment.nonpermissive.train.matrix
    .[rownames(.) %in% nonpermissive.test.species.vector == TRUE, ] -> alignment.nonpermissive.test.matrix
    print(permissive.test.species.vector)
    print(nonpermissive.test.species.vector)
    rbind(alignment.permissive.test.matrix, alignment.nonpermissive.test.matrix) -> alignment.test.matrix;
    ##
    comparison.dt <- rbindlist(mapply(all.positions.vector, all.positions.on.human.ACE2.vector, FUN=function(temp.position, temp.position.on.human.ACE2){
        return(data.table(
            nonpermissive.sequence.species=rownames(alignment.nonpermissive.train.matrix),
            position=temp.position,
            position.on.human.ACE2=temp.position.on.human.ACE2,
            nonpermissive.sequence.residue=alignment.nonpermissive.train.matrix[, temp.position],
            all.possible.positive.sequence.residues=paste(collapse=",", unique(sort(alignment.permissive.train.matrix[, temp.position]))),
            is.in.positive.sequences=alignment.nonpermissive.train.matrix[, temp.position] %in% alignment.permissive.train.matrix[, temp.position]
        ))
    }, SIMPLIFY=FALSE))
    ##
    nonpermissive.informative.dt <- comparison.dt[is.in.positive.sequences==FALSE & is.na(position.on.human.ACE2) == FALSE]
    negative.informative.residue.positions.on.human.ACE2.vector <- unique(sort(nonpermissive.informative.dt[, position.on.human.ACE2]))
    ##
    rownames(alignment.human.ACE2.based.matrix) -> .;
    expand.grid(
        first.species=.,
        second.species=.,
        stringsAsFactors=FALSE
    ) -> .;
    data.table(.) -> .;
    . -> species.pairs.dt;
    ##
    species.pairs.dt -> .;
    .[
      , weighted.distance := sum((
          alignment.human.ACE2.based.matrix[first.species, negative.informative.residue.positions.on.human.ACE2.vector] !=
          alignment.human.ACE2.based.matrix[second.species, negative.informative.residue.positions.on.human.ACE2.vector]
      )), list(first.species, second.species)] -> .;
    . -> weighted.distances.dt;
    ##
    weighted.distances.dt -> .;
    .[first.species != second.species] -> .;
    merge(x=., y=temp.label.fixed.dt[, list(first.species=species.name, first.species.label=label.fixed)], by="first.species", all.x=TRUE, all.y=FALSE) -> .;
    merge(x=., y=temp.label.fixed.dt[, list(second.species=species.name, second.species.label=label.fixed)], by="second.species", all.x=TRUE, all.y=FALSE) -> .;
    .[, list(R=log10(
                 min(.SD[second.species.label=='non-permissive', weighted.distance]) /
                 min(.SD[second.species.label=='permissive', weighted.distance])
             )), list(first.species, first.species.label)] -> .;
    .[, first.source:=c("train", "test")[(first.species %in% rownames(alignment.test.matrix)) + 1] ] -> .;
    .[, first.species.label.predicted:=c("non-permissive", "permissive")[(R>0)+1] ] -> .;
    ##
    .[, variant.name:=temp.variant.name] -> .;
    . -> temp.final.prediction.dt
} -> final.prediction.dts.collection

rbindlist(final.prediction.dts.collection, use.names=TRUE) -> final.prediction.combined.dt

fwrite(final.prediction.combined.dt, "./final.prediction.combined.dt.csv")

copy(final.prediction.combined.dt) -> .;
.[, list(count=.N), list(variant.name, source=first.source, truth=first.species.label, response=first.species.label.predicted)] -> .;
.[, metric:=c(
        "permissive__permissive"="TP",
        "permissive__non-permissive"="FN",
        "non-permissive__permissive"="FP",
        "non-permissive__non-permissive"="TN"
    )[paste(sep="", truth, "__", response)]] -> .;
dcast(., variant.name + source ~ metric, value.var="count", fill=0) -> .;
.[, `:=`(
    recall=TP/(TP+FN),
    specificity=TN/(TN+FP),
    precision=TP/(TP+FP),
    npv=TN/(TN+FN)
)] -> .;
melt(., id.vars=c("variant.name", "source"), measure.vars=c("recall", "precision", "specificity", "npv"), variable.name="metric", value.name="value") -> .;
.[, source.prettified:=factor(source, levels=c("train", "test"))] -> .;
. -> to.plot.dt

to.plot.dt -> .;
ggplot(., aes(x=variant.name, y=value, fill=metric)) -> .;
. + geom_bar(stat="identity") -> .;
. + geom_hline(yintercept=0.5, linetype="dashed") -> .;
. + coord_flip() -> .;
. + facet_grid(~metric + source.prettified) -> .;
. + theme_pubr(base_size=10) -> .;
. + theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(size=8)) -> .;
. + labs(x="variant") -> .;
. + ggtitle("P:=permissive; N:=non-permissive;\nrecall:=TP/(TP+FN); specificity:=TN/(TN+FP);\nprecision:=TP/(TP+FP); npv:=TN/(TN+FN)") -> .;
ggsave("220929-CR22-evaluation.png", ., device="png", width=21*0.8, height=29*0.8, units="cm")
