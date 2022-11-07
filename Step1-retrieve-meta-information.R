library(data.table)
library(reutils)
library(XML)

parent.path <- "results/Step1/NCBI_datasets_download_ortholog_ACE2_human/"

ortholog.info.dt <- fread(paste(sep="", parent.path, "/ncbi_dataset/data/data_table.tsv"))
## the following needs manual annotation (eutils will err)
## 8673  260995 1415580 1825980

tax.info.XML <- efetch(ortholog.info.dt[, tax_id] %>% sort %>% unique %>% {setdiff(., c(8673, 260995, 1415580, 1825980))}, db="taxonomy") %>% content
saveXML(tax.info.XML, paste(sep="", parent.path, "/ortholog.info.xml"))

tax.info.dt <- xpathApply(tax.info.XML, "/TaxaSet/Taxon", fun=function(temp.node){data.table(tax_id=xmlValue(getNodeSet(temp.node, "TaxId")), scientific.name=xpathSApply(temp.node, "LineageEx/Taxon/ScientificName", fun=xmlValue)) }) %>% rbindlist


tax.name.dt <- xpathApply(tax.info.XML, "/TaxaSet/Taxon", fun=function(temp.node){data.table(tax_id=xmlValue(getNodeSet(temp.node, "TaxId")), common.name=xmlValue(getNodeSet(temp.node, "OtherNames/GenbankCommonName")) %>% paste(collapse=";")) }) %>% rbindlist


tax.info.dt[scientific.name=='Primates', category:="primates"]
tax.info.dt[scientific.name=='Rodentia', category:="rodents"]
tax.info.dt[scientific.name=='Carnivora', category:="carnivores"]
tax.info.dt[scientific.name=='Cetacea', category:="whales and dolphins"]
tax.info.dt[scientific.name=='Chiroptera', category:="bats"]
tax.info.dt[scientific.name=='Artiodactyla', category:="even-toed ungulates"]
tax.info.dt[scientific.name %in% c('Diprotodontia', 'Dasyuromorphia', 'Didelphimorphia'), category:="marsupials"]
tax.info.dt[scientific.name=='Mammalia', category:="mammals"]
tax.info.dt[scientific.name=='Aves', category:="birds"]
tax.info.dt[scientific.name=='Amphibia', category:="amphibians"]
tax.info.dt[scientific.name=='Amniota', category:="amniotes"]
tax.info.dt[scientific.name %in% c("Actinopterygii", "Sarcopterygii"), category:="bony fishes and tetrapods"]
tax.info.dt[scientific.name=='Chordata', category:="chordates"]

tax.core.info.dt <- merge(x=tax.info.dt[is.na(category)==FALSE][, temp:=1], y=tax.name.dt, by="tax_id")

fwrite(tax.core.info.dt, paste(sep="", parent.path, "/tax.core.info.dt.csv"))
