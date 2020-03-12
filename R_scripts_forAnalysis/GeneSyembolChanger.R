#################################################
#data import
############################################

raw.data <- read.table("ref_annot.gtf.gene_spans.txt",
                           sep="\t",
                           header=F,
                           stringsAsFactors=FALSE,
                           quote="",
                           comment.char="#")
#raw.data[,6]

for(i in 1:nrow(raw.data)){
  if(raw.data[i, 6] == "."){
    raw.data[i, 6] <- raw.data[i, 1]
  }else{
  }
}

write.table(raw.data,
            file="ref_annot.gtf.gene_spans",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)
