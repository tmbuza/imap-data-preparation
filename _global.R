#############
# From CRAN using install.packages()
#############
if(!require("BiocManager")){install.packages("BiocManager")}
library("BiocManager")

if(!require("remotes")){install.packages("remotes")}
library("remotes")

if(!require("devtools")){install.packages("devtools")}
library("devtools")

if(!require("tidyverse")) {install.packages("tidyverse")}
suppressPackageStartupMessages(library("tidyverse"))

if(!require("schtools")){install.packages("schtools")}
library(schtools)

if(!require("knitr")){install.packages("knitr")}
library(knitr)

if(!require("rmarkdown")){install.packages("rmarkdown")}
suppressPackageStartupMessages(library(rmarkdown))

if(!require("servr")){install.packages("servr")}
library(servr)

if(!require("bookdown")){install.packages("bookdown")}
library(bookdown)

if(!require("tools")){install.packages("tools")}
library(tools)

if(!require("yaml")){install.packages("yaml")}
library(yaml)

if(!require("Matrix")){install.packages("Matrix")}
library(Matrix)

if(!require("vegan")){install.packages("vegan")}
library(vegan)

if(!require("svglite")){install.packages("svglite")}
library(svglite)

if(!require("igraph")) {install.packages("igraph")}
library(igraph)

if(!require("ggtext")) {install.packages("ggtext")}
library(ggtext)

if(!require("schtools")) {install.packages("schtools")}
library(schtools)

if(!require("ggraph")) {install.packages("ggraph")}
library(ggraph)

if(!require("plotly")) {install.packages("plotly")}
library(plotly)

if(!require("kableExtra")) {install.packages("kableExtra")}
library(kableExtra)

if(!require("ggdendro")) {install.packages("ggdendro")}
library(ggdendro)

if(!require("dendextend")) {install.packages("dendextend")}
library(dendextend)

if(!require("ggimage")) {install.packages("ggimage")}
library(ggimage)

if(!require("ggimage")) {install.packages("ggimage")}
library(ggimage)

if(!require("TDbook")) {install.packages("TDbook")}
library(TDbook)

if(!require("ggnewscale")) {install.packages("ggnewscale")}
library(ggnewscale)

if(!require("rsvg")) {install.packages("rsvg")}
library(rsvg)

if(!require("emojifont")) {install.packages("emojifont")}
library(emojifont)

if(!require("emoji")) {install.packages("emoji")}
library(emoji)

if(!require("GGally")) {install.packages("GGally")}
library(GGally)

if(!require("reticulate")) {install.packages("reticulate")}
library(reticulate)



#############
# From GitHub using devtools or remotes
#############

if(!require("funModeling")){devtools::install_github("pablo14/funModeling")}
suppressPackageStartupMessages(library("funModeling"))

if(!require("qiime2R")){remotes::install_github("jbisanz/qiime2R", dependencies = TRUE)}
library(qiime2R)

if(!require("metagMisc")){remotes::install_github("vmikk/metagMisc", dependencies = TRUE)}
library(metagMisc)

if(!require("microViz")){devtools::install_github("david-barnett/microViz")}
library(microViz)

if(!require("gifski")){devtools::install_github("r-rust/gifski")}
library(gifski)

if(!require("MicEco")){remotes::install_github("Russel88/MicEco")}
library(MicEco)




##############
# From bioconda using BiocManager
##############

if(!require("phyloseq")){BiocManager::install("phyloseq")}
library(phyloseq)

if(!require("microbiome")){BiocManager::install("microbiome")}
library(microbiome)

if(!require("rhdf5")) {BiocManager::install("rhdf5")}
library(rhdf5)

if(!require("ggtree")) {BiocManager::install("ggtree")}
library(ggtree)

if(!require("lefser")) {BiocManager::install("lefser")}
library(lefser)

if(!require("DESeq2")) {BiocManager::install("DESeq2")}
library(DESeq2)

if(!require("DESeq2")) {BiocManager::install("DESeq2")}
library(DESeq2)

if(!require("dada2")) {BiocManager::install("dada2")}
library(dada2)

if(!require("microbiomeMarker")) {BiocManager::install("microbiomeMarker")}
library(microbiomeMarker)



plotbeta<-function(physeq,group,shape=NULL,distance="bray",method="PCoA",color=NULL,size=3,ellipse=FALSE){
    if(!taxa_are_rows(physeq)){
        physeq <- t(physeq)
    }
    beta<-betadiv(physeq,distance = distance,method=method)
    df <- as.data.frame(beta$beta)
    PCs <- beta$PCs
    tab <- as(sample_data(physeq),"data.frame")
    df <- cbind(df[,1:4],tab[rownames(df),])
    df$group<-tab[,group]
    if(is.null(color)){
        color<-distcolor[1:length(unique(df$group))]
    }
    if(!is.null(shape)){
        df$shape<-tab[,shape]
        p <- ggplot(df,aes_string("Axis.1","Axis.2",color="group",shape="shape"))
    }else{
        p <- ggplot(df,aes_string("Axis.1","Axis.2",color="group"))
    }
    p<-p+geom_point(size=size)+scale_color_manual(values=color)
    p <- p+theme_light(base_size=15)+xlab(paste0("Axis1 (",round(PCs[1]*100,2),"%)"))+ylab(paste0("Axis2 (",round(PCs[2]*100,2),"%)"))
    if(isTRUE(ellipse)){
        p <- p + stat_ellipse()
    }
    p
}

plotalpha<-function(physeq,group,method=c("Observed","Simpson", "Shannon"),color=NULL,geom="boxplot",
                    pvalue=0.05,padj=NULL,sig.only=TRUE, wilcox=FALSE,show.number=FALSE){
    if (!taxa_are_rows(physeq)) {
        physeq <- t(physeq)
    }
    rich<-richness(physeq,method = method)
    name<-levels(factor(colnames(rich)))
    tab<-as(sample_data(physeq),"data.frame")
    rich$group<-tab[rownames(rich),group]
    if(isTRUE(wilcox)){
        res<-do_wilcox(rich,"group")
    }else{
        res<-do_ttest(rich,"group")
    }
    if(sum(res$p<pvalue)<1){
        message("No significant difference between any of the groups")
        pvalue = 1
    }
    if(isTRUE(sig.only)){
        res <- subset(res,p < pvalue)
        if(!is.null(padj)){
            res <- res[res$p.adj < padj,]
        }
    }
    vals<-rich%>%gather(type,val,-group)%>%group_by(type,group)%>%summarise(ma=max(val))%>%spread(group,ma)
    pos <- apply(res, 1, function(x)max(vals[vals$type==x[1],x[3:4]]))
    mpos <- apply(res, 1, function(x)min(vals[vals$type==x[1],x[3:4]]))
    if(geom=="boxplot"){
        p<-rich%>%gather(type,val,-group)%>%ggboxplot(x="group",y="val",color="group")
    }else if(geom=="violin"){
        p<-rich%>%gather(type,val,-group)%>%ggviolin(x="group",y="val",color="group")
    }else if(geom=="dotplot"){
        p<-rich%>%gather(type,val,-group)%>%ggdotplot(x="group",y="val",color="group")
    }else{
        stop("Please specify one type of boxplot,violin,dotplot")
    }
    if(!isTRUE(show.number)){
        res$p.signif<-sapply(res$p,function(x).getstar(x))
        res$p.adj.signif<-sapply(res$p.adj,function(x).getstar(x))
    }else{
        res$p.signif<-res$p
        res$p.adj.signif<-res$p.adj
    }
    if(is.null(color)){
        color<-distcolor[1:length(unique(rich$group))]
    }
    p<-facet(p,facet.by = "type",scales = "free_y",ncol = length(method))
    if(!is.null(padj)){
        p<-p+stat_pvalue_manual(res,label = "p.adj.signif",y.position = pos+2*mpos/nrow(res))
    }else{
        p<-p+stat_pvalue_manual(res,label = "p.signif",y.position = pos+mpos/nrow(res))
    }
        p<-p+xlab("")+ylab("")+
        theme(legend.position = "none",axis.text.x=element_text(angle=90,vjust=0.5, hjust=1))+
        scale_color_manual(values=color)
    p
}

plotbar<-function(physeq,level="Phylum",color=NULL,group=NULL,top=5,return=FALSE,fontsize.x = 5, fontsize.y = 12){
    pm <- psmelt(physeq)
    if(is.null(color)){
        len<-length(unique(pm[,level]))
        color<-distcolor[1:len]
    }
    if(is.null(group)){
        group_var<-c("Sample",level)
    }else{
        group_var<-c(group,level)
    }
    d<-pm%>%group_by_at(vars(one_of(group_var)))%>%summarise(su=sum(Abundance))
    d <- as.data.frame(d)
    d[,level][is.na(d[,level])]<-"NA"
    dx <- pm%>%group_by_at(vars(one_of(level)))%>%summarise(su=sum(Abundance))
    dx <- dx[order(dx$su,decreasing = T),]
    sel <- dx%>%head(top)%>%select(!!level)%>%pull(1)
    d <- d[d[,level]%in%sel,]
    if(is.null(group)){
        p<-ggplot(d,aes_string("Sample","su",fill=level))
    }else{
        p<-ggplot(d,aes_string(group,"su",fill=level))
    }
    p<-p+geom_bar(stat = "identity",position = "fill")+scale_fill_manual(values=color)+
    theme_light()+
    scale_y_continuous(expand = c(0, 0.001)) +
    theme(axis.text.x=element_text(angle=90,size=fontsize.x, vjust=0.5, hjust=1),
              axis.text.y=element_text(size=fontsize.y),
              panel.background = element_blank(),axis.ticks.x = element_blank())+
        xlab("")+ylab("")
    if(isTRUE(return)){
        return(pm[,c("OTU","Abundance",group_var)])
    }else{
        return(p)
    }
}

plotdiff<-function(res,level="Genus",color=NULL,pvalue=0.05,padj=NULL,log2FC=0,size=3,fontsize.x=5,fontsize.y=10,horiz=TRUE){
    if(!is.null(padj)){
        pval<-padj
        sigtab <- subset(res,padj<pval&abs(log2FoldChange)>log2FC)
    }else{
        pval<-pvalue
        sigtab <- subset(res,pvalue<pval&abs(log2FoldChange)>log2FC)
    }
    x <- tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
    x <- sort(x, TRUE)
    sigtab$Phylum <- factor(as.character(sigtab$Phylum), levels=names(x))
    if(is.null(color)){
        len<-length(unique(sigtab$Phylum))
        color<-distcolor[1:len]
    }
    # Genus order
    sigtab$name<-paste0(sigtab[,level],"(",rownames(sigtab),")")
    x <- tapply(sigtab$log2FoldChange, sigtab$name, function(x) max(x))
    x <- sort(x, TRUE)
    sigtab$name <- factor(as.character(sigtab$name), levels=names(x))
    p <- ggplot(sigtab, aes_string(x="name", y="log2FoldChange", color="Phylum"))+
        geom_point(size=3) +theme_light()+xlab(level)+
        theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=fontsize.x),
              axis.text.y = element_text(size=fontsize.y))+
        scale_color_manual(values=color)
    if(isTRUE(horiz)){
        p<-p+coord_flip()+theme(axis.text.x=element_text(angle=0,size=fontsize.x))
    }
    p
}

plotLDA<-function(x,group,lda=2,pvalue=0.05,padj=NULL,color=NULL,fontsize.x=4,fontsize.y=5){
    x <- subset(x,LDAscore>lda)
    if(!is.null(padj)){
        x <- subset(x,p.adj<padj)
    }else{
        x <- subset(x,p.value<pvalue)
    }
    x <- subset(x,direction%in%group)
    x<-x %>%mutate(LDA=ifelse(direction==group[1],LDAscore,-LDAscore))
    p<-ggplot(x,aes(x=reorder(tax,LDA),y=LDA,fill=direction))+
        geom_bar(stat="identity",color="white")+coord_flip()+
        theme_light()+theme(axis.text.x = element_text(size=fontsize.x),
                            axis.text.y = element_text(size=fontsize.y))
    if(is.null(color)){
        color <- distcolor[c(2:3)]
    }
    p<-p+scale_fill_manual(values=color)+xlab("")
    p
}


plotmarker<-function(x,level="Genus",top=30,rotate=FALSE,dot.size=8,label.color="black",label.size=6){
    x <- x[1:top,]
    x <- x[order(x$Value),]
    x$label<-paste0(x[,level],"(",x$OTU,")")
    p<-ggdotchart(x,x="label",y="Value",add="segments",color=I("#00AFBB"),rotate=rotate,dot.size=dot.size,sorting="descending",
                  add.params = list(color = "#00AFBB", size = 1.5),
               label=round(x$Value,2),font.label = list(color = label.color, size = label.size,vjust=0.2))
    if(isTRUE(rotate)){
        p<-p+xlab(level)+ylab("Mean Decrease Accuracy")
    }else{
        p<-p+xlab(level)+ylab("Mean Decrease Accuracy")
    }
    p
}


plotquality<-function(file,n = 5e+05, aggregate = FALSE){
    dada2::plotQualityProfile(file,n=n,aggregate = aggregate)
}



##############
 # Installing specific version use:
##############
# devtools::install_version("package", version = "x.x.x")
# remotes::install_version("package", version = "x.x.x"