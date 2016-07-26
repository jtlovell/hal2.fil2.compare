#' @title Gene and mutation plotting
#'
#' @description
#' \code{plotGeneModel} Given a gff file and snpeff annotated snp file,
#' plot the mutations pertaining to a specific gene.
#' @param gff The gff file, must contain the following column names:
#' chr, type, start, end, orientation, info
#' @param snpEffVCF The snpeff annotated vcf file, must contain the following column names:
#' CHR, POS, INFO, where info is the snpeff annotations
#' @param upstreamBuffer How many bp upstream should be plotted
#' @param downstreamBuffer How many bp downstream should be plotted
#' @param features2plot What gff elements should be plotted?
#' @param colors Colors for exon, 3'UTR and 5'UTR respectively
#' @param mutations2annotate Partial string matching for types of differences to annotate
#' @param pchMutation Shape of points that specify annotated sites
#' @param colMutation Color of points that specify annotated sites
#' @param windowSize Sliding window (for difference densities) window size
#' @param stepSize Sliding window step size
#' @param scaleBar If NULL, do not plot a scale bar, otherwise, specify the proportion of
#' total gene space that the scale bar covers
#' @param ... additional arguments, not currently in use.
#'
#' @return a dataframe containing the annotated vcf file of Hal2 Fil2 differences
#' for the gene specified in geneID
#' @examples
#' plotGeneModel(gff = gf, snpEffVCF = sev, geneID = "Pahal.C00786", windowSize=200, stepSize = 5)
#' @export

plotGeneModel<-function(gff, snpEffVCF, geneID, upstreamBuffer = 1000, downstreamBuffer = 500,
                        features2plot = c("exon","five_prime_UTR","three_prime_UTR"),
                        colors = c("steelblue3","lightsteelblue1","lightsteelblue1"),
                        mutations2annotate = c("missense","stop","start"),
                        pchMutation = c(2,8,8), colMutation = c("darkblue","darkred","green"),
                        orientation = NULL, windowSize = 100, stepSize=10,
                        scaleBar=.1, ...){
  #1. Subset gff to the gene of interest
  gff<-gff[grep(geneID,gff$info),]
  if(nrow(gff)==0) stop("geneID not found in the gff file\n")

  #2. Grab orientation information out of gff (unless specified)
  if(is.null(orientation)){
    orientation <- ifelse(gff$orientation[1]=="+","forward","reverse")
  }
  dir <- ifelse(orientation == "forward",2,1)

  #3. subset vcf to gene of interest and features of interest
  vcf<-data.frame(snpEffVCF[grep(geneID, snpEffVCF$INFO),])
  if(nrow(vcf)==0) stop("geneID not found in the vcf file\n")
  gff<-gff[gff$type %in% features2plot,]

  #4. calculate standardized positions relative to the transcription start site
  if(orientation == "forward"){
    xrange<-c(min(c(gff$start,gff$end))-upstreamBuffer,max(c(gff$start,gff$end))+downstreamBuffer)
  }else{
    xrange<-c(min(c(gff$start,gff$end))-downstreamBuffer,max(c(gff$start,gff$end))+upstreamBuffer)

  }

  #5. make basic plot
  par(mar=c(5, 6, 4, 2) + 0.1)
  plot(1,type = "n", xlim=xrange, ylim=c(0,1.2), axes=F, bty="n", xlab = NA, ylab=NA)
  arrows(x0=xrange[1], x1=xrange[2],y0=1,y1=1, length = .1, code = dir)
  features2plot<-features2plot[features2plot %in% gff$type]
  for(j in 1:length(features2plot)){
    gff2<-gff[gff$type == features2plot[j],]
    for(i in 1:nrow(gff2)){
      with(gff2[i,], rect(xleft = min(start,end), xright=max(start, end), ytop = 1.05, ybottom=.95,
                          border = "dodgerblue4" , col = colors[j]))
    }
  }

  #6. some code to make the intron pieces from exon data only
  if(sum(gff$type=="exon")>1){
    ex<-gff[gff$type=="exon",]
    ex<-ex[order(ex$start),]
    ex$ma<-apply(ex[,c("start","end")],1,max)
    ex$mi<-apply(ex[,c("start","end")],1,min)
    ex$sizes<-abs(ex$ma-ex$mi)
    start<-ex$mi[-1]
    end<-ex$ma[-nrow(ex)]
    mids<-(start-end)/2
    mids<-mids+end
    for(i in 1:length(mids)){
      segments(x0=end[i],x1=mids[i],y0=1.05,y1=1.1, col="dodgerblue4", lwd=1)
      segments(x0=mids[i],x1=start[i],y0=1.1,y1=1.05, col="dodgerblue4", lwd=1)
    }
  }


  #7. process annotated VCF file to get the relevant snps
  anns<-lapply(vcf$INFO, function(x) strsplit(x,",", fixed=T)[[1]])
  annsC<-lapply(anns, function(x) {
    y<-x[grep(geneID,x)]
    if(length(y)>1){
      y<-y[!grepl("intergenic_region",y)]
      if(length(y)>1){
        y<-y[1]
      }
    }
    y
  })

  #8. make a column on missense mutation information
  out<-data.frame(effect=unlist(annsC), t(sapply(annsC, function(x) strsplit(x, "|", fixed=T)[[1]][c(2,3,11)])))
  names(out)[2:4]<-c("rel.pos", "effect","AA.change")

  #9. parse amino acid change information
  out$AA.change<-gsub("p.","",out$AA.change, fixed=F)
  out$AA.change<-gsub("\\d", "", out$AA.change)
  out$AA.change<-paste(substr(out$AA.change,1,3),substr(out$AA.change,4,6),sep="-")
  out$AA.change[!grepl("missense",out$rel.pos)]<-NA

  #10. add in other annotation information
  out$plotAnn<-out$AA.change
  out$color<-NA
  out$pch<-NA
  for(i in 1:length(mutations2annotate)){
    mut<-mutations2annotate[i]
    if(!grepl("missense",mut)){
      out$plotAnn[grepl(mut, out$rel.pos)]<-mut
    }
    out$color[grepl(mut, out$rel.pos)]<-colMutation[i]
    out$pch[grepl(mut, out$rel.pos)]<-pchMutation[i]
  }

  #11. combine information and prepare other columns to plot
  tp<-cbind(vcf[,-which(colnames(vcf) == "INFO")], out)
  tp<-tp[order(tp$POS),]
  tp<-tp[tp$POS>=min(xrange) & tp$POS<=max(xrange),]
  #12. add in segments for each SNP
  for(i in 1:nrow(tp)){
    with(tp, segments(x0=POS[i],x1=POS[i],y0=1.01,y1=.99, col="orange", lwd=1))
  }

  #13. Sliding window SNP density
  wind = windowSize
  step = stepSize
  xs<-seq(from = xrange[1], to = xrange[2], by=step)
  sw<-sapply(xs, function(x) sum(abs(tp$POS-x)<=(wind/2)))
  sw<-((sw/max(sw))/5)
  lines(xs,sw+.7)

  #14. Plot points for each type of annotated snp
  toan<-sapply(as.character(tp$rel.pos), function(x)
    any(sapply(mutations2annotate, function(y) grepl(y,x))))
  if(sum(toan)>0){
    mtp<-tp[toan,]
    with(mtp, points(x=POS, y=rep(.5, length(POS)),
                     pch = pch, col = color))
    with(mtp, text(x=POS, y=rep(.3, length(POS)),
                   labels = plotAnn, col = color,
                   srt=90, cex=.6, adj = c(.5,.5)))
  }

  #15. add scale bar
  if(!is.null(scaleBar)){
    seg<-round(diff(xrange)*scaleBar,-2)
    beg<-min(gff$start)
    end<-beg+seg
    segments(x0 = beg, x1 = end, y0=1.11, y1=1.11)
    segments(x0=beg,x1=beg, y0 = 1.1, y1=1.12)
    segments(x0=end,x1=end, y0 = 1.1, y1=1.12)
    text(x = beg, y = 1.11, label = paste(seg,"bp"), adj = c(0,-.5))
  }

  #16. add annotation track labels
  axis(2, at = c(1,.7,.4), labels = c("gene model","mut. density", "mut. annot."), las=2,
       line=-1, lwd=0)


  title(main = paste(geneID, ": ", paste(gff$chr[1], "...", min(gff$start)," - ",max(gff$end)), sep =""))
  return(tp)
}
