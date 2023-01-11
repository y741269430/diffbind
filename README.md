# diffbind

run1  

    vim sort_idx.sh
    #!/bin/bash
    ## sort bam and index for diffbind ##

    cat filenames | while read i; 
    do
    nohup samtools sort -@ 8 ./mapbam/${i}_mm10_bowtie2.mapped.bam -o ./mapbam/${i}_mm10_bowtie2.sorted.bam &&
    samtools index -@ 8 ./mapbam/${i}_mm10_bowtie2.sorted.bam &
    done

run2  

    library(DiffBind)
    library(tidyverse, quietly = TRUE)
    options(stringsAsFactors = F)

    sample <- c('SNI14bFGF_1', 'SNI14bFGF_2', 
                'SNI14_1', 'SNI14_2')

    Treatment <- rep(c('SNI14bFGF', 'SNI14'), each = 2)

    bamReads <- c('featurecount/bam/DRG-CUT-SNI14bFGF-1_mm10_bowtie2.mapped.bam', 
                  'featurecount/bam/DRG-CUT-SNI14bFGF-2_mm10_bowtie2.mapped.bam',
                  'featurecount/bam/DRG-CUT-SNI14-1_mm10_bowtie2.mapped.bam', 
                  'featurecount/bam/DRG-CUT-SNI14-2_mm10_bowtie2.mapped.bam')

    peak <- c('featurecount/SEACR/DRG-CUT-SNI14bFGF-1_seacr_BSNI14.peaks.stringent.bed', 
              'featurecount/SEACR/DRG-CUT-SNI14bFGF-2_seacr_BSNI14.peaks.stringent.bed',
              'featurecount/SEACR/DRG-CUT-SNI14-1_seacr_BSNI14.peaks.stringent.bed', 
              'featurecount/SEACR/DRG-CUT-SNI14-2_seacr_BSNI14.peaks.stringent.bed')

    samples <- data.frame("SampleID" = sample,
                          "Tissue" = NA, 
                          "Factor" = NA,
                          "Condition" = NA,
                          "Treatment" = Treatment,
                          "Replicate" = rep(c(1,2),2), 
                          "bamReads" = bamReads,
                          "ControlID"  = NA,
                          "bamControl" = NA,
                          "Peaks" = peak,
                          "PeakCaller" = rep("bed", 4))

    sampleDba <- dba(sampleSheet = samples)

    sampleCount <- dba.count(sampleDba, summits = 250)

    sampleCon <- dba.contrast(sampleCount, categories = DBA_TREATMENT, minMembers = 2)
    sampleDiff <- dba.analyze(sampleCon, method = DBA_DESEQ2)

    deg <- dba.report(sampleDiff, method=DBA_DESEQ2, contrast = 1, th=1)
    deg2 <- data.frame(deg)

    dba.show(sampleDiff, bContrasts = T)

    peak <- list(dba.report(sampleDiff, contrast = 3),
                 dba.report(sampleDiff, contrast = 1),
                 dba.report(sampleDiff, contrast = 2))

    names(peak) <- c('CTRL_ION7', 'ION7TP_ION7', 'ION7TP_CTRL')
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    peakAnnoList <- lapply(peak, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 0), 
                           annoDb="org.Mm.eg.db", verbose=FALSE, overlap="all")
    peakAnno_df <- lapply(peakAnnoList, function(x){x <- as.data.frame(x)})

    save(sampleCount, sampleCon, sampleDiff, peakAnnoList, peakAnno_df, file='count/all_diffbind.RData')
