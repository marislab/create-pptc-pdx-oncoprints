###create mutational burden matrix
sig.df <- pptc.merge[,c("Tumor_Sample_Barcode","Chromosome","Start_position","Reference_Allele","Tumor_Seq_Allele2")]
#require Sample, chr, pos, ref, alt
names(sig.df) <- c("Sample", "chr", "pos", "ref", "alt")
#### Convert to deconstructSigs input
sigs.input <- mut.to.sigs.input(mut.ref = sig.df, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt")

mut.sum <-apply(sigs.input,1,sum)
log.mut <- log2(mut.sum)
