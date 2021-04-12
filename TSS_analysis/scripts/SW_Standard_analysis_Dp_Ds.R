#setwd("/scratch/tswenty/GoSTRIPES_sing/STRIPES/tsr_Ds")
setwd("/scratch/scwalls/STRIPEseq_pipelines/GoSTRIPES_sing/STRIPES/figures/4_7_21")

library(TSRchitect)
library(TSRexploreR)
library(GenomicRanges)
library(DESeq2)
library(ChIPseeker)
library(ggseqlogo)
library(ggplot2)

#load("/scratch/tswenty/GoSTRIPES_sing/STRIPES/tsr_Ds/PdSTRIPE_complete.RData")
load("/scratch/scwalls/STRIPEseq_pipelines/GoSTRIPES_sing/STRIPES/tsr_Ds/PdSTRIPE_complete.RData")

#creating the annotation and assembly files
#update both paths below
Dp.annot <- "/data/LynchLabCME/Daphnia/DaphniaTissues/GoSTRIPES_sing_DpTissues/STRIPES/DpGENOME/PA42.4.0.gff"
Dp.assembly <- "/data/LynchLabCME/Daphnia/DaphniaTissues/GoSTRIPES_sing_DpTissues/STRIPES/DpGENOME/PA42.4.1.fasta"

# Generate sample sheet
sample_sheet <- data.frame(
  sample_name=stringr::str_glue("DpDevel_{seq_len(4)}"),
  file_1=NA, file_2=NA,
  condition=rep("TimePoint_D", 4)
)

#writing the tss files to the workspace
tss.1 <- PdSTRIPE@tssCountData[[1]]
tss.2 <- PdSTRIPE@tssCountData[[2]]
tss.3 <- PdSTRIPE@tssCountData[[3]]
tss.4 <- PdSTRIPE@tssCountData[[4]]

colnames(tss.1) <- c("seq","TSS", "strand", "score")
colnames(tss.2) <- c("seq","TSS", "strand", "score")
colnames(tss.3) <- c("seq","TSS", "strand", "score")
colnames(tss.4) <- c("seq","TSS", "strand", "score")

#making granges files from tss data frames
tss.1.gr <- makeGRangesFromDataFrame(tss.1,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.2.gr <- makeGRangesFromDataFrame(tss.2,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.3.gr <- makeGRangesFromDataFrame(tss.3,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.4.gr <- makeGRangesFromDataFrame(tss.4,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

Dp.tss <- list(tss.1.gr, tss.2.gr, tss.3.gr, tss.4.gr)
names(Dp.tss) <- c("D1", "D2", "D3", "D4")

#Creating the TSR explorer object
exp <- tsr_explorer(TSSs=Dp.tss, 
                    genome_annotation=Dp.annot, genome_assembly=Dp.assembly,
                    sample_sheet = sample_sheet
)

#Initial TSS processing

exp <- format_counts(exp, data_type="tss")

#Normalize TSSs
exp <- normalize_counts(exp, data_type = "tss", method = "DESeq2")

## TSS annotation
exp <- annotate_features(exp, data_type = "tss", feature_type="transcript")

##### Sequence analysis
## creating a truncated object for sequence analysis
## some intervals are too close to the edges of short scaffolds, so this was my workaround

plot_threshold_exploration(exp, samples="D1", point_size=1) + scale_color_viridis_c()
ggsave(file="threshold_D1.png")

plot_threshold_exploration(exp, samples="D2", point_size=1) + scale_color_viridis_c()
ggsave(file="threshold_D2.png")

plot_threshold_exploration(exp, samples="D3", point_size=1) + scale_color_viridis_c()
ggsave(file="threshold_D3.png")

plot_threshold_exploration(exp, samples="D4", point_size=1) + scale_color_viridis_c()
ggsave(file="threshold_D4.png")

## Correlation plot: replicates
plot_correlation(
  exp, data_type="tss",
  use_normalized=TRUE, font_size=12,
  heatmap_colors=viridis::viridis(100)
)

ggsave(file="Ds_correlation_matrix.png")

### Genomic distribution analysis:
#### To update with new UTR-less annotation

plot_genomic_distribution(exp, data_type="tss", samples=c("D1", "D2","D3", "D4")) +
  scale_fill_viridis_d(direction=-1, name="Annotation")

ggsave(file="genomic_distribution_Ds.png")

### promoter fraction plot
#### To update with new UTR-less annotation
plot_detected_features(exp, data_type="tss", samples=c("D1", "D2","D3", "D4")) +
  scale_fill_viridis_d(direction=-1)

ggsave(file="promoter_fraction_Ds.png")

### Density plot
plot_density(exp, data_type="tss", samples="D1")
ggsave(file="TSS_density_CDS_D1.png")

plot_density(exp, data_type="tss", samples="D2")
ggsave(file="TSS_density_CDS_D2.png")

plot_density(exp, data_type="tss", samples="D3")
ggsave(file="TSS_density_CDS_D3.png")

plot_density(exp, data_type="tss", samples="D4")
ggsave(file="TSS_density_CDS_D4.png")

### TSS pileup heatmap
## need to update annotation/colour
plot_heatmap(
  exp, data_type="tss", samples="D1",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="TSS_pileup_D1.png")

plot_heatmap(
  exp, data_type="tss", samples="D2",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="TSS_pileup_Ds2.png")

plot_heatmap(
  exp, data_type="tss", samples="D3",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="TSS_pileup_D3.png")

plot_heatmap(
  exp, data_type="tss", samples="D4",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="TSS_pileup_D4.png")

### Sequence logo analysis- all three replicates
plot_sequence_logo(exp, samples="D1")
ggsave(file="D1_seq_logo.png")

plot_sequence_logo(exp, samples="D2")
ggsave(file="D2_seq_logo.png")

plot_sequence_logo(exp, samples="D3")
ggsave(file="D3_seq_logo.png")

plot_sequence_logo(exp, samples="D4")
ggsave(file="D4_seq_logo.png")

### Dinucleotide frequency- all three replicates
#### TODO: repeat As1 with truncated tss data (get code from Taylor)
plot_dinucleotide_frequencies(exp, samples="D1") +
  scale_fill_viridis_c()
ggsave(file="dinucleotide_frequency_D1.png")

plot_dinucleotide_frequencies(exp, samples="D2") +
  scale_fill_viridis_c()
ggsave(file="dinucleotide_frequency_D2.png")

plot_dinucleotide_frequencies(exp, samples="D3") +
  scale_fill_viridis_c()
ggsave(file="dinucleotide_frequency_D3.png")

plot_dinucleotide_frequencies(exp, samples="D4") +
  scale_fill_viridis_c()
ggsave(file="dinucleotide_frequency_D4.png")

### TSS Sequence colour map

plot_sequence_colormap(exp, samples="D1", rasterize=TRUE)
ggsave(file="sequence_colormap_D1.png")
plot_sequence_colormap(exp, samples="D2", rasterize=TRUE)
ggsave(file="sequence_colormap_D2.png")
plot_sequence_colormap(exp, samples="D3", rasterize=TRUE)
ggsave(file="sequence_colormap_D3.png")
plot_sequence_colormap(exp, samples="D4", rasterize=TRUE)
ggsave(file="sequence_colormap_D4.png")

### identify TSRs using clustering
exp <- tss_clustering(exp, threshold=3, n_samples=4, max_distance = 25)

# Associate TSSs with TSRs
exp <- associate_with_tsr(exp)

# Annotate TSRs
exp <- annotate_features(exp, data_type="tsr", upstream=250, downstream=100,
                         feature_type="transcript")

# Mark dominant TSS per TSR
exp <- mark_dominant(exp, data_type="tss")

# Calculate TSR metrics
exp <- tsr_metrics(exp)

###### Add tsr analysis from Standard Analysis documentation

### plot selected tsr metrics
#### TODO: to re-generate using custom script
plot_tsr_metric(exp, tsr_metrics=c("score", "width"), log2_transform=TRUE, samples="all")
ggsave(file="plot_tsr_metrics_Ds.png")

# Dinucleotide motifs by TSR shape
#### TODO: examine As2 and As3, which look odd
plot_sequence_logo(exp, dominant=TRUE, samples="all",
                   data_conditions=conditionals(data_grouping=shape_class)
)
ggsave(file="dinucl_motif_plot_shape_Ds.png")

plot_sequence_logo(exp, dominant=TRUE, samples="D1",
                   data_conditions=conditionals(data_grouping=shape_class)
)

ggsave(file="dinucl_motif_plot_shape_D1.png")

plot_sequence_logo(exp, dominant=TRUE, samples="D2",
                   data_conditions=conditionals(data_grouping=shape_class)
)

ggsave(file="dinucl_motif_plot_shape_D2.png")

plot_sequence_logo(exp, dominant=TRUE, samples="D3",
                   data_conditions=conditionals(data_grouping=shape_class)
)

ggsave(file="dinucl_motif_plot_shape_D3.png")

plot_sequence_logo(exp, dominant=TRUE, samples="D4",
                   data_conditions=conditionals(data_grouping=shape_class)
)

ggsave(file="dinucl_motif_plot_shape_D4.png")

### Plot selected gene track
#### TODO: Taylor check with Bob about rendering more clearly
gene_tracks(
  exp, feature_name="dp_gene7125", use_normalized = TRUE, promoter_only=TRUE,
  samples=c(TSS="D1", TSR="D1")
)

ggsave(file="gene_track_D1_gene7125.png")

