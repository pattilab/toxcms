# install.packages
file_path <- "/Users/Lingjue/Documents/Research/PattiLab/Bioinformatics/Rcode/toxcms_1.0.3.tar.gz" #Mac OS
file_path <- "X:/Lingjue/toxcms_1.0.3.tar.gz" #Windows
install.packages(pkgs = file_path, repos = NULL, type="source")

# load packages
library(toxcms)
library(data.table)

# upload example dataset from toxcms
feature <- data.table(read.csv(system.file("extdata", "dataset.csv", package="toxcms")))

# Step 1 statistical analysis
etom_dosestat <- calcdosestat(Feature = feature, Dose_Levels = c("_0uM","_10uM","_50uM","_200uM"), multicomp = "none",
                              p.adjust.method = "none",projectName = "etomoxir_dataset")

# Step 2.1 monotonic trend filtering
etom_drreport_mono <- trendfilter(etom_dosestat, pval_cutoff = 0.05, pval_thres = 1, anova_cutoff = 0.05, trend = "mono",
                                           relChange_cutoff = 0.05, export = FALSE)
# Step 3.1 ED50 modeling and estimation
etom_drreport_fit = fitdrc(DoseResponse_report=etom_drreport_mono, Dose_values=c(0, 10, 50, 200), ED=0.5, export = TRUE,
                           mz_tag = "mzmed", rt_tag = "rtmed", plot=TRUE)

# Step 4.1 clustering
etom_drreport_clust = clusttrend(etom_drreport_mono, reference_index = NULL, sort.method =c("clust","layer"), sort.thres = 20, dist.method = "euclidean", hclust.method = "average",
                     mztag = "mzmed", rttag = "rtmed", heatmap.on = TRUE, plot.all = TRUE, filename = "testdataset_hclust_20clusters.pdf")

# Step 5.1 PCA of all features with color-coded ED50 values
plotpca(DoseResponse_report = etom_drreport_fit, DoseStat = etom_dosestat, EDrange = c(0,200))

# Step 2.2 inflection trend filtering
etom_drreport_reverse <- trendfilter(etom_dosestat, pval_cutoff = 0.05, pval_thres = 1, anova_cutoff = 0.05, trend = "reverse",
                                     relChange_cutoff = 0.05, export = TRUE)

# Step 3.2 inflection trend plotting
plottrend(etom_drreport_reverse, Dose_conditions = c("0uM","10uM","50uM","200uM"), y_transform = T,
          mz_tag = "mzmed", rt_tag = "rtmed")
