comma := ,
empty :=
space := $(empty) $(empty)

include include_methods.mk

## Define methods to include in summary plots
#MTplot := $(MT)
MTplot := edgeRLRT Wilcoxon \
edgeRQLF BPSC DESeq2 DESeq2nofilt edgeRLRTdeconv \
MASTcpm MASTcpmDetRate MASTtpm MASTtpmDetRate SCDE edgeRLRTrobust voomlimma SeuratBimod \
SeuratBimodnofilt SeuratBimodIsExpr2 SeuratTobit DESeq2census edgeRLRTcensus monoclecensus \
monocle limmatrend ROTSvoom ROTScpm ROTStpm metagenomeSeq ttest monoclecount \
edgeRQLFDetRate DESeq2betapFALSE
# D3E 

MTplotc := $(subst $(space),$(comma),$(MTplot))
