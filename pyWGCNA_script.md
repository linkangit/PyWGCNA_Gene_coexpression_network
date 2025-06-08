```python
#pip install pywgcna
```


```python
# python
###############################################################################
# PyWGCNA analysis: 16-sample 2×2 design (Genotype × Treatment)
# Author: Linkan Dash      Date: 6/2/25
###############################################################################

import os
import pandas as pd
from PyWGCNA.geneExp import GeneExp                 # expression container
from PyWGCNA.wgcna     import WGCNA                 # main class
from PyWGCNA.utils     import getGeneList           # optional gene-info helper
```


```python
###############################################################################
# 0. User paths
###############################################################################
expr_file        = "expression_TPM.csv"           # genes rows, samples cols
geneinfo_file    = "gene_info.csv"                # optional; can be None
sampleinfo_file  = "sample_info.csv"              # required

species          = "arabidopsis"                       ### <— CHANGED
level            = "gene"
output_dir       = "WGCNA_results"
figtype          = "pdf"

os.makedirs(output_dir, exist_ok=True)
```


```python
###############################################################################
# 1. Load data
###############################################################################
expr  = pd.read_csv(expr_file, index_col=0)

# Gene annotation ------------------------------------------------------------
if geneinfo_file and os.path.exists(geneinfo_file):
    ginfo = pd.read_csv(geneinfo_file, index_col=0)
else:
    # BioMart dataset for A. thaliana resides on Ensembl Plants
    ginfo = getGeneList(dataset = "athaliana_eg_gene",   ### <— CHANGED
                        server_domain = "http://plants.ensembl.org/biomart")
    # keep only rows present in the expression set
    ginfo = ginfo.set_index("gene_id").loc[expr.index].reset_index()

# Sample annotation ----------------------------------------------------------
sinfo = pd.read_csv(sampleinfo_file, index_col=0)
sinfo = sinfo.loc[expr.columns]                        # enforce same order

```


```python
###############################################################################
# 2. Wrap in GeneExp
###############################################################################
geneExp = GeneExp(species     = species,
                  level       = level,
                  geneExp     = expr.T,
                  geneInfo    = ginfo,
                  sampleInfo  = sinfo)

```


```python
###############################################################################
# 3. Create WGCNA object  (parameters tuned for n = 16)
###############################################################################
wgcna = WGCNA(name           = "Ath_32s",
              species        = species,
              level          = level,
              geneExp       = expr.T,
              TPMcutoff      = 1,
              powers         = list(range(1,21)),
              RsquaredCut    = 0.85,
              MeanCut        = 100,
              networkType    = "signed hybrid",
              TOMType        = "signed",
              minModuleSize  = 40,
              MEDissThres    = 0.25,
              save           = True,
              outputPath     = output_dir,
              figureType     = figtype)
```

    [92mSaving data to be True, checking requirements ...[0m



```python
###############################################################################
# 4. Pipeline
###############################################################################
print("STEP 1/4  Pre-processing …")
wgcna.preprocess(show=False)

print("STEP 2/4  Detecting modules …")
wgcna.findModules(kwargs_function={
        'cutreeHybrid': {'deepSplit': 2, 'pamRespectsDendro': False}
})
```

    STEP 1/4  Pre-processing …
    [1m[94mPre-processing...[0m
    	Detecting genes and samples with too many missing values...
    	Done pre-processing..
    
    STEP 2/4  Detecting modules …
    [1m[94mRun WGCNA...[0m
    [96mpickSoftThreshold: calculating connectivity for given powers...[0m
    will use block size  1814
        Power  SFT.R.sq     slope truncated R.sq      mean(k)    median(k)  \
    0       1  0.432981  0.474847        0.31164  7197.583043  7757.621859   
    1       2  0.044119 -0.159918      -0.006917  3928.711308  4116.915813   
    2       3  0.259218 -0.416978       0.340681  2463.539052  2577.862693   
    3       4  0.403141 -0.641576       0.553558  1673.235057  1692.538843   
    4       5  0.484318 -0.801255       0.676557  1198.261681  1144.967467   
    5       6  0.534536 -0.934762       0.742483   891.572452   791.543942   
    6       7  0.594596 -1.002508       0.805697   683.067526    561.56052   
    7       8   0.63665 -1.062506       0.846323   535.654861   407.735041   
    8       9  0.668406 -1.122521       0.871724   428.159751   298.332003   
    9      10   0.66732 -1.213331       0.875617   347.772936   222.074016   
    10     11  0.672484 -1.286049       0.882409   286.382123   166.997224   
    11     12  0.682908 -1.336928       0.889499   238.654447   126.567558   
    12     13  0.692914 -1.385179       0.894252   200.973756    99.289999   
    13     14  0.690917 -1.425538       0.891583    170.82273    78.301028   
    14     15  0.694045 -1.464295       0.888268   146.408984    62.545516   
    15     16  0.693265 -1.494619       0.882649   126.430994    51.073805   
    16     17  0.698184   -1.5028       0.880603   109.927199    43.367549   
    17     18  0.694703 -1.518656       0.868248     96.17623     36.35825   
    18     19  0.721902 -1.498754       0.876988    84.629438    31.077066   
    19     20  0.767924 -1.465155       0.902988    74.864352    26.584156   
    
              max(k)  
    0   12091.170441  
    1    8061.273104  
    2    5904.210163  
    3    4618.164292  
    4    3746.915382  
    5    3098.130177  
    6    2598.340082  
    7     2203.69689  
    8    1897.189291  
    9    1668.129741  
    10   1476.940221  
    11    1315.50128  
    12    1177.83617  
    13   1059.971397  
    14    958.090564  
    15    869.222271  
    16    791.247789  
    17    722.469935  
    18    661.514155  
    19    607.256222  
    [92mNo power detected to have scale free network!
    Found the best given power which is 20.[0m
    [96mcalculating adjacency matrix ...[0m
    	Done..
    
    [96mcalculating TOM similarity matrix ...[0m
    	Done..
    
    [96mGoing through the merge tree...[0m
    ..cutHeight not given, setting it to 0.993  ===>  99% of the (truncated) height range in dendro.
    	Done..
    
    [96mCalculating 32 module eigengenes in given set...[0m
    	Done..
    
    mergeCloseModules: Merging modules whose distance is less than 0.25
    fixDataStructure: data is not a Dictionary: converting it into one.
    multiSetMEs: Calculating module MEs.
      Working on set 1 ...
    [96mCalculating 32 module eigengenes in given set...[0m
    	Done..
    
    multiSetMEs: Calculating module MEs.
      Working on set 1 ...
    [96mCalculating 29 module eigengenes in given set...[0m
    	Done..
    
    multiSetMEs: Calculating module MEs.
      Working on set 1 ...
    [96mCalculating 27 module eigengenes in given set...[0m
    	Done..
    
    multiSetMEs: Calculating module MEs.
      Working on set 1 ...
    [96mCalculating 26 module eigengenes in given set...[0m
    	Done..
    
      Calculating new MEs...
    multiSetMEs: Calculating module MEs.
      Working on set 1 ...
    [96mCalculating 26 module eigengenes in given set...[0m
    	Done..
    
    [96mCalculating 26 module eigengenes in given set...[0m
    	Done..
    
    fixDataStructure: data is not a Dictionary: converting it into one.
    orderMEs: order not given, calculating using given set 0
    	Done running WGCNA..
    



    
![png](output_6_1.png)
    



    
![png](output_6_2.png)
    



```python
# python
import os
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as sch
import networkx as nx                         # only for the hub-gene graph

FIGDIR = wgcna.outputPath                     # convenience variable
```


```python
plt.figure(figsize=(5,4))
sft = wgcna.sft                                # dataframe created by pickSoftThreshold
plt.scatter(sft.Power, sft['SFT.R.sq'], c='k')
plt.axhline(0.85, ls='--', c='red')
plt.xlabel('Soft–threshold power')
plt.ylabel('Scale-free fit $R^{2}$')
plt.title('Scale-free topology fit')
plt.tight_layout()
plt.savefig(os.path.join(FIGDIR, 'soft_threshold_fit.pdf'))
plt.show()
```


    
![png](output_8_0.png)
    



```python
# after you have created 'wgcna' and already loaded sinfo
wgcna.updateSampleInfo(sampleInfo = sinfo)
```


```python
# 1) Attach categorical covariates and interactions
wgcna.datExpr.obs["Genotype"]      = sinfo["Genotype"].values          # WT / Mut
wgcna.datExpr.obs["Treatment"]     = sinfo["Treatment"].values         # Unt / Trt
wgcna.datExpr.obs["Time"]          = sinfo["Time"].values              # T0 / T1 / T2 / ...
wgcna.datExpr.obs["GenTxTrt"]      = (
    sinfo["Genotype"] + "_" + sinfo["Treatment"]).values               # WT_Unt etc.
wgcna.datExpr.obs["TrtxTime"]      = (
    sinfo["Treatment"] + "_" + sinfo["Time"]).values 
wgcna.datExpr.obs["GenTxTime"]      = (
    sinfo["Genotype"] + "_" + sinfo["Time"]).values 
wgcna.datExpr.obs["GenTxTrtxTime"]  = (
    sinfo["Genotype"] + "_" + sinfo["Treatment"] + "_" + sinfo["Time"]).values  # WT_Unt_T0 etc.

# 2) Convert to numeric codes (so Pearson r is meaningful)
for col in ["Genotype", "Treatment", "GenTxTrt", "Time", "TrtxTime", "GenTxTime", "GenTxTrtxTime"]:
    wgcna.datExpr.obs[col] = (
        wgcna.datExpr.obs[col].astype("category").cat.codes
    )
```


```python
print(wgcna.datExpr.obs.head())
# Should now list Genotype, Treatment, GenTxTrt
```

                   Genotype  Treatment  Time Batch  GenTxTrt  TrtxTime  GenTxTime  \
    WT_Unt_5d_B1          1          1     1    B1         3         3          3   
    WT_Unt_5d_B2          1          1     1    B2         3         3          3   
    WT_Unt_5d_B3          1          1     1    B3         3         3          3   
    WT_Unt_5d_B4          1          1     1    B4         3         3          3   
    Mut_Unt_5d_B1         0          1     1    B1         1         3          1   
    
                   GenTxTrtxTime  
    WT_Unt_5d_B1               7  
    WT_Unt_5d_B2               7  
    WT_Unt_5d_B3               7  
    WT_Unt_5d_B4               7  
    Mut_Unt_5d_B1              3  



```python
meta_cols = ["Genotype", "Treatment", "GenTxTrt", "Time", "TrtxTime", "GenTxTime", "GenTxTrtxTime"]
wgcna.module_trait_relationships_heatmap(
        meta_cols,
        alternative = "two-sided", figsize=(50, 10))
```


    
![png](output_12_0.png)
    



```python
# python
###############################################################################
# 5.  EXTRA VISUALISATIONS  +  EXPORT OF CORE RESULTS
###############################################################################
import os, itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy.cluster.hierarchy import dendrogram, linkage


# --------------------------------------------------------------------------- #
#  5A-3  Heat-map of module eigengenes across samples                         #
# --------------------------------------------------------------------------- #
if hasattr(wgcna, "MEs"):
    mes = wgcna.MEs.copy()
else:
    mes = wgcna.moduleEigengenes

# order samples by any trait you like – here by Genotype/Treatment
sample_order = sinfo.sort_values(["Genotype", "Treatment"]).index
mes = mes.loc[sample_order]

plt.figure(figsize=(21, 12))
sns.heatmap(mes.T, cmap="vlag", center=0, cbar_kws=dict(label="Eigengene\nexpression (z)"))
plt.yticks(rotation=0)
plt.title("Module eigengenes across samples")
plt.xlabel("Samples")
plt.tight_layout()
plt.savefig(os.path.join(FIGDIR, "eigengene_heatmap.pdf"))
plt.show()
```


    
![png](output_13_0.png)
    



```python
import pandas as pd
from pathlib import Path

out_file = Path("WGCNA_results/gene_module_assignment.csv")

# 1.  The colour (module) that each gene belongs to
# -------------------------------------------------
# • wgcna.datExpr.columns  → gene IDs (because datExpr = expr.T)
# • wgcna.moduleColors     → parallel list/array of colour labels
gene2mod = pd.DataFrame(pd.Series(wgcna.datExpr.var.moduleColors,
                      index=wgcna.datExpr.var.index,
                      name="Module"))

# (Optional) add any extra gene annotation you have already loaded
# ----------------------------------------------------------------
gene2mod = gene2mod.merge(ginfo, left_on="gene_id", right_on="gene_id", how="left")

gene2mod.to_csv(out_file, index=False)
print(f"Gene–module table written to → {out_file.resolve()}")
```

    Gene–module table written to → /home/jovyan/pyWGCNA_new/WGCNA_results/gene_module_assignment.csv



```python

```
