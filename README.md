# Data and code for Rosati, Pogorelyy, Minervina et al 2021.

This repository contains R scripts reproducing the analysis from the manuscript.

0) Download the repository.

1) To reproduce the pipeline from the merged processed single cell dataset, see _TCR_analysis.R_. 

This analysis depends on _tcrdist_ implementation from [conga](https://github.com/phbradley/conga) package.
Download conga repository, and compile _tcrdist_:
```
cd conga/tcrdist_cpp
make
```

then edit the third line in TCR_analysis.R to include path to your _tcrdist_ folder: 
```
path_tcrdist="path_to_conga_in_your_system/tcrdist_cpp/"
```

2) To reproduce all the steps to generate merged processed single cell dataset, see the scripts in _preprocessing_ folder (all done by [@erosix](https://github.com/erosix)).