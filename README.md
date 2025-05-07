# R Scripts - Biochemistry & Morfology

The analysis is divided into two main groups: **IP (injection)** and **OR (oral)**.
Mouse model - BALB/c ♀

### Data Analysis

The data is analyzed in two experimental groups:
1. **"IP" groups**: VEH_IP, 17A_IP, TUB_IP, COMB_IP, MP_IP - studies related to IP treatments of different therapy and vehicle with just DMSO.
2. **"OR" groups**: VEH_OR, GEN_OR - studies related to oral treatments of genistein and vehicle with just DMSO.

The plots show the distribution of results (the Result) across different groups (e.g., "VEH_IP", "17A_IP", etc.) for different analytes.

The plots also include statistical tests (p-values) performed for each group pair, with significant results marked by an asterisk "*" if p-value < 0.05.

### Requirements

To run the script, you need to have R installed along with the following packages:

```r
install.packages(c("dplyr", "readr", "stringr", "rstatix", "ggplot2", "ggpubr"))

```

### Running analisys

```
git clone https://github.com/your-account/r-scripts-biochemistry.git
```
