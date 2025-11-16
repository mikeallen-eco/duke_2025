# Duke Farms Grassland Bird Survey 2025

This repository contains code and data to reproduce the analysis in the report:

**_“Long-term grassland bird monitoring on Duke Farms: 2010–2025”_**  
Prepared by **Michael C. Allen, Ph.D., Tsuga Biodiversity Insights, LLC**

---

## Repository Structure & Usage

This project is organized as an RStudio-ready workflow using the project file `duke_2025.Rproj`.  
Running the analysis typically begins with `main.R`, which:

- Sources `setup.R` to install required packages and functions  
- Calls the seven analysis steps, each implemented as a function in the `R/` directory  

Additional components:

- **`bugs/`** contains the hierarchical distance sampling (HDS) JAGS model  
- **`data/`** contains raw and processed bird survey data from 2010–2025  
    - The 2025 dataset is stored in `20250520_grassland_birds_2025_Allen/`  
    - Scanned datasheets for the 2025 field season are in `scanned_datasheets/`

---

## Directory Structure

```text
duke_2025/
├── duke_2025.Rproj
├── main.R
├── README.md
│
├── bugs/
│   └── distance_model_2025.txt
│
├── R/
│   ├── setup.R
│   ├── step1A_format_2010_2013_data_for_analysis.R
│   ├── step1B_format_2018_2019_data_for_analysis.R
│   ├── step1C_format_2024_2025_data_for_analysis.R
│   ├── step1D_get_jags_data_for_distance_sampling.R
│   ├── step2_run_JAGS_mod.R
│   ├── step3_plot_HDS_density.R
│   ├── step4_get_HDS_density_function.R
│   ├── step5_jags_table_html_function.R
│   ├── step6_plot_raw_density_function.R
│   └── step7_map_mean_counts_2025_function.R
│
├── data/
│   ├── Grassland_Birds_2010-2013_data.csv
│   ├── Grassland_Birds_2018-2019_data.csv
│   ├── Grassland_Birds_2024_data.csv
│   ├── Kaufman_Skeet.geojson
│   ├── pt_coords_2010_2025.csv
│   │
│   ├── 20250520_grassland_birds_2025_Allen/
│   │   ├── Grassland_Birds_2025_data.csv
│   │   ├── Grassland_Birds_2025_methods.docx
│   │   └── Grassland_Birds_2025_ReadMe.txt
│   │
│   └── scanned_datasheets/
│       ├── Kaufman_05202025_CTBVO.pdf
│       ├── Kaufman_06062025_MA.pdf
│       ├── Kaufman_06232025_CTBVO.pdf
│       ├── Kaufman_07072025_MA.pdf
│       ├── Skeet_05202025_MA.pdf
│       ├── Skeet_06042025_CTBVO.pdf
│       ├── Skeet_06242025_MA.pdf
│       └── Skeet_07072025_CTBVO.pdf
```

## Package Versions
Analysis was conducted using:
R 4.5.1 (2025-06-13)
JAGS 4.3.1
R package versions:
```text
dplyr_1.1.4
tidyr_1.3.1
lubridate_1.9.4
ggplot2_4.0.0
jagsUI_1.6.2
sf_1.0-21
scales_1.4.0
kableExtra_1.4.0
knitr_1.50
patchwork_1.3.2
ggspatial_1.1.10
```