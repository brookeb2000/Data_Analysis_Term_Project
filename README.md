# Term Project – Mutual Information Between Gut Bacterial Abundance and Crohn’s Disease Status

## Note to self:
-work in jupyter lab and go through conda to use the environment I set up for this for juypter lab
-use VS code to open the same project and do the git stuff there
-for this project want to learn how to use branches (Create just one feature branch for your “extra work” or “new analysis.”)

## Overview
Brief 2-3 sentence summary of the project goal and context.
This project analyzes whether knowing the abundance of certain gut bacteria reduces the uncertainty about whether or not a person has Chron's disease. To achieve this, we computed the mutual information ... 

## Data Sources
This dataset from GMrepo contains multiple molecular profiling layers collected from 200 patients, including metagenomics, 16S sequencing, metaproteomics (TMT-MS3), and metabolomics. In this project, we focus specifically on the shotgun metagenomics-derived microbial relative abundance, as well as patient phenotype (Chron disease, colitis, healthy).

## Requirements
- Python 3.10+
- pandas >= 2.0
- matplotlib >= 3.9
(or: `pip install -r requirements.txt`)

## Repository Structure
- /notebooks/ … exploratory & final analysis notebooks
- /scripts/ … helper Python scripts
- /data/ … raw and cleaned data (not tracked if large)
- /results/ … output CSVs and plots
- README.md … project documentation

## Usage
1. Clone the repo:  
   `git clone https://github.com/<your-username>/<project>.git`
2. Install dependencies:  
   `pip install -r requirements.txt`
3. Open `notebooks/final_analysis.ipynb` or run  
   `python scripts/run_analysis.py`
4. Figures will appear in `/results/figures/`.

## Results
Key findings or a sample figure.

## License
MIT

## Acknowledgments
Course: CSDS 313 – Intro to Data Analysis, Instructor X
