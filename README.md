# Mutual Information Between Gut Bacterial Abundance and Crohn’s Disease Status

## Question
Does knowing the abundance of the following gut bacteria reduce uncertainty about whether or not a person has Crohn’s disease?
- Faecalibacterium
- Faecalibacterium prausnitzii 
- Roseburia
- Coprococcus
- Bifidobacterium
- Escherichia coli
- Veillonella
- Bacteroides fragilis
- Akkermansia muciniphila

**Note:** These microbes were chosen because they are well-documented in IBD research

## Hypothesis
Knowing the levels of the bacteria listed above in a patient will significantly reduce uncertainty about whether or not they have Crohn’s disease. 

## Motivation
While plenty of research has been done to fully characterize the microbiome, we want to look at highly studied bacteria species and genera to see if there was a statistical connection between presence and diagnosis. A connection such as this could help improve diagnostic tools clinically.

## Data Information

**Dataset Source**
This project uses data from GMrepo (Project PRJEB42155), a curated multi-omics repository containing metagenomics, 16S rRNA sequencing, metaproteomics (TMT-MS3), and metabolomics profiles from 200 patients.

**Variables Used for Our Analysis From the Raw Data**
- tube_id (unique identifier for each sample taken from a patient)
- Diagnosis (Crohn’s disease, healthy)
- Shotgun metagenomics–derived microbial relative abundance of the 9 selected bacteria 
- Final_Metabolomics_Weight (used for validation steps later)

**Sample Size**
- After filtering the raw data, we have 132 unique tube_ids (18 Healthy, 114 Crohn’s Disease)

## Methods
1. clean the raw data to get the variables we are interested in for our analysis
   - tube_id
   - Diagnosis
   - Shotgun metagenomics–derived microbial relative abundance of the 9 selected bacteria in the sample
   - Final_Metabolomics_Weight (used for validation steps later)
2. Binned relative abundance into 3 bins (zero, low, high)
   Note: computing mutual information requires that continuous data is grouped into bins
3. calculated the mutual information between Diagnosis and each bacteria 
4. ran permutation tests
   Note: shuffeled the disease and healthy labels for each bacteria data 10,000 times and claculated MI each time to obtain a null distribution
5. got p-values
   Note: compared MI value we got to the null distribution in the above step to get the p-value
6. did multiple hypothesis testing correction 
   Note: because we ran 9 hypothesis tests, we need to correct for multiple hypothesis testing, we used false discovery rate of .05 to determine what tests we should keep as significant

## Results
Mutual information values were low across all taxa, indicating substantial overlap between Crohn’s disease and healthy samples, with only a subset remaining significant after FDR correction.

![MI values (bits) for 9 IBD-associated taxa; red dots indicate significant after FDR correction](results/figures/MI_rdotplot.png)











## Requirements
- Python 3.10+
- pandas >= 2.0
- matplotlib >= 3.9
(or: `pip install -r requirements.txt`)

## Repository Structure
- /data/ … 
- /figures/ … 
- /results/ … 
- /scripts/ … 
- README.md … 

## Usage
1. Clone the repo:  
   `git clone https://github.com/<your-username>/<project>.git`
2. Install dependencies:  
   `pip install -r requirements.txt`
3. Open `notebooks/final_analysis.ipynb` or run  
   `python scripts/run_analysis.py`
4. Figures will appear in `/results/figures/`.

## License
MIT

## Acknowledgments
Course: CSDS 313 – Intro to Data Analysis, Instructor Dr. Mehmet Koyuturk
