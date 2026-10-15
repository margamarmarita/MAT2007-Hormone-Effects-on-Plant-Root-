# MAT2007-Final-Project---Hormone-Effects-on-Plant-Root-Lengths
Final project using C++ for statistical analysis and Python for visualizations, answering the **research question**: "How strongly does each hormone change total root length compared to control plants, and do the three hormones differ from each other?"

**Project Description**

This project investigates the effects of different phytohormones, or plant-specific hormones, on the total root length in different samples of _Arabidopsis thaliana_ (model organism for plant genetic studies).

The datasets contain mean root lengths of each treatment: _control_ (untreated), _auxin_ (main regulator of growth and development), _cytokinin_ (cell-division regulator), and _abscisic acid_ (stress-response hormone). The average root lengths of each condition were calculated, alongside standard errors, with the intent of determining the one with the greatest effect on root growth. Additionally, correlation was calculated to determine any significant relationship between the hormones.

**Dataset** 

The dataset comes from AraPheno, a website contianing data on different phenotypic values (expressed traits) of  _Arabidopsis thaliana_. 

- **Links**: https://arapheno.1001genomes.org/study/27/ (Control), https://arapheno.1001genomes.org/study/28/ (Auxin), https://arapheno.1001genomes.org/study/30/ (Cytokinin), and https://arapheno.1001genomes.org/study/29/ (Abscisic acid)
- **Original Article:** Ristova, D., Giovannetti, M., Metesch, K., & Busch, W. (2018). Natural genetic variation shapes root system responses to phytohormones in Arabidopsis. The Plant Journal, 96(2), 468–481. https://doi.org/10.1111/tpj.14034

The columns utilized in the research paper include _accession_id_ (ID of each accession) and _phenotype_value_ (mean total root length of each accession). The attached C++ script downloads the needed .CSV files automatically. 

**Requirements**

This code requires a C++17 compiler (e.g. g++, clang++, or MSVC), curl, and Python 3. The Python requires the packages numpy, pandas, and matplotlib. 

**Running the Code**

1. Download the C++ and Python script. Make sure they are in the same folder.
2. Open the C++ code in VS Code. 
3. Make sure to have a working compiler and curl. When sure of this, compile the C++ code.
4. Run the program.
5. Check that the .CSV files are downloaded: results_summary.csv, pairwise_comparisons.csv, correlations.csv, ratios_per_accession.csv, and ratio_histogram.csv.
6. Install the required packages for Python in a new terminal, with the command below: pip install numpy pandas matplotlib
7. Run the Python script. It should download the figures in a folder titled "figures".

In terms of data processing, the C++ code cleans up the datasets itself: it removes any values that are NA (no value saved) and selects only the columns necessary. It also skips empty lines or rows that cannot be converted to numbers to avoid error or interference with the calculations. 

It calculates statistics such as mean (both overall and subsample) root length, relative root length for each treatment, standard error, and subsampling uncertainties, as well as Pearson's correlation (_r_ value).

**Output**

After running both scripts, the user should have the .CSV files (specified above) from the C++ code and the figures from the Python code saved to their device.

**Report**

The report in this repository explains the project in more detail. Resources are linked in the bibliography for more information. 


