# **Annovar Analysis Tool**

### Installation

The analysis tool can be installed from GitHub by directly downloading the raw file content of [`Annovar_Analysis_Sort`](https://github.com/NicoHadas/Annovar_Tool/blob/main/Python/Annovar_Analysis_Sort.py)

### Quick Guide

`Annovar_Analysis_Sort` analyzes Annovar ouput data and identifies the pathogenic mutations

```{python}
# Import Pandas
import pandas as pd

# Import Annovar_Analysis_Sort.py
import Annovar_Analysis_Sort.py as AA

# Run the Annovar Tool
AA.Annovar_Analysis_Sort("Input_File_Path", "Output_File_Path")
```

## Dependencies 

[Pandas](https://pandas.pydata.org/docs/getting_started/install.html)
