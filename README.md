# **Annovar Analysis Tool**

### Installation

The analysis tool can be installed from GitHub by directly downloading the raw file content of [`Annovar_Analysis_Sort`](https://github.com/NicoHadas/Annovar_Tool/blob/main/Python/Annovar_Analysis_Sort.py)

### Quick Guide

`Annovar_Analysis_Sort` analyzes Annovar ouput data and identifies the pathogenic mutations

```{bash}
# Run the Annovar Tool in Command Line
Annovar_Analysis_Sort.py Annovar_Analysis_Sort Input_File_Path Output_File_Path
```
For Example
```{bash}
Annovar_Analysis_Sort.py Annovar_Analysis_Sort /home/User/Annovar_Input.xlsx /home/User/Annovar_Output.xlsx
```
## Dependencies 

[Pandas](https://pandas.pydata.org/docs/getting_started/install.html)

[argparse](https://docs.python.org/3/library/argparse.html)
