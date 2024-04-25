# Genomic Analysis Tool

## Overview
The Genomic Analysis Tool is a Python script made to identify exon edges and measuring transcription abundance. 

## Features
- **Exon Edge Finder**: Identifies exon edges within a specified genomic coordinate range based on junction data retrieved from Snaptron.
- **Genomic Abundance Analyzer**: Measures transcription abundance within a specified genomic coordinate range using BigWig files.
- **GUI And CLI**: A User-Friendly interface and a Command Line Interface. 
## Requirements
- Python 3.x
- Requests library (`pip install requests`)
- pyBigWig library (`pip install pyBigWig`)

## Usage
### Command-line Interface (CLI) Mode
Run the script with the following command:
```python Genomic_Analysis_Tool.py cmd```
### Graphical User Interface (GUI) Mode
Run the script without any arguments:
```python Genomic_Analysis_Tool.py```
## Usage Examples
### Exon Edge Finder
```
python Genomic_Analysis_Tool.py cmd
Choose Analysis Option (Exon Edge Finder / Genomic Abundance Analyzer): Exon Edge Finder
Enter Genomic Coordinate (chr#:start-end_strand): chr1:1000-2000_+
Enter Filter Condition (optional): coverage>0.5
``` 
### Genomic Abundance Analyzer
```
python Genomic_Analysis_Tool.py cmd
Choose Analysis Option (Exon Edge Finder / Genomic Abundance Analyzer): Genomic Abundance Analyzer
Enter Genomic Coordinate (chr#:start-end_strand): chr1:1000-2000_+
Enter BigWig File Path: sample.bw
```
## Authors
- Meghana Sripathi

## Acknowledgments
- Snaptron API for providing junction data
- PyBigWig library for handling BigWig files
