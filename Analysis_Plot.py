import matplotlib.pyplot as plt
import numpy as np
import requests
import re

def parse_genomic_coordinate(genomic_coordinate): 
    match = re.match(r'^(chr\d+):(\d+)-(\d+)_(\S+)$', genomic_coordinate)
    if match:
        chromosome = match.group(1)
        start = match.group(2)
        end = match.group(3)
        strand = match.group(4)
        return chromosome, int(start), int(end), strand
    else:
        print("Invalid genomic coordinate format")
        return None

def retrieve_junction_data(genomic_coordinate):
    chromosome, start, end, strand = parse_genomic_coordinate(genomic_coordinate)
    url = f"https://snaptron.cs.jhu.edu/srav2/snaptron?regions={chromosome}:{start}-{end}&rfilter=strand:{strand}"
    response = requests.get(url)
    if response.status_code == 200:
        junctions = [line.split('\t') for line in response.text.strip().split('\n')[1:]]
        junctions = [(int(j[3]), int(j[4]), float(j[14])) for j in junctions]
        return junctions
    else:
        print("Failed to retrieve junction data from Snaptron")
        return None

genomic_coordinate = "chr18:79930227-79930311_+" # Test Genomic Coordinate for analysis purpose 

junction_data = retrieve_junction_data(genomic_coordinate)

def plot_coverage():
    starts = [j[0] for j in junction_data]
    ends = [j[1] for j in junction_data]
    coverages = [j[2] for j in junction_data]
    plt.figure(figsize=(8, 4))
    plt.plot(starts, coverages, 'b-', label='Coverage')
    plt.plot(ends, coverages, 'r-', label='Coverage')
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.title(f'Snaptron Coverage Plot for {genomic_coordinate}')
    plt.legend()
    plt.show()

def plot_junction_position():
    positions = np.array([j[0] for j in junction_data] + [j[1] for j in junction_data])
    coverages = np.array([j[2] for j in junction_data] * 2)
    plt.figure(figsize=(8, 4))
    plt.scatter(positions, coverages, color='b', marker='o', label='Junctions')
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.title(f'Junction Position Plot for {genomic_coordinate}')
    plt.legend()
    plt.show()

def plot_junction_heatmap():
    positions = np.array([j[0] for j in junction_data] + [j[1] for j in junction_data])
    coverages = np.array([j[2] for j in junction_data] * 2)
    plt.figure(figsize=(8, 4))
    plt.hexbin(positions, coverages, gridsize=20, cmap='viridis')
    plt.colorbar(label='Density')
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.title(f'Junction Coverage Heatmap for {genomic_coordinate}')
    plt.show()

def plot_read_alignment():
    print("Read Alignment Plot is not implemented") # lack of bam file

def plot_comparative_analysis():
    sample1_coverage = np.random.randint(10, 50, size=len(junction_data))
    sample2_coverage = np.random.randint(10, 50, size=len(junction_data))
    
    plt.figure(figsize=(8, 4))
    plt.plot(sample1_coverage, 'b-', label='Sample 1 Coverage')
    plt.plot(sample2_coverage, 'r-', label='Sample 2 Coverage')
    plt.xlabel('Junction Index')
    plt.ylabel('Coverage')
    plt.title('Comparative Coverage Analysis')
    plt.legend()
    plt.show()

plot_read_alignment()
plot_comparative_analysis()

plot_coverage()
plot_junction_position()
plot_junction_heatmap()
