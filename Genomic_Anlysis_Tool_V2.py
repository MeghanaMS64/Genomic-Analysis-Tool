import sys
import requests
import re
import pyBigWig

# Function to parse genomic coordinates
def parse_genomic_coordinate(genomic_coordinate):
    # Using regular expression to match the genomic coordinate pattern
    match = re.match(r'^(chr\d+):(\d+)-(\d+)_(\S+)$', genomic_coordinate)
    if match:
        # Extracting chromosome, start, end, and strand from the match
        chromosome = match.group(1)  # Chromosome identifier (e.g., "chr1")
        start = match.group(2)       # Start position of the genomic region
        end = match.group(3)         # End position of the genomic region
        strand = match.group(4)      # Strand orientation (e.g., "+" or "-")
        return chromosome, int(start), int(end), strand
    else:
        print("Invalid genomic coordinate format")
        sys.exit(1)

# Function to retrieve junction data from Snaptron
def retrieve_junction_data(genomic_coordinate, filter_condition=None):
    chromosome, start, end, strand = parse_genomic_coordinate(genomic_coordinate)
    # Constructing the Snaptron API URL based on the genomic coordinate and optional filter condition
    url = f"https://snaptron.cs.jhu.edu/srav2/snaptron?regions={chromosome}:{start}-{end}&rfilter=strand:{strand}"
    if filter_condition:
        url += f"&rfilter={filter_condition}"
    # Sending GET request to the Snaptron API
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        print("Failed to retrieve junction data from Snaptron")
        sys.exit(1)

# Function to identify exon edges from junction data
def identify_exon_edges(genomic_coordinate, junction_data):
    # Parsing junction data and extracting relevant information
    junctions = [line.split('\t') for line in junction_data.strip().split('\n')[1:]]
    junctions = [(int(j[3]), int(j[4]), float(j[14])) for j in junctions]
    junctions.sort(key=lambda x: x[0])
    exon_edges = []
    current_start = None
    current_end = None
    # Identifying exon edges based on coverage threshold
    for start, end, coverage in junctions:
        if coverage > 0.5:
            if current_start is None:
                current_start = start
            current_end = end
        elif current_start is not None:
            exon_edges.append((current_start, current_end))
            current_start = None
            current_end = None
    if current_start is not None:
        exon_edges.append((current_start, current_end))
    return exon_edges

# Function to measure transcription abundance from BigWig files
def measure_transcription_abundance(genomic_coordinate, bigwig_files):
    chromosome, start, end, _ = parse_genomic_coordinate(genomic_coordinate)
    abundance = []
    # Looping through each BigWig file
    for bw_file in bigwig_files:
        bw = pyBigWig.open(bw_file)
        if bw is not None:
            # Retrieving values from the BigWig file for the given genomic region
            values = bw.values(chromosome, start, end)
            if values:
                # Summing up non-None values to calculate abundance
                abundance.append(sum(filter(lambda x: x is not None, values)))
            else:
                abundance.append(0)
        else:
            print(f"Failed to open BigWig file: {bw_file}")
            sys.exit(1)
    return abundance

# Function to analyze genomic data based on user-selected option
def analyze_genomic_data(genomic_coordinate, analyze_option, filter_condition=None, bigwig_file=None):
    if analyze_option == "Exon Edge Finder":
        try:
            # Retrieving junction data from Snaptron and identifying exon edges
            junction_data = retrieve_junction_data(genomic_coordinate, filter_condition)
            exon_edges = identify_exon_edges(genomic_coordinate, junction_data)
            print("Exon Edges:")
            for edge in exon_edges:
                print(f"Start: {edge[0]}, End: {edge[1]}")
        except Exception as e:
            print(f"An error occurred: {str(e)}")
    elif analyze_option == "Genomic Abundance Analyzer":
        try:
            # Measuring transcription abundance from BigWig files
            abundance = measure_transcription_abundance(genomic_coordinate, [bigwig_file])
            print(f"Transcription Abundance:\n{bigwig_file}: {abundance[0]}")
        except Exception as e:
            print(f"An error occurred: {str(e)}")

# Main function for command-line mode
if len(sys.argv) == 2 and sys.argv[1] == "cmd":
    analyze_option = input("Choose Analysis Option (Exon Edge Finder / Genomic Abundance Analyzer): ")
    genomic_coordinate = input("Enter Genomic Coordinate (chr#:start-end_strand): ")
    if analyze_option == "Exon Edge Finder":
        filter_condition = input("Enter Filter Condition (optional): ")
        analyze_genomic_data(genomic_coordinate, analyze_option, filter_condition)
    elif analyze_option == "Genomic Abundance Analyzer":
        bigwig_file = input("Enter BigWig File Path: ")
        analyze_genomic_data(genomic_coordinate, analyze_option, bigwig_file=bigwig_file)
else:
    import tkinter as tk
    from tkinter import messagebox
    from tkinter import filedialog
    import requests
    import re
    import sys
    import pyBigWig
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    import csv

        # Function to update UI based on selected analysis option
    def update_ui(selected_option):
        if selected_option == "Exon Edge Finder":
            filter_label.grid(row=2, column=0, padx=5, pady=5)
            filter_entry.grid(row=2, column=1, padx=5, pady=5)
            analyze_button.grid(row=3, column=0, columnspan=2, padx=5, pady=5)
            csv_button.grid(row=4, column=0, columnspan=2, padx=5, pady=5)
            bigwig_button.grid_forget()
            bigwig_entry.grid_forget()
        elif selected_option == "Genomic Abundance Analyzer":
            filter_label.grid_forget()
            filter_entry.grid_forget()
            csv_button.grid_forget()
            analyze_button.grid(row=5, column=0, columnspan=2, padx=5, pady=5)
            bigwig_button.grid(row=3, column=0, columnspan=2, padx=5, pady=5)
            bigwig_entry.grid(row=4, column=0, columnspan=2, padx=5, pady=5)


    def parse_genomic_coordinate(genomic_coordinate): 
        match = re.match(r'^(chr\d+):(\d+)-(\d+)_(\S+)$', genomic_coordinate)
        if match:
            chromosome = match.group(1)
            start = match.group(2)
            end = match.group(3)
            strand = match.group(4)
            return chromosome, int(start), int(end), strand
        else:
            messagebox.showerror("Error", "Invalid genomic coordinate format")
            sys.exit(1)

    def retrieve_junction_data(genomic_coordinate, filter_condition=None):
        chromosome, start, end, strand = parse_genomic_coordinate(genomic_coordinate)
        url = f"https://snaptron.cs.jhu.edu/srav2/snaptron?regions={chromosome}:{start}-{end}&rfilter=strand:{strand}"
        if filter_condition:
            url += f"&rfilter={filter_condition}"
        response = requests.get(url)
        if response.status_code == 200:
            return response.text
        else:
            messagebox.showerror("Error", "Failed to retrieve junction data from Snaptron")
            sys.exit(1)

    def identify_exon_edges(genomic_coordinate, junction_data):
        junctions = [line.split('\t') for line in junction_data.strip().split('\n')[1:]]
        junctions = [(int(j[3]), int(j[4]), float(j[14])) for j in junctions]
        junctions.sort(key=lambda x: x[0])
        exon_edges = []
        current_start = None
        current_end = None
        for start, end, coverage in junctions:
            if coverage > 0.5:
                if current_start is None:
                    current_start = start
                current_end = end
            elif current_start is not None:
                exon_edges.append((current_start, current_end))
                current_start = None
                current_end = None
        if current_start is not None:
            exon_edges.append((current_start, current_end))
        return exon_edges

    def plot_snaptron_data(genomic_coordinate, ax):
        junction_data = retrieve_junction_data(genomic_coordinate)
        junctions = [line.split('\t') for line in junction_data.strip().split('\n')[1:]]
        junctions = [(int(j[3]), int(j[4]), float(j[14])) for j in junctions]
        starts = [j[0] for j in junctions]
        ends = [j[1] for j in junctions]
        coverages = [j[2] for j in junctions]
        ax.plot(starts, coverages, 'b-', label='Coverage')
        ax.plot(ends, coverages, 'r-', label='Coverage')
        ax.set_xlabel('Position')
        ax.set_ylabel('Coverage')
        ax.set_title(f'Snaptron Coverage Plot for {genomic_coordinate}')
        ax.legend()

    def analyze_genomic_data(genomic_coordinates):
        output_window = tk.Toplevel(root)
        output_window.title("Analysis Results")
        output_text = tk.Text(output_window, wrap=tk.WORD, width=80, height=20)
        output_text.pack(fill=tk.BOTH, expand=True)

        # Create window for plots
        plots_window = tk.Toplevel(root)
        plots_window.title("Snaptron Coverage Plots")

        # Create subplots
        num_plots = len(genomic_coordinates)
        num_cols = min(num_plots, 3)
        num_rows = (num_plots - 1) // 3 + 1
        fig, axs = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows))

        for i, genomic_coordinate in enumerate(genomic_coordinates):
            try:
                if analyze_option.get() == "Exon Edge Finder":
                    junction_data = retrieve_junction_data(genomic_coordinate)
                    exon_edges = identify_exon_edges(genomic_coordinate, junction_data)
                    
                    output_text.insert(tk.END, f"Exon Edges for {genomic_coordinate}:\n")
                    for edge in exon_edges:
                        output_text.insert(tk.END, f"Start: {edge[0]}, End: {edge[1]}\n")
                elif analyze_option.get() == "Genomic Abundance Analyzer":
                    file_path = filedialog.askopenfilename(title="Select BigWig file", filetypes=[("BigWig files", "*.bw")])
                    abundance = measure_transcription_abundance(genomic_coordinate, [file_path])
                    output_text.insert(tk.END, f"Transcription Abundance for {genomic_coordinate}:\n{file_path}: {abundance[0]}\n")
                
                # Plot Snaptron data
                row = i // num_cols
                col = i % num_cols
                plot_snaptron_data(genomic_coordinate, axs[row, col])
            except Exception as e:
                output_text.insert(tk.END, f"Error analyzing {genomic_coordinate}: {str(e)}\n")

        # Adjust layout of subplots
        plt.tight_layout()
        
        # Embed plots in Tkinter window
        canvas = FigureCanvasTkAgg(fig, master=plots_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # Create GUI
    root = tk.Tk()
    root.title("Genomic Analysis Tool")

    # Analyze Option
    analyze_option = tk.StringVar()
    analyze_option.set("Exon Edge Finder")
    analyze_option_label = tk.Label(root, text="Choose Analysis Option:")
    analyze_option_label.grid(row=0, column=0, padx=5, pady=5)
    analyze_option_menu = tk.OptionMenu(root, analyze_option, "Exon Edge Finder", "Genomic Abundance Analyzer", command=update_ui)
    analyze_option_menu.grid(row=0, column=1, padx=5, pady=5)

    # Genomic Coordinate
    genomic_label = tk.Label(root, text="Genomic Coordinate (chr#:start-end_strand):")
    genomic_label.grid(row=1, column=0, padx=5, pady=5)
    genomic_entry = tk.Entry(root, width=50)
    genomic_entry.grid(row=1, column=1, padx=5, pady=5)

    # Filter Condition
    filter_label = tk.Label(root, text="Filter Condition (optional):")
    filter_entry = tk.Entry(root, width=50)

    bigwig_button = tk.Button(root, text="Choose BigWig File", command=lambda: bigwig_entry.insert(tk.END, filedialog.askopenfilename(title="Select BigWig file", filetypes=[("BigWig files", "*.bw")])))
    bigwig_entry = tk.Entry(root, width=50)

    # Analyze Button
    analyze_button = tk.Button(root, text="Analyze", command=lambda: analyze_genomic_data([genomic_entry.get()]))
    analyze_button.grid(row=3, column=0, columnspan=2, padx=5, pady=5)

    # CSV File Button
    def analyze_csv_file():
        file_path = filedialog.askopenfilename(title="Select CSV file", filetypes=[("CSV files", "*.csv")])
        if file_path:
            with open(file_path, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                genomic_coordinates = [f"{row['chromosome']}:{row['start']}-{row['end']}_{row['strand']}" for row in reader]
                analyze_genomic_data(genomic_coordinates)

    csv_button = tk.Button(root, text="Analyze CSV File", command=analyze_csv_file)
    csv_button.grid(row=4, column=0, columnspan=2, padx=5, pady=5)

    update_ui("Exon Edge Finder")

    root.mainloop()
