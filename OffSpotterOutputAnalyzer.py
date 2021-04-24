from collections import Counter
import csv

""" This suite of code analyzes OffSpotter_remote outputs.
It requires a list of all text file names, one in each line, in gRNA_list.txt.
It creates OffSpotter_summary.csv containing summary information for each
candidate gRNA for further analysis.

SIMPLY RUN:
offSpotterAnalyzer()
"""


def offSpotterAnalyzer():
    # Make sure to include gRNA_list.txt in the same directory as OffSpotterOutputAnalyzer.py
    
    def read_file(file_name):
        # Reads gRNA text files and returns each line as an element in a list
        f = open(file_name, "r")
        copy = f.read()
        f.close()

        return copy.split("\n")

    
    def columnToList(data, col):
        # Takes a list with rows of sRNA data as its elements and extract the information in column col
        ret_list = []
        
        for row in data:
            if row != "":
                # This conditional is necessary to avoid reading the last empty line in the gRNA files
                split_row = row.split("\t")
                ret_list.append(split_row[col])

        return ret_list
    

    gRNA_list = read_file("gRNA_list.txt")
    gRNA_list.pop()
    summary = []

    for gRNA in gRNA_list:
        # Scans through each file individually
        gRNA_data = read_file(gRNA + ".txt")
        gRNA_data = gRNA_data[4:]

        # Extracts chromosome # information for each hit
        chromo = columnToList(gRNA_data, 0)
        # Extracts the number of mismatched bases for each hit
        numMismatch = columnToList(gRNA_data, 6)
        
        summary.append((gRNA,
                         numMismatch.count("0"),
                         numMismatch.count("1"),
                         numMismatch.count("2"),
                         numMismatch.count("3"),
                         numMismatch.count("4"),
                         numMismatch.count("5"),
                         len(gRNA_data),
                         len(Counter(chromo))))

    # Write CSV file
    csvwrite = open("OffSpotter_summary.csv", "wb")
    csv_writer = csv.writer(csvwrite, delimiter = ',')
    # n MM = total number of hits with n mismatches
    csv_writer.writerow(["Sequence", "0 MM", "1 MM", "2 MM", "3 MM", "4 MM",
                         "5 MM", "Total MM", "Unique chromosome hits"])

    for row in summary:
        csv_writer.writerow(row)

    print("OffSpotter_summary.csv PRINTED.")
                         
