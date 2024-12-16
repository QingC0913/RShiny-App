"""
1. download series matrix files from NCBI
2. extract, convert/save as CSV
3. use this code to process CSV to get only sample characteristics
4. use the output file in R for further processing and analysis
"""
# selecting rows 
# file = "GSE64810_series_matrix.csv"
# output_file = "processed_series_matrix.csv"
# csv = open(file, "r")
# output_csv = open(output_file, "w")
rowname = "!Sample_characteristics_ch1"
# for line in csv: 
#     vals = line.split(",")
#     if vals[0] == rowname:
#         output_csv.write(line)
# csv.close()
# output_csv.close()

# editing row values
file = "processed_series_matrix.csv"
output_file = "processed_samples.csv"
csv = open(file, "r")
output_csv = open(output_file, "w")
for line in csv: 
    vals = line.split(",")
    i = 50 if (vals[1] == "") else 1
    value = vals[i]
    print(value)
    key = value.split(":")[0]
    vals = ["NA" if v == "" else v for v in vals][1:]
    line = ",".join(vals).replace(str(key) + ": ", "")
    line = (key + "," + line).replace(" ", "_").replace("-", "_")
    print(line)
    output_csv.write(line)

csv.close()
output_csv.close()
