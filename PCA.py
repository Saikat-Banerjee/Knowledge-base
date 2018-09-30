
data_path = "C:\Users\sktba\Desktop\Elucidata\Assignment-2_gene_data- (2).csv"
import csv
with open(data_path, "rb") as f:
    reader = csv.reader(f)
    i = reader.next()
    rest = [row for row in reader]


