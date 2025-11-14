import sys
import csv

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py input.tsv output.tsv")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    with open(input_file, newline='', encoding='utf-8') as infile:
        reader = list(csv.reader(infile, delimiter='\t'))

    # Insert 'x' at the beginning of the first row
    if reader:
        reader[0].insert(0, 'x')

    # Write back to output file
    with open(output_file, 'w', newline='', encoding='utf-8') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerows(reader)

if __name__ == "__main__":
    main()
