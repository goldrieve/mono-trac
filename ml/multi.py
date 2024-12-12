import csv
import sys
import argparse

def extract_averages(file_path):
    cleaned_header = file_path.split(':')[0]
    results = {}
    try:
        with open(file_path, "r") as file:
            lines = file.readlines()

        for i in range(len(lines)):
            if lines[i].startswith("#"):
                name = lines[i].strip().split(":")[1].strip()
                if i + 4 < len(lines) and "Average" in lines[i + 4]:
                    average_line = lines[i + 4].replace("Average", "").strip()
                    average_values = [value.split(":")[1] for value in average_line.split()]
                    results[name] = average_values
    except FileNotFoundError:
        print(f"File {file_path} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    return cleaned_header, results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract averages from files and output to CSV.")
    parser.add_argument("file_paths", nargs="+", help="Paths to the input files.")
    parser.add_argument("-o", "--output", default="output.csv", help="Path to the output CSV file.")
    args = parser.parse_args()

    file_paths = args.file_paths
    output_file = args.output
    all_results = {}
    headers = ["isolate"]

    for file_path in file_paths:
        cleaned_header, results = extract_averages(file_path)
        if not all_results:
            # Initialize the headers and results dictionary
            headers.extend([f"{cleaned_header}_T", f"{cleaned_header}_C", f"{cleaned_header}_A", f"{cleaned_header}_G"])
            for isolate, averages in results.items():
                all_results[isolate] = averages
        else:
            # Extend the headers and update the results dictionary
            headers.extend([f"{cleaned_header}_T", f"{cleaned_header}_C", f"{cleaned_header}_A", f"{cleaned_header}_G"])
            for isolate, averages in results.items():
                if isolate not in all_results:
                    all_results[isolate] = [""] * (len(headers) - 5)  # Fill with empty strings for previous columns
                all_results[isolate].extend(averages)

    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        for isolate, averages in all_results.items():
            writer.writerow([isolate] + averages)