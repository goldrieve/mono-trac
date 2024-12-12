import csv
import sys

def extract_averages(file_path):
    cleaned_header = [file_path.split(':')[0]]
    results = [["isolate", f"{cleaned_header[0]}_T", f"{cleaned_header[0]}_C", f"{cleaned_header[0]}_A", f"{cleaned_header[0]}_G"]]
    try:
        with open(file_path, "r") as file:
            lines = file.readlines()

        for i in range(len(lines)):
            if lines[i].startswith("#"):
                name = lines[i].strip().split(":")[1].strip()
                if i + 4 < len(lines) and "Average" in lines[i + 4]:
                    average_line = lines[i + 4].replace("Average", "").strip()
                    average_values = [value.split(":")[1] for value in average_line.split()]
                    results.append([name] + average_values)
    except FileNotFoundError:
        print(f"File {file_path} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    return results

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python extract.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    results = extract_averages(file_path)
    if results:
        with open("output.csv", "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(results)
