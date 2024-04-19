import pandas as pd
from Bio import SeqIO
import glob
import os
import ast
import re

def get_files(dir_name):
    joined_files = os.path.join(dir_name, "*.csv")
    joined_list = glob.glob(joined_files)
    return joined_list

def show_files(lst):
    for file in lst:
        print(file)

def get_list_count(lst):
    count = 0
    for file in lst:
        count += 1
    return count

def get_new_out_file_name(lst):
    s = lst[0]
    lst = re.findall(f'\(.*?\)', s)
    new_file_name = str(lst[0])
    return new_file_name


def merge_files(dir_name):
    joined_files = os.path.join(dir_name, "*.csv")
    joined_list = glob.glob(joined_files)
    df = pd.concat(map(pd.read_csv, joined_list), ignore_index=True)
    df.sort_values(['miRNA ID', 'Peak Start'], ascending=True)
    df.to_csv(f'{dir_name}/t1.csv', index=False)

# Define a function to check if two values are within a given range
def within_range(val1, val2, range_value):
    return abs(val1 - val2) <= range_value

def init_group(dir_name):
    csv_file = f'{dir_name}/t1.csv'
    df = pd.read_csv(csv_file)

    # Sort the DataFrame by 'Peak Start'
    df.sort_values(by=['miRNA ID', 'Peak Start'], inplace=True)
    #g = df.groupby('miRNA ID')

    # Initialize variables
    grouped_data = []
    current_group = []
    previous_start = None
    previous_ID = None
    current_group_start = 0

    # Iterate through the DataFrame rows
    for index, row in df.iterrows():
        peak_start = row['Peak Start']
        peak_end = row['Peak End']
        id = row['miRNA ID']

        if (previous_ID is None) or (id != previous_ID) or (current_group_start is None) or (not within_range(peak_start, previous_start, 10)) or (peak_end - current_group_start > 40):
            if current_group:
                grouped_data.append(current_group)
            current_group = [row]
            current_group_start = peak_start
        else:
            current_group.append(row)
            
        previous_start = peak_start
        previous_ID = id

    # Append the last group
    if current_group:
        grouped_data.append(current_group)

    output_csv = f'{dir_name}/grouped_by_peak.csv'
    with open(output_csv, 'w') as f:
        f.write('Group,miRNA ID,Peak Start,Peak End,Reads per Million,File Name\n')
        for idx, group in enumerate(grouped_data, start=1):
            for item in group:
                f.write(f'Group {idx},{item["miRNA ID"]},{item["Peak Start"]},{item["Peak End"]}, {item["Reads per Million"]}, {item["File Name"]}\n')


def get_orig_strings(db):
    records = {}
    with open(db, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            s = str(record.seq)
            if len(s) > 0:
                records[record.id] = s
            else:
                continue
    return records


def final_group(group, arg2):
    
    min_start = group['Peak Start'].min()
    max_end = group['Peak End'].max()


    peak_data = arg2[group['miRNA ID'].iloc[0]][int(min_start):int(max_end)]
    
    reads_and_files = []
    for index, row in group.iterrows():
        reads_and_files.append({
            row['File Name']: row['Reads per Million']
        })
    
    return pd.DataFrame({
            'miRNA ID': [group['miRNA ID'].iloc[0]],
            'Peak Start': [min_start],
            'Peak End': [max_end],
            'Files : Reads': [reads_and_files],
            'Peak Data': [peak_data]
        })


def get_final_group(dir_name, hmap):
    file = f'{dir_name}/grouped_by_peak.csv'
    df2 = pd.read_csv(file)
    # set final_frame grouped_df and call second group function
    grouped = df2.groupby(['Group', 'miRNA ID']).apply(final_group, arg2 = hmap)
    grouped.rename_axis(index={'miRNA ID': 'ID'}, inplace=True)

    grouped = grouped.sort_values(by=['miRNA ID', 'Peak Start'])

    # final_frame = grouped_df.groupby(['Peak Group', 'miRNA ID']).apply(final_group)
    grouped.to_csv(f'{dir_name}/files_and_reads.csv', index=False)

def export_final_csv(dir_name, count, fname):
    # Read the CSV file into a DataFrame
    csv_file = f'{dir_name}/files_and_reads.csv'  # Replace with your actual CSV file path
    original = pd.read_csv(csv_file)

    #original_df = final_frame

    # Convert the string representations of lists into actual lists
    original['Files : Reads'] = original['Files : Reads'].apply(ast.literal_eval)

    # Create a set of all keys from the dictionaries
    all_keys = set()
    for index, row in original.iterrows():
        reads_list = row['Files : Reads']
        for dictionary in reads_list:
            all_keys.update(dictionary.keys())

    # Initialize a dictionary to store the data
    data_dict = {key: [] for key in all_keys}

    # Iterate through each row
    for index, row in original.iterrows():
        reads_list = row['Files : Reads']
        
        # Initialize a dictionary for the row's data
        row_data = {key: 0 for key in all_keys}
        
        # Update the dictionary with actual read values
        for dictionary in reads_list:
            for key, value in dictionary.items():
                if key in row_data.keys():
                    if float(row_data[key]) < float(value):
                        row_data[key] = float(value)
                    else:
                        continue
                else:
                    row_data[key] = float(value)
        
        # Append values to the main data dictionary
        for key, value in row_data.items():
            data_dict[key].append(value)

    # Create a DataFrame from the extracted data
    read_df = pd.DataFrame(data_dict)
    extracted_col = original['Peak Data']
    original = original.drop(columns=['Peak Data'])
    final_df = pd.concat([original, read_df], axis=1)
    final_df = final_df.drop(columns=['Files : Reads'])
    final_df.insert(count+3, 'Peak Data', extracted_col)
    final_df = final_df.sort_values(by=['miRNA ID', 'Peak Start'])
    os.makedirs(f'{dir_name}/Joined_Results')
    final_df.to_csv(f'{dir_name}/Joined_Results/{fname}.csv', index=False)



def main():

    while True:
        directory_name = input("Enter Full Path to Directory Containing Files to Merge: ")
        if directory_name == '':
            print('Invalid - Please Enter Directory Path: ')
        if not directory_name == '':
            file_list = get_files(directory_name)
            show_files(file_list)
            break

    db_file_path = input("Enter Full Path to Database File Used for Sequence Alignment: ")

    print()
    print('Merging and Grouping Files...')
    print()

    out_file = get_new_out_file_name(file_list)
    file_count = get_list_count(file_list)

    merge_files(directory_name)

    init_group(directory_name)

    hmap = get_orig_strings(db_file_path)

    get_final_group(directory_name, hmap)

    export_final_csv(directory_name, file_count, out_file)

    print(f"Files Merged and Output Exported to 'Joined_Outputs/{out_file}.csv'")
    os.remove(f'{directory_name}/files_and_reads.csv')
    os.remove(f'{directory_name}/grouped_by_peak.csv')
    os.remove(f'{directory_name}/t1.csv')

if __name__ == "__main__":
    main()
