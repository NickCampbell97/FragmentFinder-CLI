import csv
from datetime import datetime
import os
from scipy.signal import find_peaks
import numpy as np
import zipfile
from collections import OrderedDict


# check all values surrounding true peak index to see if they should be included in the peak
def val_check(arr, index, threshold): 
    result = []
    right_indeces = []
    left_indeces = []

    for i in range((index-1), 1):
        if arr[i] > threshold:
            left_indeces.append(i)
        else:
            break
    for i in range((index+1), len(arr)):
        if arr[i] > threshold:
            if len(right_indeces) < 10:
                right_indeces.append(i)
            else:
                break
        else:
            break
    if not right_indeces:
        right_indeces.append(index)
    if not left_indeces:
        left_indeces.append(index)
        
    result.append(min(left_indeces))
    result.append(max(right_indeces))
    return result

# get the peak data for each of the alignment vectors and write to csv file
def get_peaks(results_dict, original_record_dict, fasta_read_counts, fasta_file_name, db_file_name, session):
    peaks_list = []
    out_file = f"output/{session}/results/{fasta_file_name}({db_file_name}).csv"
    fname_col = f'{fasta_file_name}({db_file_name})'
    for key, val in results_dict.items():
        standard_deviation = np.std(val)
        peaks, _ = find_peaks(val, prominence=standard_deviation*0.85, distance=18)
        t = val
        for i in peaks:
            rpm = (val[i]/fasta_read_counts)*1000000
            round_rpm = round(rpm, 4)
            if round_rpm > 5.0:
                peak_slice = val_check(val, i, (val[i]-(standard_deviation*0.1)))
                peak_start = peak_slice[0]-1
                peak_end = peak_slice[-1] + 17
                peak_string = original_record_dict[key][peak_start:peak_end]
                temp = [key, round_rpm, peak_start, peak_end, peak_string, '', fname_col]
                peaks_list.append(temp)
            else:
                continue

    with open(out_file, 'w') as f:
        fields = ['miRNA ID', 'Reads per Million', 'Peak Start', 'Peak End', 'Peak Data', 'Total Reads in NGS File', 'File Name']
        writer = csv.writer(f)
        writer.writerow(fields)
        writer.writerow(['', '', '', '', '', fasta_read_counts])
        for i in peaks_list:
            writer.writerow(i)


# get DEV in full peak width format
def regen_vector(val, fasta_read_count):
    t = val.copy()
    standard_deviation = np.std(val)
    peaks, _ = find_peaks(val, prominence=standard_deviation*0.85, distance=20)
    for peak in peaks:
        if ((val[peak] / fasta_read_count) * 1000000) >= 5.0:
            peak_slice = val_check(val, peak, (val[peak]-(standard_deviation*0.1)))
            peak_start = peak_slice[0]
            peak_end = (peak_slice[-1]) + 18
            t[peak_start:peak_end] = val[peak]
        else:
            continue
    return t

'''
    delete all of the temp image files in the temp directory
        - addresses an odd bug with adding a different jpg image with the 
          same file name to pdf with ReportLab
'''
def delete_temp_jpgs(directory_path):
    try:
        files = os.listdir(directory_path)
        for file in files:
            file_path = os.path.join(directory_path, file)
            os.remove(file_path)
    except Exception as e:
        print(f"An error occurred: {e}")

'''
    cut uid off at first instance of '|' to make fit pdf page
        - all of the id's are still unique after this operation
'''
def delimit_uid(uid_str, char):
    result = ''
    found = False
    for i in range(len(uid_str)):
        if uid_str[i] == char:
            found = True
        if not found:
            result += uid_str[i]
    return result

# get current date and format into 'month day, year'
def format_date():
    current_date = datetime.now().strftime('%B %d, %Y')
    return current_date

# format time from seconds into min:sec
def format_time(secs):
    minutes, secs = divmod(secs, 60)
    return f"{int(minutes):02d}:{int(secs):02d}"

# format an integer number into a string for readability
def format_large_ints(number):
    return '{:,}'.format(number)

# convert file size in bytes to either KB, MB, GB
def convert_bytes(byte_size, unit='auto', precision=1):
    units = {'B': 1, 'KB': 1000, 'MB': 1000 ** 2, 'GB': 1000 ** 3}
    if byte_size < 1:
        return '0 MB'
    if unit == 'auto':
        for u in ['GB', 'MB', 'KB', 'B']:
            if byte_size >= units[u]:
                unit = u
                break
    if unit not in units:
        print('Invalid file size. Input files must be less than 1 TB')
    size_in_unit = byte_size / units[unit]
    formatted_size = f"{size_in_unit:.{precision}f}"

    return f"{formatted_size} {unit}"

# get the size of file
def f_size(path):
    num_bytes = os.path.getsize(path)
    s_string = convert_bytes(num_bytes)
    return s_string

# return file name from file path
def fname_from_path(path):
    base = os.path.basename(path)
    split_name = os.path.splitext(base)
    return split_name[0]

# compress pdf file to zip for storage space
def compress_pdf(output_zip, output_pdf, f_name):
    with zipfile.ZipFile(output_zip, 'w', zipfile.ZIP_DEFLATED) as zipf:
        zipf.write(output_pdf, arcname=f"{f_name}.pdf")
        os.remove(output_pdf)

# compress pickle file to zip for storage space
def compress_pickle(output_z, pick_f, pick_name):
    with zipfile.ZipFile(output_z, 'w', zipfile.ZIP_DEFLATED) as f:
        f.write(pick_f, arcname=f'{pick_name}.pickle')
        os.remove(pick_f)

# unzip pickle
def unzip_pickle(zip_path, out_folder):
    with zipfile.ZipFile(zip_path, 'r') as f:
        f.extractall(out_folder)

# separate alignment vectors by whether they have a single or multiple peaks
def split_by_peak_count(results, count):
    single_peaks = {}
    mult_peaks = {}
    for uid, vector in results.items():
        standard_deviation = np.std(vector)
        peaks, _ = find_peaks(vector, prominence=standard_deviation*0.9, distance=20)
        temp = []
        for x in range(len(peaks)):
            if ((vector[peaks[x]] / count) * 1000000.0) >= 5.0:
                temp.append(peaks[x])
            else:
                pass
        if len(temp) == 1:
            single_peaks[uid] = [vector, temp]
        elif len(temp) > 1:
            mult_peaks[uid] = [vector, temp]
        else:
            continue
    return single_peaks, mult_peaks


# create new directory
def create_dir(dir_path):
    try:
        os.mkdir(dir_path)
    except FileExistsError:
        print(f"Directory '{dir_path}' already exists.")

# create necessary subdirectories on start of session
def config_sesh_dir(base_path):
    create_dir(f'{base_path}/results')
    create_dir(f'{base_path}/reports')
    create_dir(f'{base_path}/temp')
    create_dir(f'{base_path}/data')

# check if path is directory or not
def is_directory(file_path):
    return os.path.isdir(file_path)

# check if given NGS file is valid fasta or fastq
def is_valid_file_format(file_path, ext_list):
    _, file_extension = os.path.splitext(file_path)
    return file_extension.lower() in [ext.lower() for ext in ext_list]

# return all files with the specific extensions in a given directory
def get_file_list(directory_path, ext_list):
    file_paths = []
    try:
        all_files = os.listdir(directory_path)
        file_paths = [os.path.join(directory_path, file) for file in all_files 
                      if any(file.lower().endswith(ext) for ext in ext_list)]
    except OSError as e:
        print(f"Error: {e}")
    return file_paths

'''
    for generating reports of old session 
        - make a map of ints to data file paths for option input
'''
def make_file_map(directory_path):
    if not os.path.isdir(directory_path):
        print(f"The specified path '{directory_path}' is not a valid directory.")
        return
    file_list = os.listdir(directory_path)
    file_map = OrderedDict()
    count = 1
    for file_name in file_list:
        file_path = os.path.join(directory_path, file_name)
        parts = file_name.split('(', 1)
        f_name = os.path.splitext(file_name)[0]
        if len(parts) > 1:
            first_part = parts[0]
            file_map[count] = [first_part, file_path, f_name]
            count += 1
    return file_map


