import Procedure
import Util
import os
import time
from Reports import create_canvas, plot_mult_peaks, plot_single_peaks, cleanup
import pickle
from collections import OrderedDict
import shutil


# get count of mapped miRNA's
def get_mapped_count(known_match_file):
    count = 0
    with open(known_match_file, 'r') as f:
        for line in f:
            count += 1
    return count

# get database file
def get_db_path():
    while True:
        db_choice = input("To get started, would you like to use the \ndefault (FF_human_db.txt) database? (y/n) ").strip()
        if db_choice in {'y', 'Y', 'n', 'N'}:
            break
        else:
            print('Invalid Input. Valid Inputs Include y/n (yes or no)')
    if db_choice in {'n', 'N'}:
        while True:
            path = input('\nPlease enter the full path to the \nDatabase file you would like to use: ').strip()
            if Util.is_valid_file_format(path, ['.txt', '.fasta']):
                break
            else:
                print(f'\nInvalid Database file: {path}\nThe Database file must either be a TXT or FASTA file.\n')
    elif db_choice in {'y', 'Y'}:
        path = '.defaultdb/FF_human_db.txt'
    return path


# generate report for one NGS/DB file pair from the results list that was stored in pickle
def generate_selection_report(result_list):
    results = result_list[0]
    total_reads = result_list[1]
    db_p = result_list[2]
    ngs_p = result_list[3]
    full_db_seqs = result_list[4]
    time_string = result_list[5]
    sesh = result_list[6]
    output_fname = result_list[7]
    single_peaks, mult_peaks = Util.split_by_peak_count(results, total_reads)
    canv = create_canvas(db_p, ngs_p, sesh)
    plot_single_peaks(canv, single_peaks, full_db_seqs, total_reads, db_p, ngs_p, time_string, sesh)
    plot_mult_peaks(canv, mult_peaks, full_db_seqs, total_reads, sesh)
    cleanup(canv, sesh)
    out_z = f'output/{sesh}/reports/{output_fname}.zip'
    out_pdf = f'output/{sesh}/reports/{output_fname}.pdf'
    Util.compress_pdf(out_z, out_pdf, output_fname)


# run alignment algo and generate csv
def align(db_path, ngs_path, session):
    ngs_name = Util.fname_from_path(ngs_path)
    db_name = Util.fname_from_path(db_path)

    # separate the sequences and get read count value
    temp_matches = f'output/{session}/temp/matches.txt'
    start_preprocess = time.perf_counter()
    total_reads = Procedure.separate_seqs(db_path, ngs_path, temp_matches, 18)
    end_preprocess = time.perf_counter()
    pre_time = (end_preprocess - start_preprocess)

    # align sequences and get resulting dictionary
    start_process = time.perf_counter()
    results = Procedure.seq_aligner(db_path, temp_matches, 18)
    end_process = time.perf_counter()
    proc_time = (end_process - start_process)

    total_time = (pre_time + proc_time) # total time it took to align files
    time_string = Util.format_time(total_time)
    os.remove(temp_matches)

    full_db_seqs = Procedure.original_strings(db_path)
    
    Util.get_peaks(results, full_db_seqs, total_reads, ngs_name, db_name, session)

    output_fname = f'{ngs_name}({db_name})'

    print(f'Alignment Complete! The results file has been generated in\noutput/{session}/results/{output_fname}.csv\n')
    # return tuple => time, result list:
    #   the results list contains the results dictionary, read count, db and ngs file paths
    #   the original sequence dictionary from db file, time string, session id, and result file name
    return (total_time, [results, total_reads, db_path, ngs_path, full_db_seqs, time_string, session, output_fname])

# for aligning and generating report when there is just one NGS file in the session
def align_single_ngs(db_p, ngs_p, sesh):

    db_base = os.path.basename(db_p)
    ngs_base = os.path.basename(ngs_p)

    base_output_dir = f'output/{sesh}'
    Util.create_dir(base_output_dir)
    Util.config_sesh_dir(base_output_dir)

    print(f'\nAligning {db_base} and {ngs_base}...')

    t, res_list = align(db_p, ngs_p, sesh)

    time_string = res_list[5]
    print(f'\nTime Taken to Align (min:sec): {time_string}')

    output_fname = res_list[7]
    pickle_file = f'output/{sesh}/data/{output_fname}.pickle'

    # add the results to the pickle file in case they need to be accessed later
    with open(pickle_file, 'wb') as f:
        pickle.dump(res_list, f)

    Util.compress_pickle(f'output/{sesh}/data/{output_fname}.zip', pickle_file, output_fname)

    while True:
        pdf_choice = input('\nWould you like to generate a PDF \nreport of this FragmentFinder session? (y/n) ').strip()
        if pdf_choice in {'y', 'Y', 'n', 'N'}:
            break
        else:
            print('\nInvalid Input. Valid Inputs Include y or n (yes or no)')

    if pdf_choice in {'n', 'N'}:
        print('\nNo PDF report will be generated. \nThank you for using FragmentFinder!')
        shutil.rmtree(f'output/{sesh}/reports')
        shutil.rmtree(f'output/{sesh}/temp')


    elif pdf_choice in {'y', 'Y'}:
        print('\nReport generation in progress.\nThe PDF will be generated as a compressed ZIP File.\nPlease note that this process may take some time.')
        print('\nPress Ctrl c to stop this process at any time.')
        generate_selection_report(res_list)
        shutil.rmtree(f'output/{sesh}/temp')
        print(f'\nThe compressed report file has been sucessfully created in \noutput/{sesh}/reports/{res_list[7]}.zip\n\n')

# for aligning multiple NGS files per session
# each NGS and DB file pair gets their own result list (which is pickled)
def batch_align(db_p, ngs_dir, sesh):

    ngs_lst = Util.get_file_list(ngs_dir, ['.fasta', '.fastq'])
    ngs_count = len(ngs_lst)
    global_time = 0
    pickle_files = OrderedDict()

    base_output_dir = f'output/{sesh}'
    Util.create_dir(base_output_dir)
    Util.config_sesh_dir(base_output_dir)

    for i in range(len(ngs_lst)):

        ngs_p = ngs_lst[i]
        ngs_base = os.path.basename(ngs_lst[i])
        db_base = os.path.basename(db_p)

        print(f'\nAligning {db_base} and {ngs_base}...')
        total_time, res_list = align(db_p, ngs_p, sesh)
        print(f'Alignment ({i+1}/{ngs_count}) complete.')

        global_time += total_time
    
        output_fname = res_list[7]
        pickle_file = f'output/{sesh}/data/{output_fname}.pickle'

        # save all results list to use later rather than holding in memory
        with open(pickle_file, 'wb') as f:
            pickle.dump(res_list, f)

        Util.compress_pickle(f'output/{sesh}/data/{output_fname}.zip', pickle_file, output_fname)
        pickle_files[i+1] = [ngs_base, f'output/{sesh}/data/{output_fname}.zip', output_fname]
        
    global_time_str = Util.format_time(global_time)
    print(f'\nAlignment session completed in {global_time_str}.\nAll resulting outputs have been exported to output/{sesh}/results')

    while True:
        report_choice_init = input('Would you like to generate PDF reports for this session? (y/n) ')
        if report_choice_init in {'y', 'Y', 'n', 'N'}:
            break
        else:
            print('Invalid Input. Enter y or n (yes or no):')
    
    # user does not want to generate reports
    if report_choice_init in {'n', 'N'}:
        print(f'\nThank you for using FragmentFinder!\n\n')
        shutil.rmtree(f'output/{sesh}/temp')
        shutil.rmtree(f'output/{sesh}/reports')
        return 0
    
    # user wants to generate reports
    if report_choice_init in {'y', 'Y'}:
        print('\nBelow are the NGS files that were used during this session:')

        for key, val in pickle_files.items():
            print(f'{key}.) {val[0]}')

        report_choice_final = input('\nEnter the numbers corresponding to the alignment\noutput you would like to generate reports for\nseparated by commas (i.e. 1,2,3,...) \nor enter 0 to generate reports for all files:\n')
        
        if report_choice_final == '0':
            num_results = len(pickle_files)
            i = 1

            for key, val in pickle_files.items():
                pickle_path = val[1]
                out_f = val[2]
                Util.unzip_pickle(pickle_path, f'output/{sesh}/data')
                with open(f'output/{sesh}/data/{out_f}.pickle', 'rb') as file:
                    res_lst = pickle.load(file)
                generate_selection_report(res_lst)
                os.remove(f'output/{sesh}/data/{out_f}.pickle')
                print(f'\nPDF report ({i}/{num_results}) complete')
                i += 1

        else:
            input_list = report_choice_final.split(',')
            input_list = [int(value) for value in input_list]
            lst_len = len(input_list)
            for i in range(lst_len):
                val = pickle_files[input_list[i]]
                pickle_path = val[1]
                out_f = val[2]
                Util.unzip_pickle(pickle_path, f'output/{sesh}/data')
                with open(f'output/{sesh}/data/{out_f}.pickle', 'rb') as file:
                    res_lst = pickle.load(file)
                generate_selection_report(res_lst)
                os.remove(f'output/{sesh}/data/{out_f}.pickle')
                print(f'\nPDF report ({i+1}/{lst_len}) complete')

    print('\nAll PDF reports have been generated.\nThank you for using FragmentFinder!\n')
    shutil.rmtree(f'output/{sesh}/temp')
    return 0
        

def main():

    print('\nFragmentFinder Command Line Interface\n')

    while True:
        sesh_type = input("Enter 'n' to start a new session.\nEnter 'r' to generate reports from a previous session.\n")
        if sesh_type in {'n', 'N', 'r', 'R'}:
            break
        else:
            print('\nInvalid selection.\n')
            continue
    
    if sesh_type in {'n', 'N'}:
        session_num = input('\nPlease enter a unique Session ID.\nThis will be used to save results from this session.\n')
        print()
        db_path = get_db_path()

        while True:
            choice = input('\nTo align one NGS file, enter a single file path.\nTo align multiple NGS files, enter the path to a \ndirectory containing all of the NGS files.\n\n')
            if not Util.is_directory(choice) and Util.is_valid_file_format(choice, ['.fasta', 'fastq']):
                ngs_path = choice
                align_single_ngs(db_path, ngs_path, session_num)
                return 0
            elif not Util.is_directory(choice) and not Util.is_valid_file_format(choice, ['.fasta', '.fastq']):
                print(f'\nInvalid NGS file: {choice}.\nPlease select a valid FASTA or FASTQ type file.\n')
                continue
            if Util.is_directory(choice):
                ngs_lst = Util.get_file_list(choice, ['.fasta', '.fastq'])
                if len(ngs_lst) != 0:
                    print('\nThe following NGS files will be aligned to the Database file:\n')
                    for i in range(len(ngs_lst)):
                        f = os.path.basename(ngs_lst[i])
                        print(f'{i+1}.) {f}')
                    cont = input('\nWould you like to continue (y/n): ')
                    if cont in {'y', 'Y'}:
                        batch_align(db_path, choice, session_num)
                        return 0
                    else:
                        print("\nSelection: 'n'... Starting Over\n")
                        continue
                else:
                    print(f"Invalid directory: {choice}.\nPlease select a valid directory containing \nat least one FASTA or FASTQ type file.\n")
                    continue
        
    elif sesh_type in {'r', 'R'}:
        while True:
            sesh_id = input('Please enter ID of the Session you would like to generate reports for.\n')
            if not Util.is_directory(f'output/{sesh_id}'):
                print(f'Invalid Session ID: {sesh_id}.\nPlease enter a valid Session ID.\n')
                continue
            else:
                d_path = f'output/{sesh_id}/data'
                df_map = Util.make_file_map(d_path)
                print(f'\nBelow are the NGS files used during session {sesh_id}:')
                for key, val in df_map.items():
                    print(f'{key}.) {val[0]}')

                report_path = f'output/{sesh_id}/reports'
                os.makedirs(report_path)
                os.mkdir(f'output/{sesh_id}/temp')

                report_selection = input('\nEnter the numbers corresponding to the alignment\noutput you would like to generate reports \nfor separated by commas (i.e. 1,2,3,...) \nor enter 0 to generate reports for all files:\n')
                print()
                if report_selection == '0':
                    num_results = len(df_map)
                    i = 1
                    for key, val in df_map.items():
                        pickle_path = val[1]
                        out_f = val[2]
                        Util.unzip_pickle(pickle_path, f'output/{sesh_id}/data')
                        with open(f'output/{sesh_id}/data/{out_f}.pickle', 'rb') as file:
                            res_lst = pickle.load(file)
                        generate_selection_report(res_lst)
                        os.remove(f'output/{sesh_id}/data/{out_f}.pickle')
                        print(f'PDF report ({i}/{num_results}) complete\n')
                        i += 1

                else:
                    input_list = report_selection.split(',')
                    input_list = [int(value) for value in input_list]
                    lst_len = len(input_list)
                    for i in range(lst_len):
                        val = df_map[input_list[i]]
                        pickle_path = val[1]
                        out_f = val[2]
                        Util.unzip_pickle(pickle_path, f'output/{sesh_id}/data')
                        with open(f'output/{sesh_id}/data/{out_f}.pickle', 'rb') as file:
                            res_lst = pickle.load(file)
                        generate_selection_report(res_lst)
                        os.remove(f'output/{sesh_id}/data/{out_f}.pickle')
                        print(f'PDF report ({i+1}/{lst_len}) complete\n')

                shutil.rmtree(f'output/{sesh_id}/temp')
                print(f"\nAll PDF reports have been successfully generated.\nThey are available in {report_path}")
                return 0
                    


if __name__ == "__main__":
    main()
    


    
    






