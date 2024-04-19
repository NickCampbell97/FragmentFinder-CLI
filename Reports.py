import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import Util
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
import numpy as np

# create instance of report lab canvas
def create_canvas(db_path, ngs_path, session):
    db_exp = Util.fname_from_path(db_path)
    ng_exp = Util.fname_from_path(ngs_path)
    output_pdf = f'output/{session}/reports/{ng_exp}({db_exp}).pdf'
    c = canvas.Canvas(output_pdf, pagesize=letter)
    return c

# set style of title page
def style_title_page(c, page_width, page_height):
    c.setFillColor(HexColor('#FFFFFF'))
    c.rect(0,0,page_width,page_height,fill=1)
    c.setFillColor(HexColor('#071330'))
    c.rect(0, 710, page_width, page_height, fill=1)

# set style of the rest of pages
def style_peak_page(c, page_width, page_height):
    c.setFillColor(HexColor('#FFFFFF'))
    c.rect(0,0,page_width,page_height,fill=1)

# save the pdf report and delete all of the temp image files
def cleanup(pdf_canvas, session):
    pdf_canvas.save()
    Util.delete_temp_jpgs(f'output/{session}/temp')

# create title page
def write_title_page(c, reads, db_path, ngs_path, time_str, page_width, page_height, session):
    db_name = Util.fname_from_path(db_path)
    db_size = Util.f_size(db_path)
    ng_name = Util.fname_from_path(ngs_path)
    ng_size = Util.f_size(ngs_path)

    r_str = Util.format_large_ints(reads)
    date = Util.format_date()
    half = page_width/2.0
    sesh_str = f'Session ID: {session}'
    title = 'FragmentFinder Session Results'
    dbf = f'Database File  |  {db_name}  |  {db_size}'
    ngsf = f'NGS File  |  {ng_name}  |  {ng_size}'
    total_reads = f'Total Reads in NGS File  |  {r_str}'
    time = f'Time Taken (min:sec)  |  {time_str}'
    prompt = 'Below are all derived RNA fragments that have a minimum Read per Million value of 5.'

    title_width = c.stringWidth(title, 'Helvetica-Bold', 25)
    t_start = ((page_width / 2.0) - (title_width / 2.0)) 
    t_end = ((page_width / 2.0) + (title_width / 2.0))

    text_obj = c.beginText((half-t_start), (page_height-25))
    text_obj.setFillColor(HexColor('#071330'))
    text_obj.setFont('Helvetica-Bold', 25)
    c.drawText(text_obj)

    c.setFont('Helvetica-Bold', 25)
    c.setFillColor(HexColor('#faf0e6'))
    c.drawCentredString((half), (page_height-50), title)
    c.setStrokeColor(HexColor('#071330'))
    c.setFillColor(HexColor('#071330'))
    c.setFont('Helvetica-Bold', 16)
    c.drawCentredString(half, 650, sesh_str)
    c.drawCentredString(half, 610, dbf) 
    c.drawCentredString(half, 570, ngsf)
    c.drawCentredString(half, 530, total_reads) 
    c.drawCentredString(half, 490, time) 
    c.drawCentredString(half, 450, date) 
    c.setFont('Helvetica', 13)
    c.drawCentredString(half, 410, prompt) 
    c.line(50, 390, (page_width-50), 390)
    c.showPage()

# write data for a single peak of a DEV
def write_pdf_peak(c, h1, h2, h3, h4, l1, l2, l3, peak_string):
    text_obj = c.beginText(50, h1)
    text_obj.setFillColor(HexColor('#071330'))
    text_obj.setFont('Helvetica-Bold', 11)
    c.drawText(text_obj)
    c.setStrokeColor(HexColor('#071330'))
    c.setFillColor(HexColor('#071330'))
    c.setFont('Helvetica-Bold', 12)
    c.drawString(50, h1, peak_string)
    c.setFont('Helvetica', 11)
    c.drawString(50, h2, l1)
    c.drawString(50, h3, l2)
    c.drawString(50, h4, l3)

# get width and height of pdf page
width, height = letter

# create title page and write all of the single-peak alignment vectors to pdf file - two per page 
def plot_single_peaks(pdf_canvas, res, original_record_dict, read_count, db_path, ngs_path, time_str, session):
    
    style_title_page(pdf_canvas, width, height)
    write_title_page(pdf_canvas, read_count, db_path, ngs_path, time_str, width, height, session)

    fig = plt.figure()

    i = 1
    header = 750
    im_y = 521
    h1 = 510 
    h2 = h1 - 15 
    h3 = h2-15 
    h4 = h3-15 
    h5 = h4 - 20 

    for uid, lst in res.items():
        vector = lst[0]
        temp = lst[1]

        standard_deviation = np.std(vector)

        ax = fig.add_subplot(1,1,1)
        z = vector.copy()
        v = Util.regen_vector(z, read_count)

        ax.plot(v, label='Expression')
        
        ax.set_xlabel('Position on miRNA')
        ax.set_ylabel('Total Reads')

        for peak in temp:

            ax.plot(peak, z[peak], 'ko', label='Peak Start')
            rpm = round(((vector[peak] / read_count) * 1000000.0), 2)
            peak_slice = Util.val_check(z, peak, (z[peak]-(standard_deviation*0.1)))
            peak_start = peak_slice[0]
            peak_end = (peak_slice[-1]) + 18
            z[peak_start:peak_end] = vector[peak]

            ax.vlines(x=[peak_start, peak_end], ymin=0, ymax=z[peak], color='g', linestyle='dashed')
            ax.fill_betweenx([0, z[peak]], peak_start, peak_end, color='g', alpha=0.2)

            id_str = Util.delimit_uid(uid, '|')
            line1 = f'Reads per Million: {(rpm)}'
            line2 = f'Peak Start: {(peak_start+1)}'
            line3 = f'Peak End: {(peak_end+1)}'
            peak_string = f'Peak: {(original_record_dict[uid][peak_start:peak_end])}'
            
            pdf_canvas.setFillColor(HexColor('#071330'))
            pdf_canvas.setFont('Helvetica-Bold', 13)
            pdf_canvas.drawCentredString((width/2.0), header, id_str)
            
            write_pdf_peak(pdf_canvas, h1, h2, h3, h4, line1, line2, line3, peak_string)

        ax.legend()
        fig.savefig(f'output/{session}/temp/{uid}.jpg', dpi=72)
        fig.clf()
        pdf_canvas.drawImage(f'output/{session}/temp/{uid}.jpg', 145, im_y, width=325, height=225)
        pdf_canvas.line(50, h5, (width-50), h5)

        if i % 2 == 1:
            header = h5-30
            im_y = header - 229
            h1 = header - 250
            h2 = h1 - 15
            h3 = h2 - 15 
            h4 = h3 - 15
        elif i % 2 == 0:
            header = 750
            im_y = 521
            h1 = 510 
            h2 = h1 - 15 
            h3 = h2-15 
            h4 = h3-15 
            h5 = h4 - 10 
            pdf_canvas.showPage() 
        i += 1
    pdf_canvas.showPage()

# same things as plot_single_peaks but no title page created and writing multiple
# peaks for each uid/alignment vector combo
def plot_mult_peaks(pdf_canvas, res, original_record_dict, read_count, session):
    fig = plt.figure()

    for uid, lst in res.items():
        vector = lst[0]
        temp = lst[1]
        standard_deviation = np.std(vector)

        ax = fig.add_subplot(1,1,1)
        x = vector.copy()
        v = Util.regen_vector(x, read_count)

        ax.plot(v, label='Expression')
        
        ax.set_xlabel('Position on miRNA')
        ax.set_ylabel('Total Reads')

        z = vector.copy()

        i = 1
        h1 = 450
        h2 = 435
        h3 = 420
        h4 = 405
        h5 = 390

        for peak in temp:
            if i == 1:
                ax.plot(peak, z[peak], 'ko', label='Peak Start')
            if i > 1:
                ax.plot(peak, z[peak], 'ko')
            rpm = round(((vector[peak] / read_count) * 1000000.0), 2)
            peak_slice = Util.val_check(z, peak, (z[peak]-(standard_deviation*0.1)))
            peak_start = peak_slice[0]
            peak_end = (peak_slice[-1]) + 18
            z[peak_start:peak_end] = vector[peak]
            ax.vlines(x=[peak_start, peak_end], ymin=0, ymax=z[peak], color='g', linestyle='dashed')
            ax.fill_betweenx([0, z[peak]], peak_start, peak_end, color='g', alpha=0.2)

            id_str = Util.delimit_uid(uid, '|')
            line1 = f'Reads per Million: {(rpm)}'
            line2 = f'Peak Start: {(peak_start+1)}'
            line3 = f'Peak End: {(peak_end+1)}'
            peak_string = f'Peak: {(original_record_dict[uid][peak_start:peak_end])}'
            
            if i == 1:
                pdf_canvas.setFillColor(HexColor('#071330'))
                pdf_canvas.setFont('Helvetica-Bold', 13)
                pdf_canvas.drawCentredString((width/2.0), 750, id_str)
            
            write_pdf_peak(pdf_canvas, h1, h2, h3, h4, line1, line2, line3, peak_string)

            i += 1
            h1 = (h5-15)
            h2 = (h1-15)
            h3 = (h2-15)
            h4 = (h3-15)
            h5 = (h4-15)
        ax.legend()
        fig.savefig(f'output/{session}/temp/{uid}.jpg', dpi=72)
        fig.clf()
        pdf_canvas.drawImage(f'output/{session}/temp/{uid}.jpg', 113, 472, width=380, height=275)
        pdf_canvas.showPage() 
