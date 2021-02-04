from tkinter import *
from PIL import ImageTk, Image
import os
import ksea_fun as ksea
from tkinter import filedialog
from ttkthemes import themed_tk as tk
from tkinter import ttk
import pandas as pd

root = tk.ThemedTk()
root.get_themes()
root.set_theme("clam")
root.title("KSEA-Professional.exe")
root.iconbitmap('KSEA_icon.ico')
root.geometry('860x700')
root.resizable(0,0)

''' input values'''
ppIndex_path = ''
Experiment_title = ''
columns = ''
input_type = ''
scale = False
databases = ''
output_dir = ''
analysis_databases = {}

k_pic = ImageTk.PhotoImage(Image.open('KSEA_background_fin.png'))

k_title = ttk.Label(image=k_pic)

''' create a frame for all of the variable names. have a tick if you want them all to 
be used and an un tick then have the sample name and condition condition will be disabled
 if it is set to scale not comparison'''

'''function for getting database names'''


def remove_extensions(x):
    return [os.path.splitext(i)[0] for i in x]


''' function for updating the checkbox lists'''


def get_dbs(data_type):
    global analysis_databases
    global database_frame
    global input_type
    database_frame.grid_forget()

    if data_type == "select data type":
        return
    if data_type == "Proteomic" or data_type == "RNA-Seq":
        database_dir = "data_files/Gene_based"
    elif data_type == "Phospho-proteomic":
        database_dir = "data_files/Peptide_based"

    database_frame = VerticalScrolledFrame(root)
    database_frame.grid(row=4, column=4, rowspan=3)
    input_type = data_type
    databases = os.listdir(database_dir)
    analysis_databases = {}
    databases = remove_extensions(databases)
    '''for loop in check boxes'''
    for database in databases:
        analysis_databases[database] = BooleanVar()
        analysis_databases[database].set(0)
        check = ttk.Checkbutton(database_frame.interior,
                            text=database,
                            variable=analysis_databases[database],
                            onvalue=True,
                            offvalue=False,
                            width=20)
        check.pack(side=TOP, anchor=NW, pady=5)

    return


''' function for getting databases and columns for ksea analysis'''


def Get_info():
    out_col = []
    columns = list(ksea_columns.keys())
    for i in columns:
        x = ksea_columns[i]
        col = x.get()
        if col:
         out_col.append(str(i))
    out_dat = []

    databases = list(analysis_databases.keys())
    for i in databases:
        x = analysis_databases[i]
        db = x.get()
        if db:
         out_dat.append(str(i))

    out_dat = list(out_dat)
    out_col = list(out_col)
    if len(out_col) == len(columns):
        out_col = ':'

    return out_dat, out_col


class VerticalScrolledFrame(Frame):
    """A pure Tkinter scrollable frame that actually works!
    * Use the 'interior' attribute to place widgets inside the scrollable frame
    * Construct and pack/place/grid normally
    * This frame only allows vertical scrolling

    """
    def __init__(self, parent, *args, **kw):
        Frame.__init__(self, parent, *args, **kw)

        # create a canvas object and a vertical scrollbar for scrolling it
        vscrollbar = Scrollbar(self, orient=VERTICAL)
        vscrollbar.pack(fill=Y, side=RIGHT, expand=FALSE)
        canvas = Canvas(self, bd=0, highlightthickness=0,
                        yscrollcommand=vscrollbar.set)
        canvas.pack(side=LEFT, fill=BOTH, expand=TRUE)
        vscrollbar.config(command=canvas.yview)

        # reset the view
        canvas.xview_moveto(0)
        canvas.yview_moveto(0)

        # create a frame inside the canvas which will be scrolled with it
        self.interior = interior = Frame(canvas)
        interior_id = canvas.create_window(0, 0, window=interior,
                                           anchor=NW)

        # track changes to the canvas and frame width and sync them,
        # also updating the scrollbar
        def _configure_interior(event):
            # update the scrollbars to match the size of the inner frame
            size = (interior.winfo_reqwidth(), interior.winfo_reqheight())
            canvas.config(scrollregion="0 0 %s %s" % size)
            if interior.winfo_reqwidth() != canvas.winfo_width():
                # update the canvas's width to fit the inner frame
                canvas.config(width=interior.winfo_reqwidth())
        interior.bind('<Configure>', _configure_interior)

        def _configure_canvas(event):
            if interior.winfo_reqwidth() != canvas.winfo_width():
                # update the inner frame's width to fill the canvas
                canvas.itemconfigure(interior_id, width=canvas.winfo_width())
        canvas.bind('<Configure>', _configure_canvas)


def get_file():
    global ksea_columns
    global column_frame
    global ppIndex_path
    column_frame.grid_forget()
    file = filedialog.askopenfilename(initialdir='~/Documents',
                               title='Select an icon',
                               filetypes=(('csv files', '*.csv'), ('xlsx files', '*.xlsx'), ("all files", '*.*')))
    ppIndex_path = file
    file_type = os.path.splitext(file)[1]
    if file_type == ".csv":
        dat_columns = pd.read_csv(file, index_col=0).columns.values
    elif file_type == ".xlsx":
        dat_columns = pd.read_excel(file, sheet_name='ppIndex', index_col=0).columns.values


    ksea_columns = {}
    column_frame = VerticalScrolledFrame(root)
    column_frame.grid(row=3, column=5, rowspan=4)
    for column in dat_columns:
        ksea_columns[column] = BooleanVar()
        ksea_columns[column].set(True)
        check2 = ttk.Checkbutton(column_frame.interior,
                             text=column,
                             variable= ksea_columns[column],
                             onvalue=True,
                             offvalue=False,
                             width=20)
        check2.pack(anchor=E)

    File_path_entry.delete(0, END)
    File_path_entry.insert(0, file)


def get_dir():
    global output_dir
    output_dir = filedialog.askdirectory(initialdir='~/Documents', title='Select an output directory')
    Output_dir_entry.delete(0, END)
    Output_dir_entry.insert(0, output_dir)


def AnalyseKsea():
    output_databases, output_columns = Get_info()
    Experiment_title = str(Exp_title_entry.get())
    print(output_databases)
    if input_type == "Phospho-proteomic":
        analysis = ksea.KseaAnalysis(ppIndex_path=ppIndex_path,
                                     Experiment_title=Experiment_title,
                                     columns=output_columns,
                                     input_type=input_type,
                                     scale=scale,
                                     databases=output_databases,
                                     output_dir=output_dir
                                     )
    else:
        analysis = ksea.KseaAnalysisGene(ppIndex_path=ppIndex_path,
                                         Experiment_title=Experiment_title,
                                         columns=output_columns,
                                         input_type=input_type,
                                         scale=scale,
                                         databases=output_databases,
                                         output_dir=output_dir
                                         )

    analysis.start_progress(root=Progress_frame, max=800)

    analysis.progbar.grid(row=8, column=1, columnspan=5)

    analysis.run_ksea_all()


'''entries'''
Exp_title_entry = Entry(root, text='Experiment Title', width=50)
File_path_entry = Entry(root, text='File Path', width=30)
Output_dir_entry = Entry(root, text='Output Path', width=30)

'''progress bar pane'''
'''scale checkbox'''
Scale_Checkbox = ttk.Checkbutton(root, text= "scale", variable = scale, onvalue=True, offvalue=False)

'''data type drop down menu'''
dt_click = StringVar()
data_types = ['select data type', 'Proteomic', 'Phospho-proteomic', 'RNA-Seq']
dt_click.set(data_types[0])
data_type_drop = ttk.OptionMenu(root, dt_click, *data_types, command=get_dbs)

'''Find file'''
find_button = ttk.Button(root, text='Find File', command=get_file, width=15)
find_output = ttk.Button(root, text='Find Folder', command=get_dir, width=15)

'''run KSEA ttk.Button'''
run_button = ttk.Button(text="Run KSEA", command=AnalyseKsea, width=40)

'''labels'''
file_label = ttk.Label(root, text='Input File', width=12)
output_label = ttk.Label(root, text='Output Folder', width=12)
title_label = ttk.Label(root, text='Experiment Title', width=15)
database_label = ttk.Label(root, text='Databases', width=25)
Columns_label = ttk.Label(root, text='Columns', width=25)

'''Frames'''
column_frame = VerticalScrolledFrame(root, width=25)
database_frame = VerticalScrolledFrame(root, width=25)
Progress_frame = Frame(root)
column_frame.pack_propagate(0)
database_frame.pack_propagate(0)

'''padding canvases'''
top_canv = Canvas(root, bd=0, height=10)
frame_canv = Canvas(root, bd=0, height=40)
left_canv = Canvas(root, bd=0, width=30)
right_canv = Canvas(root,bd=0, width=30)

'''labels'''
k_title.grid(row=0, column=1, columnspan=5)
title_label.grid(row=2, column=1)
file_label.grid(row=3, column=1)
output_label.grid(row=4, column=1)
database_label.grid(row=2, column=4)
Columns_label.grid(row=2, column=5)
'''frames'''
column_frame.grid(row=3, column=5, rowspan=4)
database_frame.grid(row=4, column=4, rowspan=3)
Progress_frame.grid(row=8, column=1, columnspan=5, pady=15)

'''Entries'''
Exp_title_entry.grid(row=2, column=2, columnspan=2, pady=30)
File_path_entry.grid(row=3, column=2)
Output_dir_entry.grid(row=4, column=2)

'''checkbox'''
Scale_Checkbox.grid(row=5, column=1)

'''Buttons'''
run_button.grid(row=5, column=2, columnspan=2, pady=10)
find_button.grid(row=3, column=3, pady=10)
find_output.grid(row=4, column=3, pady=10)

'''Dropdown'''
data_type_drop.grid(row=3, column=4)

'''padding canvases'''
top_canv.grid(row=1, column=1, columnspan=5)
frame_canv.grid(row=6, column=1, columnspan=3)
left_canv.grid(row=0, column=0, rowspan=8)
right_canv.grid(row=0, column=6, rowspan=8)



root.mainloop()