import pandas as pd
import numpy as np
import pathlib as paths
from scipy.stats import ks_2samp
import os
import math
from tkinter import ttk

class KseaAnalysis:

    def __init__(self, ppIndex_path, Experiment_title, columns, scale, input_type, databases, output_dir):
        self.ppIndex = ppIndex_path
        self.Experiment_title = Experiment_title
        self.scale = scale
        self.output_dir = output_dir
        self.databases = databases
        self.input_type = input_type
        self.columns = columns
        self.var_dict = []
        self.ksea_output = {}
        self.database_dict = []
        self.dat_edges = []
        self.dataprog = 170
        self.progress = 0
        self.progbar =0
        self.root=0
        self.max_prog = 0
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        '''start progress'''

    def start_progress(self, root, max):
        self.root = root
        self.progbar = ttk.Progressbar(root, orient='horizontal', length=max, max=max, mode='determinate')
        self.update_progress(0)
        self.max_prog = max

    def update_progress(self, value):
        self.progress += value
        self.progbar['value'] = float(self.progress)
        self.progbar.update()

    def import_ppindex(self):
        file = self.ppIndex
        file_type = os.path.splitext(file)[1]
        if file_type == ".csv":
            dat = pd.read_csv(file, index_col=0)
        elif file_type == ".xlsx":
            dat = pd.read_excel(file, sheet_name='ppIndex', index_col=0)
        else:
            return 'file type not compatable'

        '''If sites column exists make sites index and remove from DataFrame'''
        if 'sites' in dat.columns.values:
            dat.index(dat.loc[:, 'sites'])
            del dat.loc[:, 'sites']

        '''filter DataFrame by selected columns'''
        if self.columns != ":":
            dat = dat.loc[:, self.columns]

        '''Scale DataFrame if stated'''
        if self.scale:

            dat = dat.apply(lambda x: x-min(x)/max(x)-min(x),axis=0)
            '''Change back to dataframe after scaling'''


        ''' 
        Build dictionary of rownames and corresponding variables
        this will be used to merge db edges
        '''
        #var_dict = dict.fromkeys(dat.index.values)
        var_dict = {}

        '''split edges using ;'''
        for i in dat.index:
            variables = str(i).split(';')
            for ii in variables:
                if ii not in ['', ' ']:
                    var_dict[ii] = i


        print(dat.columns)
        '''return all of the databases created'''
        self.ppIndex = dat
        self.var_dict = var_dict
        self.update_progress(value=(0.1*self.max_prog))


    def database_to_path(self, database):
        if self.input_type == 'Phospho-proteomic':
            dat_dir = 'Peptide_based'

        elif self.input_type in ['Proteomic', 'RNA-Seq']:
            dat_dir = 'Gene_based'

        print(dat_dir)
        database_file = str(database + '.csv')
        database_path = paths.PurePath('data_files', dat_dir, database_file)
        return database_path

    def import_database(self, database):
        database_path = self.database_to_path(database)
        database = pd.read_csv(database_path, index_col=0)
        nodes = database.index.values
        edge_col = database.columns[1]
        database_dict = dict.fromkeys(nodes)
        '''get list of all rownames'''
        pp_vars = self.var_dict.keys()

        '''split list to edges'''
        for i in nodes:
            edges = database.loc[i][edge_col]
            edges = str(edges).split(';')
            filt_edges = []
            ''' filter edges and return variable as ppindex rowname'''
            for ii in edges:
                if ii != '' and ii in pp_vars and ii not in filt_edges:
                    edge = self.var_dict[ii]
                    filt_edges.append(edge)

            if len(filt_edges) >= 3:
                database_dict[i] = filt_edges
            else:
                database_dict.pop(i, 'none')

        '''update progress bar'''
        self.update_progress(value=(0.1*self.dataprog))
        ''' might want to change this into a binary file with nodes as columns and marker presence in rows'''
        self.database_dict = database_dict

    def ksea(self, database):
        nodes = self.database_dict.keys()
        max_nodes = len(nodes)
        ksea_z = []
        ksea_p = []
        ksea_m = []
        ksea_D = []
        '''
        generate background values for reference background values are for every peptide which isnt 0 in the column 
        '''
        samples = self.ppIndex.columns.values
        bg = dict.fromkeys(samples)
        bg_mean = dict.fromkeys(samples)
        bg_q3 = dict.fromkeys(samples)
        bg_med = dict.fromkeys(samples)

        self.update_progress(value=(0.15*self.dataprog))
        for i in samples:
            x = self.ppIndex.loc[:, i]
            x = np.array(x[x!=0])
            bg[i] = x
            bg_mean[i] = np.mean(x)
            bg_q3[i] = np.percentile(x, 75)
            bg_med[i] = np.median(x)

        nodes_prog = 0.75*self.dataprog
        for node in nodes:
            self.update_progress(value=((1/max_nodes)*(nodes_prog)))
            edges = self.database_dict[node]
            df_edges = self.ppIndex.loc[edges, :]
            edge_D = {}
            edge_z = {}
            edge_p = {}
            edge_m = {'node':node,
                      'm':len(edges)}
            for i in samples:
                x = df_edges.loc[:, i]
                x = np.array(x[x != 0])

                if len(x) >=3:
                    edge_sd = np.std(x)
                    edge_q3 = np.percentile(x, 75)
                    edge_med = np.median(x)
                    '''calculate zscore using values'''
                    edge_z[i] = (edge_med - bg_med[i])*(math.sqrt(len(x)) / edge_sd)
                    edge_D[i] = (edge_med - bg_med[i]) + (edge_q3 - bg_q3[i])
                    edge_p[i] = ks_2samp(data1=x,
                                         data2=bg[i],
                                         alternative='two-sided',
                                         mode='auto')[1]

            ksea_p.append(edge_p)
            ksea_D.append(edge_D)
            ksea_m.append(edge_m)
            ksea_z.append(edge_z)

        '''convert output lists to dataframes'''

        df_p = pd.DataFrame(ksea_p, index=nodes, columns=samples)
        df_d = pd.DataFrame(ksea_D, index=nodes, columns=samples)
        df_m = pd.DataFrame(ksea_m)
        df_z = pd.DataFrame(ksea_z, index=nodes, columns=samples)

        'make keys for output'
        output_keys = [str(database + "_zscore"),
            str(database + "_distance"),
            str(database + "_pvalue"),
            str(database + "_m")]

        self.ksea_output[output_keys[0]] = df_z
        self.ksea_output[output_keys[1]] = df_d
        self.ksea_output[output_keys[2]] = df_p
        self.ksea_output[output_keys[3]] = df_m

        print(self.ksea_output.keys())

    def ksea_to_xlsx(self):
        xls_path = str(self.output_dir+"/"+str(self.Experiment_title +'KSEA_output.xlsx'))
        writer = pd.ExcelWriter(xls_path, engine='xlsxwriter')
        sheets = self.ksea_output.keys()
        for i in sheets:
            self.update_progress(value=(0.25*self.max_prog))
            x = self.ksea_output[i]
            x.to_excel(writer, sheet_name=i)
        writer.save()
        print(str(xls_path)+" made")
        open(xls_path, "r")
        os.startfile(xls_path)

    def run_ksea(self, database):
        print('single database')
        self.import_ppindex()
        self.import_database(str(database))
        self.ksea(str(database))
        self.ksea_to_xlsx()


    def run_ksea_all(self):
        print('multiple database')
        self.import_ppindex()
        self.dataprog = (0.8*self.max_prog)/len(self.databases)
        for i in self.databases:
            self.import_database(database=str(i))
            self.ksea(database=str(i))
        self.ksea_to_xlsx()


class KseaAnalysisGene:


    def __init__(self, ppIndex_path, Experiment_title, columns, scale, input_type, databases, output_dir):
        self.ppIndex = ppIndex_path
        self.Experiment_title = Experiment_title
        self.scale = scale
        self.output_dir = output_dir
        self.databases = databases
        self.input_type = input_type
        self.columns = columns
        self.var_dict = []
        self.ksea_output = {}
        self.database_dict = []
        self.dat_edges = []
        self.dataprog = 170
        self.progress = 0
        self.progbar = 0
        self.root = 0
        self.max_prog = 0
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        '''start progress'''

    def start_progress(self, root, max):
        self.root = root
        self.progbar = ttk.Progressbar(root, orient='horizontal', length=max, max=max, mode='determinate')
        self.update_progress(0)
        self.max_prog = max

    def update_progress(self, value):
        self.progress += value
        self.progbar['value'] = float(self.progress)
        self.progbar.update()

    def import_ppindex(self):
        file = self.ppIndex
        file_type = os.path.splitext(file)[1]
        if file_type == ".csv":
            dat = pd.read_csv(file, index_col=0)
        elif file_type == ".xlsx":
            dat = pd.read_excel(file, sheet_name='ppIndex', index_col=0)
        else:
            return 'file type not compatable'

        '''If sites column exists make sites index and remove from DataFrame'''
        if 'sites' in dat.columns.values:
            dat.index(dat.loc[:, 'sites'])
            del dat.loc[:, 'sites']

        '''filter DataFrame by selected columns'''
        if self.columns != ":":
            dat = dat.loc[:, self.columns]

        '''Scale DataFrame if stated'''
        if self.scale:
            dat = dat.apply(lambda x: x-min(x)/max(x)-min(x),axis=0)

        self.ppIndex = dat
        self.var_dict = list(dat.index.values)
        self.update_progress(value=(0.1 * self.max_prog))

    def database_to_path(self, database):
        if self.input_type == 'Phospho-proteomic':
            dat_dir = 'Peptide_based'

        elif self.input_type in ['Proteomic', 'RNA-Seq']:
            dat_dir = 'Gene_based'

        print(dat_dir)
        database_file = str(database + '.csv')
        database_path = paths.PurePath('data_files', dat_dir, database_file)
        return database_path

    def import_database(self, database):
        database_path = self.database_to_path(database)
        database = pd.read_csv(database_path, index_col=False).dropna()
        nodes = database.iloc[:,  0]
        print(nodes)
        database_dict = dict.fromkeys(nodes)
        '''get list of all rownames'''
        pp_vars = self.var_dict

        '''split list to edges'''
        for i, node in enumerate(nodes):
            edges = database.iloc[i,2]
            edges = str(edges).split(';')
            filt_edges = []
            ''' filter edges and return variable as ppindex rowname'''
            for ii in edges:
                if ii != '' and ii in pp_vars and ii not in filt_edges:
                    filt_edges.append(ii)

            if len(filt_edges) >= 3:
                database_dict[node] = filt_edges
            else:
                database_dict.pop(node, 'none')

        '''update progress bar'''
        self.update_progress(value=(0.1 * self.dataprog))
        ''' might want to change this into a binary file with nodes as columns and marker presence in rows'''
        self.database_dict = database_dict

    def ksea(self, database):
        nodes = self.database_dict.keys()
        max_nodes = len(nodes)
        ksea_z = []
        ksea_p = []
        ksea_m = []
        ksea_D = []
        '''
        generate background values for reference background values are for every peptide which isnt 0 in the column 
        '''
        samples = self.ppIndex.columns.values
        bg = dict.fromkeys(samples)
        bg_mean = dict.fromkeys(samples)
        bg_q3 = dict.fromkeys(samples)
        bg_sd = dict.fromkeys(samples)
        bg_med = dict.fromkeys(samples)

        self.update_progress(value=(0.15 * self.dataprog))
        for i in samples:
            x = self.ppIndex.loc[:, i]
            x = np.array(x[x != 0])
            bg[i] = x
            bg_mean[i] = np.mean(x)
            bg_sd[i] = np.std(x)
            bg_q3[i] = np.percentile(x, 75)
            bg_med[i] = np.median(x)

        nodes_prog = 0.75 * self.dataprog
        for node in nodes:
            self.update_progress(value=((1 / max_nodes) * nodes_prog))
            edges = self.database_dict[node]
            df_edges = self.ppIndex.loc[edges, :]
            edge_D = {}
            edge_z = {}
            edge_p = {}
            edge_m = {'node': node,
                      'm': len(edges)}
            for i in samples:
                x = df_edges.loc[:, i]
                x = np.array(x[x != 0])

                if len(x) >= 3:
                    edge_sd = np.std(x)
                    edge_q3 = np.percentile(x, 75)
                    edge_med = np.median(x)
                    '''calculate zscore using values'''
                    edge_z[i] = (edge_med - bg_med[i])*(math.sqrt(len(x)) / edge_sd)
                    edge_D[i] = (edge_med - bg_med[i]) + (edge_q3- bg_q3[i])
                    edge_p[i] = ks_2samp(data1=x,
                                         data2=bg[i],
                                         alternative='two-sided',
                                         mode='auto')[1]

            ksea_p.append(edge_p)
            ksea_D.append(edge_D)
            ksea_m.append(edge_m)
            ksea_z.append(edge_z)

        '''convert output lists to dataframes'''

        df_p = pd.DataFrame(ksea_p, index=nodes, columns=samples)
        df_d = pd.DataFrame(ksea_D, index=nodes, columns=samples)
        df_m = pd.DataFrame(ksea_m)
        df_z = pd.DataFrame(ksea_z, index=nodes, columns=samples)

        'make keys for output'
        output_keys = [str(database + "_zscore"),
                       str(database + "_distance"),
                       str(database + "_pvalue"),
                       str(database + "_m")]

        self.ksea_output[output_keys[0]] = df_z
        self.ksea_output[output_keys[1]] = df_d
        self.ksea_output[output_keys[2]] = df_p
        self.ksea_output[output_keys[3]] = df_m

        print(self.ksea_output.keys())

    def ksea_to_xlsx(self):
        xls_path = str(self.output_dir + "/" + str(self.Experiment_title + 'KSEA_output.xlsx'))
        writer = pd.ExcelWriter(xls_path, engine='xlsxwriter')
        sheets = self.ksea_output.keys()
        for i in sheets:
            self.update_progress(value=(0.25 * self.max_prog))
            x = self.ksea_output[i]
            x.to_excel(writer, sheet_name=i)
        writer.save()
        print(str(xls_path) + " made")
        os.startfile(xls_path)


    def run_ksea(self, database):
        self.import_ppindex()
        self.import_database(str(database))
        self.ksea(str(database))
        self.ksea_to_xlsx()

    def run_ksea_all(self):
        print('multiple database')
        self.import_ppindex()
        self.dataprog = (0.8 * self.max_prog) / len(self.databases)
        for i in self.databases:
            self.import_database(database=str(i))
            self.ksea(database=str(i))
        self.ksea_to_xlsx()
