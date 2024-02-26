# Code for figures in Buitrago Silva et al. 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.formula.api import ols
from scipy.stats import mannwhitneyu
import plotly.graph_objects as go # Import the graphical object
from sklearn.cluster import KMeans


# This should be Table S2
working_dir = '/Users/stephansanders/Dropbox/Papers_author/Dina_SLC6A1_Aug2023/Supplemental_Tables/'
input_file_name = f'{working_dir}Table_S2_SLC6A1_Variants_Submit.xlsx'
excel_file = pd.ExcelFile(input_file_name)
df = excel_file.parse('Data', header=0, index_col=0)
df.head()

# This should be Table S7
input_file_name2 = f'{working_dir}Table_S7_All_SLC6A1_SNVs.xlsx'
excel_file2 = pd.ExcelFile(input_file_name2)
df2 = excel_file2.parse('Data', header=0, index_col=0)
df2.head()

# Columns
uptake_values = 'Combined_mean_reuptake_rescale'
uptake_categories = 'Combined_Category'
surface_expression = 'UCSF_Surface'


def is_float(element) -> bool:
    try:
        float(element)
        return True
    except ValueError:
        return False


def make_figure_1 (df, uptake_values):
	
	df_up = df[(df[uptake_values] != '.')]
	df_sort = df_up.sort_values(uptake_values, axis=0, ascending=True, inplace=False)

	fig = plt.figure(figsize=(6,35))
	ax = fig.add_subplot(111)
	plt.rcParams['pdf.use14corefonts'] = True  # PDF export makes text as text
	plt.rcParams['axes.unicode_minus'] = False  # Make sure the minus sign prints out correctly

	y_pos = np.arange(len(df_sort[uptake_values]))
	ax.barh(y_pos, df_sort[uptake_values], xerr=df_sort['Combined_mean_reuptake_sem'], align='center', capsize=2)
	ax.set_yticks(y_pos)
	ax.set_yticklabels(df_sort['ProteinFull'].values.tolist())
	ax.invert_yaxis() 
	plt.axvline(0, linestyle='--')
	plt.axvline(0.524693499, linestyle='--')
	plt.axvline(-0.492047675, linestyle='--')
	plt.axvline(-0.820503901, linestyle='--')

	out_name = f'{working_dir}Buitrago_Silva_Fig1.pdf'  # Name of image
	plt.savefig(out_name)


def make_figure_2a (df, uptake_values):

	df_fv = df[(df[uptake_categories] != '.') & (df['Variant_impact_cat'] != '.')]
	df_fv.loc[:, "Variant_impact_cat"]=df_fv["Variant_impact_cat"].apply(str)

	metric_name = [uptake_values]
	class_name = ["PTV", "Missense/IF", "Missense/IF_CT", "Synonymous", ]

	list_of_scores = []
	list_of_classes = []
	list_of_groups = []

	for metric in metric_name:

	    # Put none null scores into a list with a matching list of gene names
	    scores = df_fv[metric][df_fv[metric] != '.'].values.tolist()
	    classes = df_fv['Variant_impact_cat'][df_fv[metric] != '.'].values.tolist()
	    groups = df_fv['Combined_Category'][df_fv[metric] != '.'].values.tolist()
	    for i, score in enumerate(scores):
	        if is_float(score):
	            list_of_scores.append(float(score))
	            list_of_classes.append(classes[i])
	            list_of_groups.append(groups[i])
	        
	# Make a swarm plot
	dict_metric = {'Class':list_of_classes,'Score':list_of_scores,'Group':list_of_groups}
	new_df = pd.DataFrame(dict_metric)

	fig, axes = plt.subplots()
	plt.rcParams['pdf.use14corefonts'] = True  # PDF export makes text as text
	plt.rcParams['axes.unicode_minus'] = False  # Make sure the minus sign prints out correctly
	plt.rcParams['figure.figsize']=(6,6)
	plt.rcParams.update(_get_seaborn_axes_style())
	image = sns.swarmplot(x='Class', y='Score', hue='Group', data=new_df, ax=axes, size=3,
	                      order=class_name)
	image = sns.boxplot(x='Class', y='Score', data=new_df, ax=axes, order=class_name)
	image.axhline(0, linestyle='--')
	image.axhline(0.4753065014345206, linestyle='--')
	image.axhline(-0.4920476745946898, linestyle='--')
	image.axhline(-0.8205039012005264, linestyle='--')
	image.set_ylim(-1.05,0.60)
	axes.set_xlabel('Class')
	axes.set_ylabel('Variant type')
	plt.legend([],[], frameon = False)

	out_name = f'{working_dir}Buitrago_Silva_Fig2a.pdf'  # Name of image
	plt.savefig(out_name)
    

def make_figure_2b (df, uptake_values):
	
	df_fg = df[(df[uptake_values] != '.') & (df['gnomADv2_Count_Cat_non_neuro'] != '.')]
	df_fg.loc[:, "gnomADv2_Count_Cat_non_neuro"]=df_fg["gnomADv2_Count_Cat_non_neuro"].apply(str)

	metric_name = [uptake_values]
	class_name = ["0", "1", "2 to 10", "≥10"]

	list_of_scores = []
	list_of_classes = []
	list_of_groups = []

	for metric in metric_name:

	    # Put none null scores into a list with a matching list of gene names
	    scores = df_fg[metric][df_fg[metric] != '.'].values.tolist()
	    classes = df_fg['gnomADv2_Count_Cat_non_neuro'][df_fg[metric] != '.'].values.tolist()
	    groups = df_fg['Combined_Category'][df_fg[metric] != '.'].values.tolist()
	    for i, score in enumerate(scores):
	        if is_float(score):
	            list_of_scores.append(float(score))
	            list_of_classes.append(classes[i])
	            list_of_groups.append(groups[i])
	        
	# Make a swarm plot
	dict_metric = {'Class':list_of_classes,'Score':list_of_scores,'Group':list_of_groups}
	new_df = pd.DataFrame(dict_metric)

	fig, axes = plt.subplots()
	plt.rcParams['pdf.use14corefonts'] = True  # PDF export makes text as text
	plt.rcParams['axes.unicode_minus'] = False  # Make sure the minus sign prints out correctly
	plt.rcParams['figure.figsize'] = (6,6)
	plt.rcParams.update(_get_seaborn_axes_style())
	image = sns.swarmplot(x='Class', y='Score', hue='Group', data=new_df, size=3, ax=axes, 
	                      order=class_name)
	image = sns.boxplot(x='Class', y='Score', data=new_df, ax=axes, order=class_name)
	image.axhline(0, linestyle='--')
	image.axhline(0.4753065014345206, linestyle='--')
	image.axhline(-0.4920476745946898, linestyle='--')
	image.axhline(-0.8205039012005264, linestyle='--')
	image.set_ylim(-1.05,0.60)
	axes.set_xlabel('Class')
	axes.set_ylabel('gnomAD Allele Count')
	plt.legend([],[], frameon = False)

	out_name = f'{working_dir}Buitrago_Silva_Fig2b.pdf'  # Name of image
	plt.savefig(out_name)


def make_figure_2c (df, uptake_values):

	df_fc = df[(df[uptake_categories] != '.') & (df['ClinVar_call_simp'] != '.')]

	metric_name = [uptake_values]
	class_name = ["Pathogenic", "Likely pathogenic", 
	              "Uncertain significance", "Likely benign", 
	              "Benign", "Conflicting interpretations", 
	              "No score"]

	list_of_scores = []
	list_of_classes = []
	list_of_groups = []

	for metric in metric_name:

	    # Put none null scores into a list with a matching list of gene names
	    scores = df[metric][df[metric] != '.'].values.tolist()
	    classes = df['ClinVar_call_simp'][df[metric] != '.'].values.tolist()
	    groups = df['Combined_Category'][df[metric] != '.'].values.tolist()
	    for i, score in enumerate(scores):
	        if is_float(score):
	            list_of_scores.append(float(score))
	            list_of_classes.append(classes[i])
	            list_of_groups.append(groups[i])
	        
	# Make a swarm plot
	dict_metric = {'Class':list_of_classes,'Score':list_of_scores,'Group':list_of_groups}
	new_df = pd.DataFrame(dict_metric)

	fig, axes = plt.subplots()
	plt.rcParams['pdf.use14corefonts'] = True  # PDF export makes text as text
	plt.rcParams['axes.unicode_minus'] = False  # Make sure the minus sign prints out correctly
	plt.rcParams['figure.figsize'] = (10,6)
	plt.rcParams.update(_get_seaborn_axes_style())
	image = sns.swarmplot(x='Class', y='Score', hue='Group', data=new_df, ax=axes, size=3,
	                      order=class_name)
	image = sns.boxplot(x='Class', y='Score', data=new_df, ax=axes, order=class_name)
	image.axhline(0, linestyle='--')
	image.axhline(0.4753065014345206, linestyle='--')
	image.axhline(-0.4920476745946898, linestyle='--')
	image.axhline(-0.8205039012005264, linestyle='--')
	image.set_ylim(-1.05,0.60)
	axes.set_xlabel('Class')
	axes.set_ylabel('GABA Reuptake')
	plt.legend([],[], frameon = False)

	out_name = f'{working_dir}Buitrago_Silva_Fig2c.pdf'  # Name of image
	plt.savefig(out_name)
    

def make_figure_2d (df, uptake_categories):

	# Numbers for sankey plot
	df_fc = df[(df[uptake_categories] != '.')]

	# Level one
	all_var = df_fc.shape[0]
	all_var_to_path = df_fc[(df_fc["ClinVar_call"].str.contains("athogenic")) & 
	                        (~df_fc["ClinVar_call"].str.contains("Conflicting"))].shape[0]
	all_var_to_benign = df_fc[(df_fc["ClinVar_call"].str.contains("enign"))].shape[0]
	all_var_to_vus = all_var - all_var_to_path - all_var_to_benign

	# Level two
	path_to_path2 = df_fc[(df_fc["ClinVar_call"].str.contains("athogenic")) & 
	                        (~df_fc["ClinVar_call"].str.contains("Conflicting")) &
	                        (df_fc[uptake_categories].str.contains("LoF"))].shape[0]
	path_to_benign2 = df_fc[(df_fc["ClinVar_call"].str.contains("athogenic")) & 
	                        (~df_fc["ClinVar_call"].str.contains("Conflicting")) &
	                        (~df_fc[uptake_categories].str.contains("LoF"))].shape[0]

	benign_to_path2 = df_fc[(df_fc["ClinVar_call"].str.contains("enign")) & 
	                          (df_fc[uptake_categories].str.contains("LoF"))].shape[0]
	benign_to_benign2 = df_fc[(df_fc["ClinVar_call"].str.contains("enign")) & 
	                          (~df_fc[uptake_categories].str.contains("LoF"))].shape[0]

	all_path2 = df_fc[(df_fc[uptake_categories].str.contains("LoF"))].shape[0]
	all_benign2 = df_fc[(~df_fc[uptake_categories].str.contains("LoF"))].shape[0]

	vus_to_path2 = all_path2 - path_to_path2 - benign_to_path2
	vus_to_benign2 = all_benign2 - benign_to_benign2 - path_to_benign2


	print(all_var_to_path, all_var_to_vus, all_var_to_benign, path_to_path2, path_to_benign2, \
		vus_to_path2, vus_to_benign2, benign_to_path2, benign_to_benign2)

	# make the Sankey plot

	node_label = ["all_var1", "all_var2", "Path", "VUS", "Benign", "Path2", "Benign2"]
	node_dict = {y:x for x, y in enumerate(node_label)}

	source = ['all_var1','all_var1','all_var2','all_var2','Path','Path','VUS','VUS','Benign','Benign']
	target = ['Path','VUS','VUS','Benign','Path2','Benign2','Path2','Benign2','Path2','Benign2'] 
	values = [ 32, 100, 63, 18, 30, 2, 94, 69, 3, 15 ]

	source_node = [node_dict[x] for x in source]
	target_node = [node_dict[x] for x in target]

	fig = go.Figure( 
	    data=[go.Sankey( 
	        # This part is for the node information
	        node = dict( 
	            label = node_label
	        ),
	        # This part is for the link information
	        link = dict(
	            source = source_node,
	            target = target_node,
	            value = values
	        ))])

	# With this save the plots 
	# plot(fig,
	#      image_filename='sankey_plot_2', 
	#      image='png', 
	#      image_width=1000, 
	#      image_height=600
	# )
	# And shows the plot
	fig.show()

	# In the browser click 'File', 'Export to PDF' and it is vectorised
	# Reaaranged the categories in Illustrator


def make_figure_2e (df, uptake_values):
	
	# Make dataframe for de novo by phenotype
	df_fd = df[(df[uptake_values] != '.') & (df['Inheritance_simp'] == 'De_novo')]

	df_fd_seizure = df_fd[(df_fd['Seizures'] == 'Yes') ]
	df_fd_autism = df_fd[(df_fd['Autism'] == 'Yes') ]
	df_fd_delay = df_fd[(df_fd['DD_or_ID'] == 'Yes') ]
	df_fd_sch = df_fd[(df_fd['Schizophrenia'] == 'Yes') & 
	                   (df_fd['ProteinFull'] != 'p.Ala305Thr') & 
	                   (df_fd['ProteinFull'] != 'p.Ala334Thr') ]

	df_fd_seizure.insert(0, 'Pheno', 'Seizures')
	df_fd_autism.insert(0, 'Pheno', 'ASD')
	df_fd_delay.insert(0, 'Pheno', 'NDD')
	df_fd_sch.insert(0, 'Pheno', 'Schizophrenia')

	df_pd = pd.concat([df_fd_seizure, df_fd_autism, df_fd_delay, df_fd_sch])
	df_pd.head()

	metric_name = [uptake_values]
	class_name = ["Seizures", "ASD", "NDD", "Schizophrenia"]
	list_of_scores = []
	list_of_classes = []
	list_of_groups = []

	for metric in metric_name:

	    # Put none null scores into a list with a matching list of gene names
	    scores = df_pd[metric][df_pd[metric] != '.'].values.tolist()
	    classes = df_pd['Pheno'][df_pd[metric] != '.'].values.tolist()
	    groups = df_pd['Combined_Category'][df_pd[metric] != '.'].values.tolist()
	    for i, score in enumerate(scores):
	        if is_float(score):
	            list_of_scores.append(float(score))
	            list_of_classes.append(classes[i])
	            list_of_groups.append(groups[i])
	        
	# Make a swarm plot
	dict_metric = {'Class':list_of_classes,'Score':list_of_scores,'Group':list_of_groups}
	new_df = pd.DataFrame(dict_metric)

	fig, axes = plt.subplots()
	plt.rcParams['pdf.use14corefonts'] = True  # PDF export makes text as text
	plt.rcParams['axes.unicode_minus'] = False  # Make sure the minus sign prints out correctly
	plt.rcParams['figure.figsize']=(6,6)
	plt.rcParams.update(_get_seaborn_axes_style())
	image = sns.swarmplot(x='Class', y='Score', hue='Group', data=new_df, ax=axes, size=3, 
	                      order=class_name)
	image = sns.boxplot(x='Class', y='Score', data=new_df, ax=axes, order=class_name)
	image.axhline(0, linestyle='--')
	image.axhline(0.4753065014345206, linestyle='--')
	image.axhline(-0.4920476745946898, linestyle='--')
	image.axhline(-0.8205039012005264, linestyle='--')
	image.set_ylim(-1.05,0.60)
	axes.set_xlabel('Class')
	axes.set_ylabel('Phenotype')
	plt.legend([],[], frameon = False)

	out_name = f'{working_dir}Buitrago_Silva_Fig2e.pdf'  # Name of image
	plt.savefig(out_name)

	# Test significance of schiozphrenia vs other
	sch_rees_uptake = df_pd[uptake_values][(df_pd['Pheno'] == 'Schizophrenia') & 
	                                       (df_pd['ProteinFull'] != 'p.Ala305Thr') & 
	                                       (df_pd['ProteinFull'] != 'p.Ala334Thr')]
	other_uptake = df_pd[uptake_values][(df_pd['Pheno'] != 'Schizophrenia')]
	stat, pval = mannwhitneyu(pd.to_numeric(sch_rees_uptake), pd.to_numeric(other_uptake), alternative='two-sided')
	print(f'P = {pval}')


def make_figure_2f (df, uptake_values):
	df_fum = df[(df[uptake_values] != '.') & \
               (df['Variant_impact_cat'] == 'Missense/IF') & \
               (df['Inheritance_simp'] == 'De_novo')].copy()
	df_fum.loc[:, "Unique_count_str"]=df_fum["Unique_count"].apply(str)

	metric_name = [uptake_values]
	class_name = ["1", "2", "3", "4", "5", "6", "10", "11"]

	list_of_scores = []
	list_of_classes = []
	list_of_groups = []

	for metric in metric_name:

	    # Put none null scores into a list with a matching list of gene names
	    scores = df_fum[metric][df_fum[metric] != '.'].values.tolist()
	    classes = df_fum['Unique_count_str'][df_fum[metric] != '.'].values.tolist()
	    groups = df_fum['Combined_Category'][df_fum[metric] != '.'].values.tolist()
	    for i, score in enumerate(scores):
	        if is_float(score):
	            list_of_scores.append(float(score))
	            list_of_classes.append(classes[i])
	            list_of_groups.append(groups[i])
	        
	# Make a swarm plot
	dict_metric = {'Class':list_of_classes,'Score':list_of_scores,'Group':list_of_groups}
	new_df = pd.DataFrame(dict_metric)

	fig, axes = plt.subplots()
	plt.rcParams['pdf.use14corefonts'] = True  # PDF export makes text as text
	plt.rcParams['axes.unicode_minus'] = False  # Make sure the minus sign prints out correctly
	plt.rcParams['figure.figsize']=(6,6)
	plt.rcParams.update(_get_seaborn_axes_style())
	image = sns.swarmplot(x='Class', y='Score', hue='Group', data=new_df, ax=axes, size=3,
	                      order=class_name)
	image = sns.boxplot(x='Class', y='Score', data=new_df, ax=axes, 
	                      order=class_name)
	image.axhline(0, linestyle='--')
	image.axhline(0.4753065014345206, linestyle='--')
	image.axhline(-0.4920476745946898, linestyle='--')
	image.axhline(-0.8205039012005264, linestyle='--')
	image.set_ylim(-1.05,0.60)
	axes.set_xlabel('Class')
	axes.set_ylabel('Number of individuals with recurrent variants')
	plt.legend([],[], frameon = False)

	out_name = f'{working_dir}Buitrago_Silva_Fig2f.pdf'  # Name of image
	plt.savefig(out_name)
    

def make_figure_3a (df, uptake_values, surface_expression):
	
	df_surf = df[(df[uptake_values] != '.') & (df[surface_expression] != '.')]
	df_surf = df_surf.reset_index(level=None, drop=True, inplace=False)
	df_a = df_surf[[uptake_values, surface_expression]]
	empty_vector = 0.468765210820772

	# K means
	kmeans = KMeans(3)
	kmeans.fit(df_a)
	identified_clusters = kmeans.fit_predict(df_a)
	df_surf['Combo_Clusters'] = identified_clusters 

	# Figure of clusters
	fig, ax = plt.subplots()
	plt.rcParams['pdf.use14corefonts'] = True  # PDF export makes text as text
	plt.rcParams['axes.unicode_minus'] = False  # Make sure the minus sign prints out correctly
	plt.rcParams['figure.figsize']=(8,8)
	plt.rcParams.update(_get_seaborn_axes_style())
	plt.scatter(df_surf[surface_expression],df_surf[uptake_values],c=df_surf['Combo_Clusters'],cmap='rainbow')
	plt.ylabel('GABA reuptake')
	plt.xlabel('Surface expression')

	# Add a legend
	scatter = ax.scatter(df_surf[surface_expression],df_surf[uptake_values],c=df_surf['Combo_Clusters'],cmap='rainbow')

	# produce a legend with the unique colors from the scatter
	legend1 = ax.legend(*scatter.legend_elements(),
	                    loc="upper left", title="Classes")
	ax.add_artist(legend1)

	# plt.vlines(0, min(df_surf[uptake_values]), max(df_surf[uptake_values]), linestyle='--')
	plt.vlines(1, min(df_surf[uptake_values]), max(df_surf[uptake_values]), linestyle='--')
	plt.vlines(empty_vector, min(df_surf[uptake_values]), max(df_surf[uptake_values]), linestyle='--')

	plt.hlines(-0.4920476745946898, min(df_surf[surface_expression]), max(df_surf[surface_expression]), linestyle='--')
	plt.hlines(0.4753065014345206, min(df_surf[surface_expression]), max(df_surf[surface_expression]), linestyle='--')
	plt.hlines(-0.8205039012005264, min(df_surf[surface_expression]), max(df_surf[surface_expression]), linestyle='--')

	def label_point(x, y, val, ax):
	    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
	    for i, point in a.iterrows():

	        if point['y'] > 0.2 and point['x'] < -0.1:
	            ax.text(point['x'], point['y'], str(point['val']))
	        elif point['y'] > 0.6 and point['x'] < 0.45:
	            ax.text(point['x'], point['y'], str(point['val']))
	        elif point['y'] > 0.2 and point['y'] < 0.6 and point['x'] > 0.5:
	            ax.text(point['x'], point['y'], str(point['val']))
	            

	label_point( df_surf[surface_expression], df_surf[uptake_values], df_surf['ProteinFull'], ax)

	out_name = f'{working_dir}Buitrago_Silva_Fig3a.pdf'  # Name of image
	plt.savefig(out_name)


def make_figure_4a (df, uptake_values):
	
	df_meta = df[(df[uptake_values] != '.') & (df["ClinPred_rankscore"] != '.')].copy()

	df_meta.loc[:, uptake_values]=df_meta[uptake_values].apply(float)
	df_meta.loc[:, "ClinPred_rankscore"]=df_meta["ClinPred_rankscore"].apply(float)

	fig, axes = plt.subplots()
	plt.rcParams['pdf.use14corefonts'] = True  # PDF export makes text as text
	plt.rcParams['axes.unicode_minus']=False  # Make sure the minus sign prints out correctly
	plt.rcParams['figure.figsize']=(10,10)
	plt.rcParams.update(_get_seaborn_axes_style())

	# Convert x and y columns to numeric data types
	df_meta['ClinPred_rankscore'] = pd.to_numeric(df_meta['ClinPred_rankscore'])
	df_meta[uptake_values] = pd.to_numeric(df_meta[uptake_values])

	image = sns.regplot(x='ClinPred_rankscore', y=uptake_values, data=df_meta, scatter_kws={'s':10}, line_kws={"color": "red"})
	image.axhline(0, linestyle='--')
	image.axhline(0.4753065014345206, linestyle='--')
	image.axhline(-0.4920476745946898, linestyle='--')
	image.axhline(-0.8205039012005264, linestyle='--')
	image.set_ylim(-1.05,0.60)


	axes.set_ylabel('GABA uptake')
	axes.set_xlabel('ClinPred Rankscore')

	for line in range(0,df_meta.shape[0]):
	     image.text(df_meta['ClinPred_rankscore'].iloc[line]+0.01, df_meta[uptake_values].iloc[line],
	     df_meta['ProteinFull'].iloc[line], horizontalalignment='left',
	     size='xx-small', color='black')

	out_name = f'{working_dir}Buitrago_Silva_Fig4a.pdf'  # Name of image
	plt.savefig(out_name)

	res = ols(f"{uptake_values} ~ ClinPred_rankscore", data=df_meta).fit()
	print(res.summary())
	print(res.params)


def make_figure_5 (df):

	fig = plt.figure(figsize=(10, 10))
	plt.rcParams['pdf.use14corefonts'] = True  # PDF export makes text as text
	plt.rcParams['axes.unicode_minus'] = False  # Make sure the minus sign prints out correctly
	plt.rcParams.update(_get_seaborn_axes_style())

	df_num = df[(df["Mutability"] != '.') & (df["Pred_uptake"] != '.')]
	df_num = df_num.astype({'Mutability':'float', 'Pred_uptake':'float'})
	df_num_zero = df_num[(df_num["Unique_count"] == 0)]
	df_num_one = df_num[(df_num["Unique_count"] == 1)]
	df_num_two = df_num[(df_num["Unique_count"] == 2)]
	df_num_three = df_num[(df_num["Unique_count"] == 3)]
	df_num_four = df_num[(df_num["Unique_count"] == 4)]
	df_num_five = df_num[(df_num["Unique_count"] >= 5)]

	ax1 = fig.add_subplot(111)
	ax1.axhline(0, linestyle='--')
	ax1.axhline(-0.4920476745946898, linestyle='--')
	ax1.axhline(-0.8205039012005264, linestyle='--')
	ax1.scatter(x = df_num_zero['Mutability'], y = df_num_zero['Pred_uptake'], s = 0.5, c = "lightgrey")
	ax1.scatter(x = df_num_one['Mutability'], y = df_num_one['Pred_uptake'], s = 5, c = "black")
	ax1.scatter(x = df_num_two['Mutability'], y = df_num_two['Pred_uptake'], s = 20, c = "blue")
	ax1.scatter(x = df_num_three['Mutability'], y = df_num_three['Pred_uptake'], s = 40, c = "purple")
	ax1.scatter(x = df_num_four['Mutability'], y = df_num_four['Pred_uptake'], s = 60, c = "red")
	ax1.scatter(x = df_num_five['Mutability'], y = df_num_five['Pred_uptake'], s = 200, c = "red")

	ax1.set_xscale('log')  # set_xscale is a function, not a string

	out_name = f'{working_dir}Buitrago_Silva_Fig5.pdf'  # Name of image
	plt.savefig(out_name)


def _get_seaborn_axes_style(style="ticks"):
	# Use the properties from `sns.set_style("ticks")`
	params = sns.axes_style(style)
	# Remove invalid matplotlib parameter (if present)
	params.pop("font.sans-serif", None)
	return params


make_figure_1 (df, uptake_values)
make_figure_2a (df, uptake_values)
make_figure_2b (df, uptake_values)
make_figure_2c (df, uptake_values)
make_figure_2d (df, uptake_categories)
make_figure_2e (df, uptake_values)
make_figure_2f (df, uptake_values)
make_figure_3a (df, uptake_values, surface_expression)
make_figure_4a (df, uptake_values)
make_figure_5 (df2)

