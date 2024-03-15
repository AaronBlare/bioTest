import pandas as pd

data_path = 'E:/YandexDisk/pydnameth/draft/10_MetaEPIClock/MetaEpiAge/GPL13534/GSE75704/'
original_pheno_data = f'{data_path}pheno.csv'
supplementary_data = f'{data_path}41467_2018_5325_MOESM4_ESM.xlsx'

orig = pd.read_csv(original_pheno_data, index_col='subject.id')
supp = pd.read_excel(supplementary_data)
supp['sample_ID'] = supp['sample_ID'].replace({'Cntr':'B'}, regex=True)
supp.set_index('sample_ID', inplace=True)
supp['Age at death (years)'] = supp['Age at death (years)'].round()
supp['Age at death (years)'] = supp['Age at death (years)'].astype('Int64')

res = supp['Age at death (years)'].isin(orig['Age'])
ololo=0