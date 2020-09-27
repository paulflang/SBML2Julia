import csv
import os
import pandas as pd
from natsort import natsorted, index_natsorted

basedir = os.path.dirname(os.path.abspath(__file__))

## Write condition table
out_path = os.path.join(basedir, 'condition_table.tsv')
with open(out_path, 'wt') as f:
    tsv_writer = csv.writer(f, delimiter='\t')
    tsv_writer.writerow(['conditionId', 'mu_a', 'mu_b', 'alpha_ab', 'alpha_ba', 'alpha_aa', 'alpha_bb'])

    # Write single species conditions
    for i in range(1, 13):
        simulationConditionId = f'c_{i}'
        mu_a = f'mu_{i}'
        mu_b = '0'
        alpha_ab = '0'
        alpha_ba = '0'
        alpha_aa = f'alpha_{i}_{i}'
        alpha_bb = '0'
        tsv_writer.writerow([simulationConditionId, mu_a, mu_b, alpha_ab, alpha_ba, alpha_aa, alpha_bb])

    # Write double species conditions
    for i in range(1, 13):
        for j in range(i+1, 13):
            mu_a = f'mu_{i}'
            mu_b = f'mu_{j}'
            alpha_ab = f'alpha_{i}_{j}'
            alpha_ba = f'alpha_{j}_{i}'
            alpha_aa = f'alpha_{i}_{i}'
            alpha_bb = f'alpha_{j}_{j}'
            for rep in ['i', 'ii', 'iii']:
                simulationConditionId = f'c_{i}_{j}_{rep}'
                tsv_writer.writerow([simulationConditionId, mu_a, mu_b, alpha_ab, alpha_ba, alpha_aa, alpha_bb])

print(f'Wrote to `{out_path}`.')


## Write parameter table
out_path = os.path.join(basedir, 'parameter_table.tsv')
with open(out_path, 'wt') as f:
    parameterScale = 'lin'
    estimate = 1
    objectivePriorType = 'normal'
    objectivePriorParameters = f'0.0; {1/50**0.5}' # Todo: @Sungho: What L2 priors and bounds have you used?

    tsv_writer = csv.writer(f, delimiter='\t')
    tsv_writer.writerow(['parameterId', 'parameterScale', 'lowerBound', 'upperBound', 'nominalValue', 'estimate', 'objectivePriorType', 'objectivePriorParameters'])

    # Write mu
    lowerBound = 0
    upperBound = 2
    nominalValue = 0.1
    for i in range(1, 13):
        parameterId = f'mu_{i}'
        tsv_writer.writerow([parameterId, parameterScale, lowerBound, upperBound, nominalValue, estimate, objectivePriorType, objectivePriorParameters])

    # Write alpha_ab and alpha_ba
    lowerBound = -2
    upperBound = 0
    nominalValue = -0.1
    for i in range(1, 13):
        for j in range(i+1, 13):
            parameterId = f'alpha_{i}_{j}'
            tsv_writer.writerow([parameterId, parameterScale, lowerBound, upperBound, nominalValue, estimate, objectivePriorType, objectivePriorParameters])
            parameterId = f'alpha_{j}_{i}'
            tsv_writer.writerow([parameterId, parameterScale, lowerBound, upperBound, nominalValue, estimate, objectivePriorType, objectivePriorParameters])


    # Write alpha_aa and alpha_bb
    lowerBound = -2
    upperBound = 0
    nominalValue = -0.1
    for i in range(1, 13):
        parameterId = f'alpha_{i}_{i}'
        tsv_writer.writerow([parameterId, parameterScale, lowerBound, upperBound, nominalValue, estimate, objectivePriorType, objectivePriorParameters])

    # Write sigma
    lowerBound = 0
    upperBound = 2
    nominalValue = 1
    estimate = 0
    objectivePriorType = ''
    objectivePriorParameters = ''
    parameterId = 'sigma'
    tsv_writer.writerow([parameterId, parameterScale, lowerBound, upperBound, nominalValue, estimate, objectivePriorType, objectivePriorParameters])

print(f'Wrote to `{out_path}`.')


## Write measurement table
out_path = os.path.join(basedir, 'measurement_table.tsv')
M_dir = os.path.join(basedir, 'raw_data', 'M')
data_dir = os.path.join(basedir, 'raw_data', 'data')

# Write single species experiments
m_exp = {'observableId': [], 'simulationConditionId': [], 'measurement': [], 'time': []}
for file in os.listdir(M_dir):
    df = pd.read_csv(os.path.join(M_dir, file))
    simulationConditionId = 'c_'+os.path.splitext(os.path.basename(file))[0] 
    time = list(df.iloc[:, 0])
    measurement = list(df.iloc[:, 1])
    m_exp['observableId'] = m_exp['observableId'] + ['A_OD600']*len(time)
    m_exp['simulationConditionId'] = m_exp['simulationConditionId'] + [simulationConditionId]*len(time)
    m_exp['measurement'] = m_exp['measurement'] + measurement
    m_exp['time'] = m_exp['time'] + time
df1 = pd.DataFrame(m_exp)
s = index_natsorted(df1['simulationConditionId'])
df1 = df1.iloc[s, :]

# Write double species experiments
m_exp = {'observableId': [], 'simulationConditionId': [], 'measurement': [], 'time': []}
for file in os.listdir(data_dir):
    df = pd.read_csv(os.path.join(data_dir, file))
    if '1.1' in df.columns:
        df = df.rename(columns={'1': '0', '1.1': '1'})
    s = natsorted(df.columns[1:])
    i = s[0]
    j = s[1]
    time = list(df.iloc[:, 0])
    measurement_a = list(df[i])
    measurement_b = list(df[j])
    simulationConditionId = f'c_{i}_{j}_i'
    if simulationConditionId in m_exp['simulationConditionId']:
        simulationConditionId = simulationConditionId+'i'
        if simulationConditionId in m_exp['simulationConditionId']:
            simulationConditionId = simulationConditionId+'i'
            if simulationConditionId in m_exp['simulationConditionId']:
                raise Exception('Found more than 3 replicates of same condition. This is not expected.')
    m_exp['observableId'] = m_exp['observableId'] + ['A_OD600']*len(time) + ['B_OD600']*len(time)
    m_exp['simulationConditionId'] = m_exp['simulationConditionId'] + [simulationConditionId]*2*len(time)
    m_exp['measurement'] = m_exp['measurement'] + measurement_a + measurement_b
    m_exp['time'] = m_exp['time'] + time*2

df2 = pd.DataFrame(m_exp)
s = index_natsorted(df2['simulationConditionId'])
df2 = df2.iloc[s, :]

df_final = df1.append(df2)
df_final.to_csv(out_path, sep='\t', index=False)

print(f'Wrote to `{out_path}`.')