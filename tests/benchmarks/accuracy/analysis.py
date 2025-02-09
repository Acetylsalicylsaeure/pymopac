import pandas as pd
import numpy as np

df = pd.read_csv('./result.csv')
df = df.dropna(subset=['SMILES'])
n = len(df)

# Calculate metrics
df['diff_API_MopacInput'] = abs(df['API'] - df['MopacInput'])
df['API_error'] = abs(df['API'] - df['Exp.'])
df['MopacInput_error'] = abs(df['MopacInput'] - df['Exp.'])

# Calculate statistics
metrics = {
    'diff_API_MopacInput': df['diff_API_MopacInput'],
    'API_error': df['API_error'],
    'MopacInput_error': df['MopacInput_error']
}

stats = {}
for name, metric in metrics.items():
    stats[name] = {
        'mean': metric.mean(),
        'std': metric.std(),
        'ci': 1.96 * metric.std() / np.sqrt(n)
    }

# Get top 10 differences
top10 = df.nlargest(10, 'diff_API_MopacInput')[
    ['Molecules', 'SMILES', 'diff_API_MopacInput']]

# Generate report
report = f"""# Error in Calculated Heat of Formation (kcal/mol) Predictions

Dataset: [PDDG](https://doi.org/10.1002/jcc.10162)

### Mean Average Error ± 95% CI
"""

for name, stat in stats.items():
    report += f"- {name}: {stat['mean']:.2f} ± {stat['ci']:.2f}\n\n"

report += "reference Error: 3.34, [source](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b01047)\n\n"

# Save report
with open('README.md', 'w') as f:
    f.write(report)
