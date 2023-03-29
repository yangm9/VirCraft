import pandas as pd

# Load the VIBRANT results file
vibrant_results = pd.read_csv()

# Filter for AMGs with T, B flag and auxiliary scores > 3
amgs_to_remove = vibrant_results[(vibrant_results['Classification'] == 'AMG') &                                 (vibrant_results['T flag'] == 'True') & 
                                 (vibrant_results['B flag'] == 'True') & 
                                 (vibrant_results['Auxiliary score'] > 3)]

# Remove the identified AMGs from the VIBRANT results
vibrant_results = vibrant_results[~vibrant_results.index.isin(amgs_to_remove.index)]

# Save the filtered results to a new file
vibrant_results.to_csv('filtered_vibrant_results.csv', index=False)
