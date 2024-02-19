import pandas as pd
import matplotlib.pyplot as plt
import math

# Step 2: Load OTU Data
otu_data = pd.read_csv('../ ')

# Step 3: Explore the Data
print(otu_data.head())
print(otu_data.describe())
print(otu_data.info())

# Step 4: Data Cleaning
otu_data = otu_data.dropna()
otu_data = otu_data.drop_duplicates()
# Additional data cleaning tasks as needed

# Step 5: Data Transformation
otu_data['log_abundance'] = otu_data['abundance'].apply(lambda x: math.log(x) if x > 0 else 0)
# Additional data transformation tasks as needed

# Step 6: Data Visualization (Optional)
plt.hist(otu_data['abundance'], bins=20, color='blue', alpha=0.7)
plt.xlabel('Abundance')
plt.ylabel('Frequency')
plt.title('Distribution of OTU Abundances')
plt.show()
# Additional visualization tasks as needed

# Step 7: Save Processed Data
otu_data.to_csv('processed_otu_data.csv', index=False)
