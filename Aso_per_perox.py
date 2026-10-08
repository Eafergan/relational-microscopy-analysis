import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu


#
#
#Step 1: load the "Nuclei" database for 647 signal and...
#.. and the "549 Perox" database for the peroxide measurment
#... change the folder name to your specific folder
#
#
#

file_path_Nuclei = r"D:\Work\Einat\Evaluation2\Objects_Population - Nuclei.txt"
file_path_perox = r"D:\Work\Einat\Evaluation2\Objects_Population - 594 Perox.txt"

# Since the text file from the Opera has metadata section before the data block
# That Metadata needs to be ignored so pandas will focus only on the database
# So, first we find the line containing [Data]
with open(file_path_Nuclei, encoding="utf-8-sig") as file:
    skiprows = next(i for i, line in enumerate(file) if line.strip() == "[Data]") + 1

#then we load the database starting from that row
df = pd.read_csv(file_path_Nuclei, sep="\t", skiprows=skiprows)

#now for the Perox file
with open(file_path_perox, encoding="utf-8-sig") as file:
    skiprows_perox = next(i for i, line in enumerate(file) if line.strip() == "[Data]") + 1

df_perox = pd.read_csv(file_path_perox, sep="\t", skiprows=skiprows_perox)


#
#
#
#Step 2: identify transfected cells
#
#
#

df_no_aso = df[ (df["Column"].isin([2, 3,7,8]))].copy()
df_Aso40 = df[(df["Column"] == 9) & (df["Row"].isin([2, 3, 4]))].copy()
df_Aso120 = df[(df["Column"] == 9) & (df["Row"].isin([5, 6, 7]))].copy()

ColName = "Nuclei - Intensity Cell Alexa 647 Mean" #lets focus on this specific column


plt.boxplot(
    [df_no_aso[ColName].dropna(),
     df_Aso40[ColName].dropna(),
     df_Aso120[ColName].dropna()],
    
    tick_labels=["df_no_aso", "df_Aso40", "df_Aso120"]
)
plt.ylabel("Nuclei - Intensity Cell Alexa 647 Mean")
plt.show() #this plot clearly shows Aso transfected cells has higher 647 values

#we now draw a threshold at a 647 signal value where 99% of the non transfected cells are below that value
threshold = df_no_aso[ColName].quantile(0.99) 

Num_untransfected_over_threshold = (df_no_aso[ColName] > threshold).sum()
Num_transfected_over_threshold = (df_Aso120[ColName] > threshold).sum()
Percent_over_threshold_non_transfected = Num_untransfected_over_threshold / (Num_untransfected_over_threshold + (df_no_aso[ColName] < threshold).sum()) 
Percent_transfected_over_threshold = Num_transfected_over_threshold / (Num_transfected_over_threshold + + (df_Aso120[ColName] < threshold).sum() )

#now we build the same box plot but zoom in on the threshold area
plt.boxplot(
    [df_no_aso[ColName].dropna(),
     df_Aso40[ColName].dropna(),
     df_Aso120[ColName].dropna()],
    
    tick_labels=["No DNM1L Aso#3", "DNM1L Aso#3 40nM", "DNM1L Aso#3 120nM"]
)
plt.axhline(y=threshold, color='red', linestyle='--', linewidth=2)
plt.ylim(100, 300)
plt.ylabel("Nuclei - Intensity Cell Alexa 647 Mean")
plt.show()

#here is a list of un transfected cells above the threshold to look at in FIJI. Simply checkout row, column, field and XY in that database and use it to search the correct file
df_untransfected_over_threshold = df_no_aso[df_no_aso[ColName] > threshold].copy()
#and here is a list of transfected cells 3-fold higher than threshold
df_Transfected_3xover_threshold = df_Aso120[df_Aso120[ColName] > 3*threshold].copy()



#
#
#
#step 3: measure peroxide properties of transfected vs non transfected
#If needed- we can connect between the two datasets by "Object No" column in the "Nulcei" database with the "594 Perox - Object No in Nuclei" column in the "594 Perox" database
#
# 

df_positive = df_Aso120[df_Aso120[ColName] > threshold].copy()
df_negative = df_no_aso[df_no_aso[ColName] < threshold].copy()
df_negative_Aso = df_Aso120[df_Aso120[ColName] < threshold].copy()



test_col=21 #start with 21 for "Number of 594 Perox" until 27 "Perox Ratio Width to Length" run the entire block together
print("col name: " +df_positive.columns[test_col] )
Col_test_name=df_positive.columns[test_col]

plt.boxplot(
    [df_negative_Aso[Col_test_name].dropna(),
     df_positive[Col_test_name].dropna()],
    tick_labels=["Negative", "Positive"]
)
plt.ylabel(Col_test_name[9:-1])
plt.show() 

all_values = pd.concat([
    df_negative_Aso[Col_test_name],
    df_positive[Col_test_name]
]).dropna()


u_stat, p_value = mannwhitneyu(
    df_negative_Aso[Col_test_name].dropna(),
    df_positive[Col_test_name].dropna(),
    alternative="two-sided"
)

print("p-value:", p_value)

bins = np.histogram_bin_edges(all_values, bins=100)

plt.hist(
    df_negative_Aso[Col_test_name].dropna(),
    bins=bins,
    histtype="step",
    linewidth=2,
    density=True,
    label="Negative"
)

plt.hist(
    df_positive[Col_test_name].dropna(),
    bins=bins,
    histtype="step",
    linewidth=2,
    density=True,
    label="Positive"
)

plt.text(
    0.95, 0.95,
    f"Mann–Whitney p = {p_value:.2e}",
    transform=plt.gca().transAxes,
    ha="right",
    va="top"
)
plt.xlabel(Col_test_name[9:-1])
plt.ylabel("Density")
plt.legend()
plt.show()

print(
    f"Positive: {df_positive[Col_test_name].mean():.3f} ± "
    f"{df_positive[Col_test_name].sem():.3f} SEM\n"
    f"Negative: {df_negative_Aso[Col_test_name].mean():.3f} ± "
    f"{df_negative_Aso[Col_test_name].sem():.3f} SEM"
)



# step 4: single peroxisome tests
#   Here we merge the peroxisome table with the column of the 647 threshold measurment
#   That column will be used to seperate peroxisome from treated cells to non treated cells
#   to varify: 2 columns were added to the peroxisome df:
#   one column for the 647 measurment and te second is the object id in the nuclei table
#   use both new columns to varify the merge, the object number in the 647 df along with the row column and field shows the same numbers as in the nuclei table


marker_col = "Nuclei - Intensity Cell Alexa 647 Mean"

df_perox_labeled = df_perox.merge(
    df[
        ["Row", "Column", "Field", "Object No", marker_col]
    ],
    left_on=[
        "Row",
        "Column",
        "Field",
        "594 Perox - Object No in Nuclei"
    ],
    right_on=[
        "Row",
        "Column",
        "Field",
        "Object No"
    ],
    how="left",
    validate="many_to_one"
)

df_perox_labeled_aso120=df_perox_labeled[(df_perox_labeled["Column"] == 9) & (df_perox_labeled["Row"].isin([5, 6, 7]))].copy()
df_perox_labeled_No_Aso=df_perox_labeled[(df_perox_labeled["Column"] == 7) & (df_perox_labeled["Row"].isin([5, 6, 7]))].copy()


df_positive_perox = df_perox_labeled_aso120[
    df_perox_labeled_aso120[marker_col] > threshold
].copy()

df_negative_perox = df_perox_labeled_aso120[
    df_perox_labeled_aso120[marker_col] <= threshold
].copy()

#alternativly use this
df_negative_perox = df_perox_labeled_No_Aso[
    df_perox_labeled_No_Aso[marker_col] <= threshold
].copy()



test_col=27 #start with 25 for "Number of 594 Perox" until 30 "Perox Ratio Width to Length" run the entire block together
print("col name: " +df_perox_labeled.columns[test_col] )
Col_test_name=df_perox_labeled.columns[test_col]

print(df_perox[Col_test_name].nunique()) #for some columns, there are suspiciously low number of uniqu values

plt.boxplot(
    [df_positive_perox[Col_test_name].dropna(),
     df_negative_perox[Col_test_name].dropna()],
    
    tick_labels=["df_positive_perox", "df_negative_perox"]
)
plt.ylabel(Col_test_name)
plt.show() 

all_values_perox = pd.concat([
    df_negative_perox[Col_test_name],
    df_positive_perox[Col_test_name]
]).dropna()


u_stat, p_value_perox = mannwhitneyu(
    df_negative_perox[Col_test_name].dropna(),
    df_positive_perox[Col_test_name].dropna(),
    alternative="two-sided"
)


print("p-value:", p_value_perox)

bins_perox = np.histogram_bin_edges(all_values_perox, bins=200)

plt.hist(
    df_negative_perox[Col_test_name].dropna(),
    bins=bins_perox,
    histtype="step",
    linewidth=2,
    density=True,
    label="Negative"
)

plt.hist(
    df_positive_perox[Col_test_name].dropna(),
    bins=bins_perox,
    histtype="step",
    linewidth=2,
    density=True,
    label="Positive"
)

plt.text(
    0.95, 0.95,
    f"Mann–Whitney p = {p_value_perox:.2e}",
    transform=plt.gca().transAxes,
    ha="right",
    va="top"
)

plt.xlabel(Col_test_name[9:-1])
plt.ylabel("Density")
plt.legend()
plt.show()

print(
    f"Positive: {df_positive_perox[Col_test_name].mean():.3f} ± "
    f"{df_positive_perox[Col_test_name].sem():.3f} SEM\n"
    f"Negative: {df_negative_perox[Col_test_name].mean():.3f} ± "
    f"{df_negative_perox[Col_test_name].sem():.3f} SEM"
)

threshold_test= 5
df_positive_perox_test = df_positive_perox[
    df_positive_perox[Col_test_name] < threshold_test
].copy()

df_negative_perox_test = df_negative_perox[
    df_negative_perox[Col_test_name] < threshold_test
].copy()



all_values_perox_test = pd.concat([
    df_positive_perox_test[Col_test_name],
    df_negative_perox_test[Col_test_name]
]).dropna()



bins_perox_test = np.histogram_bin_edges(all_values_perox_test, bins=100)


plt.hist(
    df_negative_perox_test[Col_test_name].dropna(),
    bins=bins_perox_test,
    histtype="step",
    linewidth=2,
    density=True,
    label="Negative"
)

plt.hist(
    df_positive_perox_test[Col_test_name].dropna(),
    bins=bins_perox_test,
    histtype="step",
    linewidth=2,
    density=True,
    label="Positive"
)

plt.text(
    0.95, 0.95,
    f"Mann–Whitney p = {p_value_perox:.2e}",
    transform=plt.gca().transAxes,
    ha="right",
    va="top"
)

plt.xlabel(Col_test_name[9:-1])
plt.ylabel("Density")
plt.legend()
plt.show()


#use this to pinpoint specific cells or object
df_test2 = df[(df["Row"]==5) & (df["Column"] == 9) & (df["Field"] == 3) & (df["Object No"] == 1)  ].copy()
