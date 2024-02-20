import numpy as np

from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd


def are_ttest_conditions_satisfied(group1, group2, group3=None, alpha=0.05):
    # Check for Normality
    _, p_value1 = stats.shapiro(group1)
    _, p_value2 = stats.shapiro(group2)

    p_value3 = 1
    if group3 is not None:
        _, p_value3 = stats.shapiro(group3)

    # Check for Homogeneity of Variances
    if group3 is None:
        _, p_value = stats.levene(group1, group2)
    else:
        _, p_value = stats.levene(group1, group2, group3)

    return p_value1 > alpha and p_value2 > alpha and p_value3 > alpha and p_value > alpha


def test_multiple_groups(data, class_column, alpha = 0.05):
    groups = []
    for group in data[class_column].unique():
        groups.append(data[data[class_column] == group])

    columns_of_interest = data.drop(class_column, axis=1).columns

    # Perform ANOVA for each MRI feature
    results = {}
    for column in columns_of_interest:

        test_type = 'ANOVA'
        if are_ttest_conditions_satisfied(*[group[column] for group in groups]):
            _, p_value = stats.f_oneway(*[group[column] for group in groups])
        else:
            test_type = 'Kruskal-Wallis'
            _, p_value = stats.kruskal(*[group[column] for group in groups])

        results[column] = p_value, test_type

    # Sort the results by p-value in ascending order
    sorted_results = sorted(results.items(), key=lambda x: x[1])

    print(f"Most relevant features distinguishing groups: {data[class_column].unique()}")
    for feature, (p_value, test_type) in sorted_results:
        if p_value < alpha:
            print('-'*100)
            print(f"{feature}: [{test_type}] p-value = {p_value}")
            print('-'*100)

            tukey_results = pairwise_tukeyhsd(data[feature], data[class_column])
            print(f"\nTukey's test results for {feature}:")
            print(tukey_results)


def test_two_groups(data, group_column, alpha = 0.05):
    results = {}
    columns_of_interest = data.drop(group_column, axis=1).columns
    classes = np.sort(data[group_column].unique())
    
    for column in columns_of_interest:
        group1_data = data[column].loc[data[group_column] == classes[0]].dropna()
        group2_data = data[column].loc[data[group_column] == classes[1]].dropna()

        if len(group1_data) < 3 or len(group2_data) < 3:
            continue

        test_type = 't-test'
        if are_ttest_conditions_satisfied(group1_data, group2_data):
            _, p_value = stats.ttest_ind(group1_data, group2_data)
        else:
            test_type = 'Mann-Whitney U test'
            _, p_value = stats.mannwhitneyu(group1_data, group2_data)

        meandiff = group2_data.mean() - group1_data.mean()

        results[column] = p_value, meandiff, test_type

    sorted_results = sorted(results.items(), key=lambda x: x[1])

    found = False
    print(f'Test results for {group_column}:')
    print(f'Group counts: {len(group1_data)} | {len(group2_data)}')
    for feature, (p_value, meandiff, test_type) in sorted_results:
        if p_value < alpha:
            found = True
            print(f"- {feature}: [{test_type}] p-value = {np.round(p_value, 4)}; mean difference = {np.round(meandiff, 4)}")

    if not found:
        print('No significant differences found')

    return sorted_results


def test_paired_data(data, class_column, alpha = 0.05):
    results = {}
    columns_of_interest = data.drop(class_column, axis=1).columns

    ordered_classes = np.sort(data[class_column].unique())
    if len(ordered_classes) != 2:
        print('Only two classes are supported')
        return

    data_1 = data.loc[data[class_column] == ordered_classes[0]]
    data_2 = data.loc[data[class_column] == ordered_classes[1]]

    for column in columns_of_interest:
        differences = data_2[column] - data_1[column]
        if data_1[column].var() == 0 or data_2[column].var() == 0 or all(differences == 0):
            continue

        # Perform the t-test
        _, p_value = stats.ttest_rel(data_1[column], data_2[column])

        meandiff = data_2[column].mean() - data_1[column].mean()

        results[column] = p_value, meandiff

    sorted_results = sorted(results.items(), key=lambda x: x[1])

    print('t-test results:')
    for feature, (p_value, meandiff) in sorted_results:
        if p_value < alpha:
            print(f"- {feature}: p-value = {np.round(p_value, 4)}; mean difference = {np.round(meandiff, 4)}")