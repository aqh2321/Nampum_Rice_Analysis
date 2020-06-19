import pandas as pd
import numpy as np
import skbio
import csv
import os
from datetime import date

class Diversity_Calcs(object):

    @staticmethod
    def print_avail_calcs()-> None:
        """
        prints out all available diviersity calculations
        :return:
        """
        print(skbio.diversity.get_alpha_diversity_metrics())
        print(skbio.diversity.get_beta_diversity_metrics())

    @staticmethod
    def output_diversity_results(current_file: str) -> pd.DataFrame:
        """

        :param current_file: CSV file matching our standard format
        :return: dataframe with useful diversity measurements
        additional columns can be added from the list output using Diversity_Calcs.print_avail_calcs.
        """
        div=pd.read_csv(current_file)
        with open(current_file) as file:
            header = file.readline().split(",")
        uniqueindex = header[0]

        #Preparing Diverity Calculations. Obtaining a clean dataframe without frequencies.
        div = div.drop(uniqueindex, axis=1)
        div = div.transpose()

        diversity_results= pd.DataFrame()
        diversity_results["Observed_Otus"]=skbio.diversity.alpha_diversity("observed_otus",div, ids=div.index)
        diversity_results["Simpsons_Index"]=skbio.diversity.alpha_diversity("simpson",div, ids=div.index)
        diversity_results["Simspons_Evenness_Measure_E"]=skbio.diversity.alpha_diversity("simpson_e",div, ids=div.index)
        diversity_results["Shannon_Entropy"]=skbio.diversity.alpha_diversity("shannon",div, ids=div.index)
        diversity_results["Fisher_Alpha"]=skbio.diversity.alpha_diversity("fisher_alpha",div, ids=div.index)
        return diversity_results

    @staticmethod
    def reyni_results(current_file: str) -> pd.DataFrame:
        """
        accepts a CSV file and returns a Dataframe with Rényi entropy values from 0 to 4
        :param current_file:
        :return: Dataframe
        """
        reyni = pd.read_csv(current_file)
        with open(current_file) as file:
            header = file.readline().split(",")
        uniqueindex = header[0]
        # Convert reyni dataframes to frequencies of each column
        reyni.set_index(uniqueindex, inplace=True)
        reynisums = pd.Series(reyni.sum())
        reynifreq = reyni / reynisums
        alphalist = [0, .99999, 2, 3, 4]
        reyniresults = {}
        for column in reynifreq:
            # print(reynifreq[column])
            listresults = []
            for alpha in alphalist:
                columnsum = 0
                for row in reynifreq[column]:
                    if row > 0:
                        columnsum += row ** alpha
                        # print(columnsum)
                    else:
                        continue
                columnsum = columnsum ** (1 / (1 - alpha))
                listresults.append(columnsum)
            reyniresults[column] = listresults
        reynidfl = pd.DataFrame.from_dict(data=reyniresults, orient='index', columns=alphalist)
        return reynidfl

    @staticmethod
    def kl_distance_cyto(current_file: str) -> None:
        """
        Accepts a CSV file in the CWD and calculates pair-wise KL distance for each sample pair,
        and then generates a format which can be interpreted by Cytoscape for visualization.

        :param current_file: CSV file of raw counts
        :return: creates a FSI file in the CWD
        """
        print("preparing output...")
        todaydate = date.today()
        outputname = ("kl_cyto_" + current_file[:-4] + str(todaydate) + ".fsi")
        with open(current_file) as file:
            header = file.readline().split(",")
        uniqueindex = header[0]

        KLdf = pd.read_csv(current_file)
        KLdf.set_index(uniqueindex, inplace=True)
        # adding one to every field for psuedocounting purposes
        KLdf = KLdf.astype("float64")
        for key, row in KLdf.iteritems():
            for key2, row2 in row.iteritems():
                KLdf.loc[key2][key] += 1

        for key, row in KLdf.iteritems():
            for key2, row2 in row.iteritems():
                if KLdf.loc[key2][key] <= 0:
                    print("error in", key2, key)

        # Convert dataframes to frequencies of each column
        for key, row in KLdf.iteritems():
            for key2, row2 in row.iteritems():
                value = KLdf.loc[key2][key]
                KLdf.loc[key2][key] = (value / sum(row))

        # calculating KL distance between two distrubtions. Distance(P||Q) = summation[Pi*Log2(Pi/Qi)]
        KLresults = {}
        for key, row in KLdf.iteritems():
            for key2, row2 in KLdf.iteritems():
                if key == key2:
                    continue
                data = row * np.log2(row / row2)
                resultkey = key + " " + key2
                reversekey = key2 + " " + key
                if reversekey in KLresults:
                    KLresults[reversekey] += sum(data)
                else:
                    KLresults[resultkey] = sum(data)
        resultframe = pd.DataFrame.from_dict(KLresults, orient='index', columns=["KL_Distance"])
        resultframe.reset_index(inplace=True)

        # splitting dictionary key column into two seperate columns:
        # new data frame with split index columns
        new = resultframe["index"].str.split(" ", n=1, expand=True)
        # making separate columns from new data frame
        resultframe["Node_1"] = new[0]
        resultframe["Node_2"] = new[1]

        # Dropping old index columns
        resultframe.drop(columns=["index"], inplace=True)
        resultframe = resultframe[['Node_1', 'Node_2', "KL_Distance"]]

        f = open(outputname, "w")
        f.write(resultframe.to_csv(index=False, sep='\t'))
        f.close()
        print(f"File {outputname} has been succesfully created")