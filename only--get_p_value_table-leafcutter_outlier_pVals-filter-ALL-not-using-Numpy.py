def search_string_in_file(file_name, string_to_search):
    """Search for the given string in file and return lines containing that string 
    along with line numbers"""
    line_number = 0
    list_of_results = []
    # Open the file in read only mode
    with open(file_name, 'r') as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            line_number += 1
            if string_to_search in line:
                # If yes, then add the line number & line as a tuple in the list
                list_of_results.append((line_number, line.rstrip()))
                # Return list of tuples containing line numbers and lines where string is found
    return list_of_results

def get_p_value_effect_size_table(data_dir, p_value_file, effect_size_file, 
                                  p_value_threshold, effect_size_thresold, save_dir, save_name):
    
    import csv
    
    # p_value_file
    f1 = open(data_dir + "/" + p_value_file, 'r')
    p_value_matrix = f1.readlines()
    
    # effect_size_file
    #f2 = open(data_dir + "/" + effect_size_file, 'r')
    #effect_size_matrix = f2.readlines()
    
    # colunm names
    p_value_column_names = p_value_matrix[0].split()
    #effect_size_column_names = effect_size_matrix[0].split()
    #print(p_value_column_names) 
    
    # check the column names are the same
    #for i in range(len(p_value_column_names)):
    #   if p_value_column_names[i] != effect_size_column_names[i]:
    #        print("ERROR: the column names in two files are NOT the same!!!")
    
    #
    # 1st row is the column, the matrix/list start from 0
    # but the file lines start from 1
    
    #res = [("cluster_id", "genome_coordinate", "sample", "p_value", "effect_size")]
    res = []
    print(len(p_value_matrix))

    for n in range(1,len(p_value_matrix)):
    #for n in range(18,21):
        #print(str(n) + " in " + str(len(p_value_matrix)))

        p_value_line = p_value_matrix[n]

        # 1st column is the string of "cluster_id"
        # check each column/sample
        for i in range(2, len(p_value_line.split())):
            if float(p_value_line.split()[i]) < p_value_threshold:
                # matrix number from 0 and file line # from 1
                # target_lines[j][0] is the file line #
                # target_lines[j][0]-1 is the index in the matrix
                cluster_id = p_value_line.split()[1]
                start = float(cluster_id.split(":")[1])
                end = float(cluster_id.split(":")[2])
                res.append((p_value_line.split()[0], p_value_line.split()[1], "sample-" + p_value_column_names[i],
                            p_value_line.split()[i], end - start))
    print("before remove duplicated: " + str(len(res)))
    
    res = list(set(res))
	
    print("AFTER remove duplicated: " + str(len(res)))

    # write as csv file
    with open(save_dir + "/" + save_name, 'w', newline='') as save_res_file:
        res_writer = csv.writer(save_res_file, delimiter=',', 
                                quotechar='"', quoting=csv.QUOTE_MINIMAL)
        res_writer.writerow(("Gene_name", "cluster_id", "sample", "p_value", "effect_size"))
        for i in range(len(res)):
            res_writer.writerow(res[i])
        
    #print(res[1])
    #print(res)

    f1.close()
    #f2.close()
    return (len(res))

data_dir = "/projectsp/foran/Team_Chan/STAR_FUSION_results/FFPE-new-data-10-2019/code-EXON-outliers-analysis-May-3-2020/leafcutter/analysis/result/kinase_genes-TK-TKL"
save_dir = "/projectsp/foran/Team_Chan/STAR_FUSION_results/FFPE-new-data-10-2019/code-EXON-outliers-analysis-May-3-2020/leafcutter/analysis/result/kinase_genes-TK-TKL"
p_value_file = "kinase_genes-TK-TKL--leafleafcutter_outlier_pVals-all.txt"
effect_size_file = "leafcutter_outlier_effSize.txt"
save_name = "candidates---kinase_genes-TK-TKL--leafleafcutter_outlier_pVals-all.txt.csv"
p_value_threshold = 2
effect_size_thresold = -9
get_p_value_effect_size_table(data_dir, p_value_file, effect_size_file, p_value_threshold, effect_size_thresold, save_dir, save_name)
