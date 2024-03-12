############
####   Function_1 to delete the gene.IDs column from each data
#####################################################################

delet.gene.Ids=function(df){new_df = df[,-1]
return(new_df)}