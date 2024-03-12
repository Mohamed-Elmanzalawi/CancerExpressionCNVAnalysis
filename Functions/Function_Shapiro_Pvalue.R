####  A function calculates the shapiro p-values 

shapiro.pvalues=function(df){Shapiro_P=c()
n1=nrow(df)
for(i in 1:n1){Shapiro_P[i]=shapiro.test(as.numeric(df[i,]))$p.value}
return(Shapiro_P)}
