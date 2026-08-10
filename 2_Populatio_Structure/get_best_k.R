##best k (SOURCE:https://baylab.github.io/MarineGenomics/week-9-population-structure-using-ngsadmix.html)
#read in the data
data<-list.files(pattern = ".log", full.names = T)

#look at data to make sure it only has the log files
data

#use lapply to read in all our log files at once
bigData<-lapply(1:25, FUN = function(i) readLines(data[i]))
# find the line that starts with "best like=" or just "b"
library(stringr)

#this will pull out the line that starts with "b" from each file and return it as a list
foundset<-sapply(1:25, FUN= function(x) bigData[[x]][which(str_sub(bigData[[x]], 1, 1) == 'b')])
foundset

#now we need to pull out the first number in the string, we'll do this with the function sub
as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) )

#now lets store it in a dataframe
#make a dataframe with an index 1:7, this corresponds to our K values
logs<-data.frame(K = rep(1:5, each=5))

#add to it our likelihood values
logs$like<-as.vector(as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) ))

#and now we can calculate our delta K and probability
tapply(logs$like, logs$K, FUN= function(x) mean(abs(x))/sd(abs(x)))

