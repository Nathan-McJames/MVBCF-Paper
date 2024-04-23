
#Load Libraries
library(foreign)
library(naniar)
library(missRanger)
library(corrplot)
library(RColorBrewer)
library(bartCause)
library(stats)
library(jtools)
library(dplyr)
library(plyr)
library(bcf)
library(mltools)
library(stats)
library(ggplot2)
library(data.table)
library(ranger)
library(ggplot2)
library(ranger)
library(dplyr)
library(dbarts)
library(rpart)
library(rpart.plot)
library(Rcpp)

#Set Random Seed
set.seed(123)

setDTthreads(1)

#Load Data
BCG<-as.data.frame(read.spss("/YOUR_DIRECTORY/bcgirlm7.sav"))
BSG<-as.data.frame(read.spss("/YOUR_DIRECTORY/bsgirlm7.sav"))
BST<-as.data.frame(read.spss("/YOUR_DIRECTORY/bstirlm7.sav"))
BTM<-as.data.frame(read.spss("/YOUR_DIRECTORY/btmirlm7.sav"))
BTS<-as.data.frame(read.spss("/YOUR_DIRECTORY/btsirlm7.sav"))


#Each row is a student matched to a teacher.
#Some students are matched to >1 teacher.
#Some students have no teacher and don't appear.
BST_BTM<-merge(BST, BTM, by="IDTEALIN")

order1<-as.numeric(as.character(BST_BTM$IDSTUD))
order2<-as.numeric(as.character(BST_BTM$BTBM14))
order2<-ifelse(is.na(order2), 0, order2)*-1

BST_BTM<-BST_BTM[order(order1, order2),]

support_teachers<-c()

for(i in 2:(length(BST_BTM$IDSTUD)))
{
  if(BST_BTM$IDSTUD[i] == BST_BTM$IDSTUD[i-1])
  {
    support_teachers<-c(support_teachers, i)
  }
}

BST_BTM<-BST_BTM[-support_teachers,]

#Each row now has the context data added to it.
BSG_BST_BTM<-merge(BSG, BST_BTM, by="IDSTUD", all.x=T)



#Each row now has school data too.
BSG_BST_BTM_BCG<-merge(BSG_BST_BTM, BCG, by="IDSCHOOL")




BST_BTS<-merge(BST, BTS, by="IDTEALIN")

order1<-as.numeric(as.character(BST_BTS$IDSTUD))
order2<-as.numeric(as.character(BST_BTS$BTBS14))
order2<-ifelse(is.na(order2), 0, order2)*-1

BST_BTS<-BST_BTS[order(order1, order2),]

support_teachers<-c()

for(i in 2:(length(BST_BTS$IDSTUD)))
{
  if(BST_BTS$IDSTUD[i] == BST_BTS$IDSTUD[i-1])
  {
    support_teachers<-c(support_teachers, i)
  }
}

BST_BTS<-BST_BTS[-support_teachers,]

BSG_BST_BTS<-merge(BSG, BST_BTS, by="IDSTUD", all.x=T)

BSG_BST_BTS_BCG<-merge(BSG_BST_BTS, BCG, by="IDSCHOOL")




#Select Columns
maths_vars<-c("BSDAGE", "BSBG01", "BSBG03", "BSBG04", "BSBG07",
              "BSBG08A", "BSBG08B",
              "BSBG09A", "BSDGEDUP", "BSBGHER",
              "BSBGSSB", "BSBGSB", "BSBGICM",
              
              "BSBG05A", "BSBG05B", "BSBG05C", "BSBG05D", "BSBG05E", "BSBG05F", "BSBG05G", 
              "BSBG11B", "BSBG11A",
              
              "BTBG01", "BTBG02", "BTBG03", "BTBG10",
              "BTBGTJS", "BTBGSOS", "BTBGLSN", "BTBGEAS", "BTDMMME",
              
              "BCBGDAS", "BCBGEAS", "BCBGMRS",
              "BCDGSBC")


maths_other<-c("BSMMAT01.x", "BSMMAT02.x", "BSMMAT03.x", 
               "BSMMAT04.x", "BSMMAT05.x", "BSBGSCM", 
               "BSBGSVM", "IDCLASS.x", "TOTWGT", "IDSTUD")


maths_treatment<-c("BSBG10")


science_vars<-c("BSDAGE", "BSBG01", "BSBG03", "BSBG04", "BSBG07",
                "BSBG08A", "BSBG08B",
                "BSBG09A", "BSDGEDUP", "BSBGHER",
                "BSBGSSB", "BSBGSB", "BSBGICS",
                
                "BSBG05A", "BSBG05B", "BSBG05C", "BSBG05D", "BSBG05E", "BSBG05F", "BSBG05G", 
                "BSBG11B", "BSBG11A",
                
                "BTBG01", "BTBG02", "BTBG03", "BTBG10",
                "BTBGTJS", "BTBGSOS", "BTBGLSN", "BTBGEAS", "BTDSMSE",
                
                "BCBGDAS", "BCBGEAS", "BCBGSRS",
                "BCDGSBC")

science_other<-c("BSSSCI01.x", "BSSSCI02.x", "BSSSCI03.x",
                 "BSSSCI04.x", "BSSSCI05.x", "BSBGSCS", 
                 "BSBGSVS", "IDCLASS.x", "TOTWGT", "IDSTUD")

science_treatment<-c("BSBG10")



XYM<-BSG_BST_BTM_BCG[,c(maths_vars, maths_other, maths_treatment)]
XYS<-BSG_BST_BTS_BCG[,c(science_vars, science_other, science_treatment)]


#Function for converting to numeric
mymap<-function(x)
{
  x<-as.numeric(as.character(x))
  
  return(x)
}

#Categories to be mapped to numeric values
from<-c("Girl","Boy","Omitted or invalid",
        "Always","Almost always","Sometimes","Never",
        "None or very few (0–10 books)",                         
        "Enough to fill one shelf (11–25 books)",                
        "Enough to fill one bookcase (26–100 books)",            
        "Enough to fill two bookcases (101–200 books)",          
        "Enough to fill three or more bookcases (more than 200)",
        "Finish <Lower secondary education—ISCED Level 2>",                             
        "Finish <Upper secondary education—ISCED Level 3>",                             
        "Finish <Post-secondary, non-tertiary education—ISCED Level 4>",                
        "Finish <Short-cycle tertiary education—ISCED Level 5>",                        
        "Finish <Bachelor’s or equivalent level—ISCED Level 6>",                        
        "Finish <Postgraduate degree: Master’s—ISCED Level 7 or Doctor —ISCED Level 8>",
        "University or Higher",                      
        "Post-secondary but not University",         
        "Upper Secondary",                           
        "Lower Secondary",                           
        "Some Primary, Lower Secondary or No School",
        "Don't Know",
        "Yes","No","I don't know", "Not applicable",
        "Every day","Almost every day","Sometimes","Never",
        "Female","Male",
        "Under 25","25–29","30–39","40–49","50–59","60 or more",
        "Major in Mathematics and Mathematics Education",       
        "Major in Mathematics but not in Mathematics Education",
        "Major in Mathematics Education but not in Mathematics",
        "All Other Majors",                                     
        "No Formal Education Beyond Upper Secondary",
        "More Affluent",                               
        "Neither More Affluent nor More Disadvantaged",
        "More Disadvantaged",
        "Major in Science and Science Education",       
        "Major in Science but not in Science Education",
        "Major in Science Education but not in Science",
        "Once a week","Once every two weeks","Once a month","Once every two month","Never or almost never")


to=c(0, 1, NA,
     3, 2, 1, 0,
     5, 20, 50, 150, 200,
     2, 3, 4, 5, 6, 8,
     5, 4, 3, 2, 1, 0,
     1, 0, NA, NA,
     3, 2, 1, 0,
     0, 1,
     25, 29, 39, 49, 59, 66,
     2, 1, 1, 0, 0,
     3, 2, 1,
     2, 1, 1,
     4,3,2,1,0)



#Function for performing the conversion
mymap2<-function(x)
{
  x<-mapvalues(x,
               from=from,
               to=to)
  
  x<-as.numeric(as.character(x))
  
  return(x)
}

XYM<-mutate_at(XYM, c(2:9, 14:22, 24:25, 31, 35, 46), mymap2)
XYS<-mutate_at(XYS, c(2:9, 14:22, 24:25, 31, 35, 46), mymap2)

XYM<-mutate_at(XYM, c(1, 10:13, 23, 26:30, 32:34, 36:42, 44), mymap)
XYS<-mutate_at(XYS, c(1, 10:13, 23, 26:30, 32:34, 36:42, 44), mymap)


XY<-merge(XYM, XYS, by="IDSTUD")

XY<-missRanger(XY, num.threads=6, num.trees=50, pmm.k = 3)


X<-as.matrix(one_hot(as.data.table(XY[,-c(1, 37:41, 44:58, 60:68, 78:79, 81:86, 89:91)])))

Z1<-XY$BSBG10.x
Z2<-XY$BSBG10.y


Z1<-ifelse(XY$BSBG10.x>=3, 1, 0)

Z2<-ifelse(XY$BSBG10.y>=3, 1, 0)

Z<-cbind(Z1, Z2)


#Get Chain
#Change to BSMMT02.x, BSMMAT03.x etc for plausible values 2, 3, 4, 5
Y<-cbind(XY$BSMMAT05.x, XY$BSSSCI05.x)

n_tree_mu<-100
n_tree_tau<-50
n_iter<-5000
n_burn<-3000
keep_every<-2

#Propensity Score Estimate
p_mod<-bart(x.train = X, y.train = Z1, k=3)
p<-colMeans(pnorm(p_mod$yhat.train))

X2<-X[,c("BSBGHER.x", "BSDGEDUP.x", "BCDGSBC.x")]
X2<-cbind(X2, (X[,"BCBGMRS"]+X[,"BCBGSRS"])/2)
X2<-cbind(X2[,"BSDGEDUP.x"]==0, X2)
X<-cbind(X, p)
X<-cbind(X[,"BSDGEDUP.x"]==0, X)

XT<-X[,c(2:52)]

XT1<-rbind(XT,XT,XT,XT,XT,XT,XT,XT,XT,XT)
XT2<-rbind(XT,XT,XT,XT,XT,XT)

XT1[,"BSBGHER.x"]<-rep(seq(min(X[,"BSBGHER.x"]), max(X[,"BSBGHER.x"]), length.out=10), each=4118)
XT2[,"BSDGEDUP.x"]<-rep(c(0,1,2,3,4,5), each=4118)

XT<-rbind(XT1, XT2)
XT<-cbind(XT[,"BSDGEDUP.x"]==0, XT)

X2T<-XT[,c(c("BSBGHER.x", "BSDGEDUP.x", "BCDGSBC.x"))]
X2T<-cbind(X2T, (XT[,"BCBGMRS"]+XT[,"BCBGSRS"])/2)
X2T<-cbind(X2T[,"BSDGEDUP.x"]==0, X2T)

sourceCpp(file = "~/YOUR_DIRECTORY/MVBCF_RI_Code.cpp")

group_id<-as.integer(XY$IDCLASS.x.x) - 1

group_id_test<-rep(group_id, 16)

my_mod <- fast_bart(X, 
                    Y, 
                    Z, 
                    X2,
                    XT,
                    X2T,
                    0.95, 
                    2, 
                    0.25, 
                    3, 
                    diag((1)^2/n_tree_mu, 2), 
                    diag((0.3)^2/n_tree_tau, 2), 
                    1, 
                    diag(1, 2), 
                    n_iter, 
                    n_tree_mu, 
                    n_tree_tau, 
                    1, 
                    group_id,
                    group_id_test,
                    diag(0.1, 2),
                    matrix(0, nrow=2, ncol=1),
                    diag(0.01, 2),
                    1)


save(my_mod, file = "~/YOUR_DIRECTORY/Absent.RData", compress = "xz")

##########################
weights<-as.numeric(as.character(BSG$TOTWGT))
#bartcause applied to timss
colnames(X)<-paste0("V", 1:ncol(X))
bart_mod1<-bartc(Y[,1], Z[,1], X, p.scoreAsCovariate = F, n.chains=1, n.threads=1, keepTrees=T, n.trees=n_tree_mu+n_tree_tau,
                 group.by=as.factor(group_id),
                 use.ranef = T,
                 n.samples = 2000,
                 n.burn = 3000,
                 n.thin = 2)
bart_ate_preds1<-apply(predict(bart_mod1, X, type="icate", group.by=as.factor(group_id)), 1, weighted.mean, w=weights)

bart_mod2<-bartc(Y[,2], Z[,2], X, p.scoreAsCovariate = F, n.chains=1, n.threads=1, keepTrees=T, n.trees=n_tree_mu+n_tree_tau,
                 group.by=as.factor(group_id),
                 use.ranef = T,
                 n.samples = 2000,
                 n.burn = 3000,
                 n.thin = 2)
bart_ate_preds2<-apply(predict(bart_mod2, X, type="icate", group.by=as.factor(group_id)), 1, weighted.mean, w=weights)

#bcf applied to timss
my_mod <- fast_bart(X, 
                    as.matrix(Y[,1]), 
                    as.matrix(Z[,1]), 
                    X2,
                    X[1:3,],
                    X2[1:3,],
                    0.95, 
                    2, 
                    0.25, 
                    3, 
                    diag((1)^2/n_tree_mu, 1), 
                    diag((0.3)^2/n_tree_tau, 1), 
                    1, 
                    diag(1, 1), 
                    n_iter, 
                    n_tree_mu, 
                    n_tree_tau, 
                    1, 
                    group_id,
                    group_id[1:3],
                    diag(0.1, 1),
                    matrix(0, nrow=1, ncol=1),
                    diag(0.01, 1),
                    1)


bcf_ate_preds1<-apply(my_mod$predictions_tau[,1,iters_to_keep], 2, weighted.mean, w=weights)

my_mod <- fast_bart(X, 
                    as.matrix(Y[,2]), 
                    as.matrix(Z[,2]), 
                    X2,
                    X[1:3,],
                    X2[1:3,],
                    0.95, 
                    2, 
                    0.25, 
                    3, 
                    diag((1)^2/n_tree_mu, 1), 
                    diag((0.3)^2/n_tree_tau, 1), 
                    1, 
                    diag(1, 1), 
                    n_iter, 
                    n_tree_mu, 
                    n_tree_tau, 
                    1, 
                    group_id,
                    group_id[1:3],
                    diag(0.1, 1),
                    matrix(0, nrow=1, ncol=1),
                    diag(0.01, 1),
                    1)


bcf_ate_preds2<-apply(my_mod$predictions_tau[,1,iters_to_keep], 2, weighted.mean, w=weights)



tempdf<-data.frame(bart_ate_preds1, bart_ate_preds2,
                   bcf_ate_preds1, bcf_ate_preds2)

write.table(tempdf,"~/YOUR_DIRECTORY/Absent_ATES.csv", row.names = F, append=T, sep = ",", col.names=F)


