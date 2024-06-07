library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(bigmemory)
library(DoubletFinder)

#code####
ABC.data <- ABC.data$`Gene Expression`
#####
setwd("G:/work/test")
H02001=readRDS("./H02001/H02001.rds")
L02001=readRDS("./L02001/L02001.rds")
H0000101=readRDS("./H0000101/H0000101.rds")
L0000102=readRDS("./L0000102/L0000102.rds")
MCD_16=readRDS("./MCD_16/MCD_16.rds")
normal_16=readRDS("./normal_16/normal_16.rds")
merge = readRDS("./merge.rds")
#FCD*8#########################################################
ABC.data <- Read10X(data.dir = r"(G:\filterd data\5\FCD\filtered_feature_bc_matrix)")
ABC.data <- ABC.data$`Gene Expression`
A1 <- CreateSeuratObject(counts = ABC.data, project = "F5", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\5\noFCD\filtered_feature_bc_matrix)")
ABC.data <- ABC.data$`Gene Expression`
A2 <- CreateSeuratObject(counts = ABC.data, project = "nF5", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\6\tsc01\filtered_feature_bc_matrix)")
ABC.data <- ABC.data$`Gene Expression`
A3 <- CreateSeuratObject(counts = ABC.data, project = "F6", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\6\notsc\filtered_feature_bc_matrix)")
ABC.data <- ABC.data$`Gene Expression`
A4 <- CreateSeuratObject(counts = ABC.data, project = "nF6", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\7\FCD_211118\filtered_feature_bc_matrix)")
ABC.data <- ABC.data$`Gene Expression`
A5 <- CreateSeuratObject(counts = ABC.data, project = "F7", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\7\noFCD_211118\filtered_feature_bc_matrix)")
ABC.data <- ABC.data$`Gene Expression`
A6 <- CreateSeuratObject(counts = ABC.data, project = "nF7", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\8\FCD220217\2.2.filtered_feature_bc_matrix)")
A7 <- CreateSeuratObject(counts = ABC.data, project = "F8", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\8\noFCD220217\2.2.filtered_feature_bc_matrix)")
A8 <- CreateSeuratObject(counts = ABC.data, project = "nF8", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\18\FCD_220811\2.2.filtered_feature_bc_matrix)")
A9 <- CreateSeuratObject(counts = ABC.data, project = "F18", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\18\noFCD_220811\2.2.filtered_feature_bc_matrix)")
A10<- CreateSeuratObject(counts = ABC.data, project = "nF18", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\19\FCD_20220815\2.2.filtered_feature_bc_matrix)")
A11<- CreateSeuratObject(counts = ABC.data, project = "F19", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\19\noFCD_20220815\2.2.filtered_feature_bc_matrix)")
A12<- CreateSeuratObject(counts = ABC.data, project = "nF19", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\21\220830_FCD\filtered_tf_bc_matrix)")
ABC.data <- ABC.data$`Gene Expression`
A13<- CreateSeuratObject(counts = ABC.data, project = "F21", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\21\220830_noFCD\filtered_tf_bc_matrix)")
ABC.data <- ABC.data$`Gene Expression`
A14<- CreateSeuratObject(counts = ABC.data, project = "nF21", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\23\220830_FCD\filtered_feature_bc_matrix)")
ABC.data <- ABC.data$`Gene Expression`
A15<- CreateSeuratObject(counts = ABC.data, project = "F23", min.cells = 3, min.features = 200)


ABC.data <- Read10X(data.dir = r"(G:\filterd data\23\220830_noFCD\filtered_feature_bc_matrix)")
ABC.data <- ABC.data$`Gene Expression`
A16<- CreateSeuratObject(counts = ABC.data, project = "nF23", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\24\MCD\filtered_feature_bc_matrix)")
ABC.data <- ABC.data$`Gene Expression`
A17<- CreateSeuratObject(counts = ABC.data, project = "F24", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\24\noMCD\filtered_feature_bc_matrix)")
ABC.data <- ABC.data$`Gene Expression`
A18<- CreateSeuratObject(counts = ABC.data, project = "nF24", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\26\211117_FCD\filtered_feature_bc_matrix)")
A19<- CreateSeuratObject(counts = ABC.data, project = "F26", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\26\211117_noFCD\filtered_feature_bc_matrix)")
A20<- CreateSeuratObject(counts = ABC.data, project = "nF26", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\34\230220FCD)")
A21<- CreateSeuratObject(counts = ABC.data, project = "F34", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\34\230220noFCD)")
A22<- CreateSeuratObject(counts = ABC.data, project = "nF34", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\35\230221FCD)")
A23<- CreateSeuratObject(counts = ABC.data, project = "F35", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\35\230221noFCD)")
A24<- CreateSeuratObject(counts = ABC.data, project = "nF35", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\42\230315F\outs\filtered_feature_bc_matrix)")
A25<- CreateSeuratObject(counts = ABC.data, project = "F42", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\42\230315nF\outs\filtered_feature_bc_matrix)")
A26<- CreateSeuratObject(counts = ABC.data, project = "nF42", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\43\230321F\outs\filtered_feature_bc_matrix)")
A27<- CreateSeuratObject(counts = ABC.data, project = "F43", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\43\230321NOF\outs\filtered_feature_bc_matrix)")
A28<- CreateSeuratObject(counts = ABC.data, project = "nF43", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\51\230509FCD\filtered_feature_bc_matrix)")
A29<- CreateSeuratObject(counts = ABC.data, project = "F51", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\51\230509noFCD\filtered_feature_bc_matrix)")
A30<- CreateSeuratObject(counts = ABC.data, project = "nF51", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\52\230519FCD\filtered_feature_bc_matrix)")
A31<- CreateSeuratObject(counts = ABC.data, project = "F52", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\52\230519noFCD\filtered_feature_bc_matrix)")
A32<- CreateSeuratObject(counts = ABC.data, project = "nF52", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\54\fcd230704\filtered_feature_bc_matrix)")
A33<- CreateSeuratObject(counts = ABC.data, project = "F54", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\54\nofcd230704\filtered_feature_bc_matrix)")
A34<- CreateSeuratObject(counts = ABC.data, project = "nF54", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\55\KYSY-11115_221218_result\result\FCD2398\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A35<- CreateSeuratObject(counts = ABC.data, project = "F55", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\55\KYSY-11115_221218_result\result\noFCD2398\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A36<- CreateSeuratObject(counts = ABC.data, project = "nF55", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\57\KYSY-11243_221218_result\result\230921F\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A37<- CreateSeuratObject(counts = ABC.data, project = "F57", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\57\KYSY-11243_221218_result\result\230921nF\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A38<- CreateSeuratObject(counts = ABC.data, project = "nF57", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\59\KYSY-11549_221218_result\result\23111f\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A39<- CreateSeuratObject(counts = ABC.data, project = "F59", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\59\KYSY-11549_221218_result\result\23111nof\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A40<- CreateSeuratObject(counts = ABC.data, project = "nF59", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\64\KYSY-12041_221218_result\result\231226F\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A41<- CreateSeuratObject(counts = ABC.data, project = "F64", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\64\KYSY-12041_221218_result\result\231226nF\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A42<- CreateSeuratObject(counts = ABC.data, project = "nF64", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"()")
A43<- CreateSeuratObject(counts = ABC.data, project = "F66", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"()")
A44<- CreateSeuratObject(counts = ABC.data, project = "nF66", min.cells = 3, min.features = 200)

####no atac
ABC.data <- Read10X(data.dir = r"(G:\filterd data\8\FCD220217\2.2.filtered_feature_bc_matrix)")
A1 <- CreateSeuratObject(counts = ABC.data, project = "F8", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\8\noFCD220217\2.2.filtered_feature_bc_matrix)")
A2 <- CreateSeuratObject(counts = ABC.data, project = "nF8", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\18\FCD_220811\2.2.filtered_feature_bc_matrix)")
A3 <- CreateSeuratObject(counts = ABC.data, project = "F18", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\18\noFCD_220811\2.2.filtered_feature_bc_matrix)")
A4<- CreateSeuratObject(counts = ABC.data, project = "nF18", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\19\FCD_20220815\2.2.filtered_feature_bc_matrix)")
A5<- CreateSeuratObject(counts = ABC.data, project = "F19", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\19\noFCD_20220815\2.2.filtered_feature_bc_matrix)")
A6<- CreateSeuratObject(counts = ABC.data, project = "nF19", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\26\211117_FCD\filtered_feature_bc_matrix)")
A7<- CreateSeuratObject(counts = ABC.data, project = "F26", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\26\211117_noFCD\filtered_feature_bc_matrix)")
A8<- CreateSeuratObject(counts = ABC.data, project = "nF26", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\34\230220FCD)")
A9<- CreateSeuratObject(counts = ABC.data, project = "F34", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\34\230220noFCD)")
A10<- CreateSeuratObject(counts = ABC.data, project = "nF34", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\35\230221FCD)")
A11<- CreateSeuratObject(counts = ABC.data, project = "F35", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\35\230221noFCD)")
A12<- CreateSeuratObject(counts = ABC.data, project = "nF35", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\42\230315F\outs\filtered_feature_bc_matrix)")
A13<- CreateSeuratObject(counts = ABC.data, project = "F42", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\42\230315nF\outs\filtered_feature_bc_matrix)")
A14<- CreateSeuratObject(counts = ABC.data, project = "nF42", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\43\230321F\outs\filtered_feature_bc_matrix)")
A15<- CreateSeuratObject(counts = ABC.data, project = "F43", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\43\230321NOF\outs\filtered_feature_bc_matrix)")
A16<- CreateSeuratObject(counts = ABC.data, project = "nF43", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\51\230509FCD\filtered_feature_bc_matrix)")
A17<- CreateSeuratObject(counts = ABC.data, project = "F51", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\51\230509noFCD\filtered_feature_bc_matrix)")
A18<- CreateSeuratObject(counts = ABC.data, project = "nF51", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\52\230519FCD\filtered_feature_bc_matrix)")
A19<- CreateSeuratObject(counts = ABC.data, project = "F52", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\52\230519noFCD\filtered_feature_bc_matrix)")
A20<- CreateSeuratObject(counts = ABC.data, project = "nF52", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\54\fcd230704\filtered_feature_bc_matrix)")
A21<- CreateSeuratObject(counts = ABC.data, project = "F54", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\54\nofcd230704\filtered_feature_bc_matrix)")
A22<- CreateSeuratObject(counts = ABC.data, project = "nF54", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\55\KYSY-11115_221218_result\result\FCD2398\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A23<- CreateSeuratObject(counts = ABC.data, project = "F55", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\55\KYSY-11115_221218_result\result\noFCD2398\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A24<- CreateSeuratObject(counts = ABC.data, project = "nF55", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\57\KYSY-11243_221218_result\result\230921F\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A25<- CreateSeuratObject(counts = ABC.data, project = "F57", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\57\KYSY-11243_221218_result\result\230921nF\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A26<- CreateSeuratObject(counts = ABC.data, project = "nF57", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\59\KYSY-11549_221218_result\result\23111f\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A27<- CreateSeuratObject(counts = ABC.data, project = "F59", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\59\KYSY-11549_221218_result\result\23111nof\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A28<- CreateSeuratObject(counts = ABC.data, project = "nF59", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\64\KYSY-12041_221218_result\result\231226F\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A29<- CreateSeuratObject(counts = ABC.data, project = "F64", min.cells = 3, min.features = 200)
ABC.data <- Read10X(data.dir = r"(G:\filterd data\64\KYSY-12041_221218_result\result\231226nF\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A30<- CreateSeuratObject(counts = ABC.data, project = "nF64", min.cells = 3, min.features = 200)

#TLE*6####################
ABC.data <- Read10X(data.dir = r"(G:\filterd data\61\KYSY-11967_221218_result\result\Ctrl_1\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A1 <- CreateSeuratObject(counts = ABC.data, project = "Ctrl_1", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\61\KYSY-11967_221218_result\result\Ctrl_2\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A2 <- CreateSeuratObject(counts = ABC.data, project = "Ctrl_2", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\61\KYSY-11967_221218_result\result\Ctrl_3\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A3 <- CreateSeuratObject(counts = ABC.data, project = "Ctrl_3", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\61\KYSY-11967_221218_result\result\TLE_1\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A4 <- CreateSeuratObject(counts = ABC.data, project = "TLE_1", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\61\KYSY-11967_221218_result\result\TLE_2\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A5 <- CreateSeuratObject(counts = ABC.data, project = "TLE_2", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\61\KYSY-11967_221218_result\result\TLE_3\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A6 <- CreateSeuratObject(counts = ABC.data, project = "TLE_3", min.cells = 3, min.features = 200)

#HS*11####################################
ABC.data <- Read10X(data.dir = r"(G:\filterd data\28\230109HSc\filtered_feature_bc_matrix)")
A1 <- CreateSeuratObject(counts = ABC.data, project = "HS28c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\37\230224HSc\filtered_feature_bc_matrix)")
A2 <- CreateSeuratObject(counts = ABC.data, project = "HS37c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\41\230314HSc\filtered_feature_bc_matrix)")
A3 <- CreateSeuratObject(counts = ABC.data, project = "HS41c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\49\230423HSc\filtered_feature_bc_matrix)")
A4 <- CreateSeuratObject(counts = ABC.data, project = "HS49c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\53\230526HSc\filtered_feature_bc_matrix)")
A5 <- CreateSeuratObject(counts = ABC.data, project = "HS53c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\30\230203HSc)")
A6 <- CreateSeuratObject(counts = ABC.data, project = "HS30c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\33\230217HSc)")
A7 <- CreateSeuratObject(counts = ABC.data, project = "HS33c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\36\230223HSc)")
A8 <- CreateSeuratObject(counts = ABC.data, project = "HS36c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\39\230307HSc\filtered_feature_bc_matrix)")
A9 <- CreateSeuratObject(counts = ABC.data, project = "HS39c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\40\2303091HSc\filtered_feature_bc_matrix)")
A10<- CreateSeuratObject(counts = ABC.data, project = "HS40c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\44\230324HSc\outs\filtered_feature_bc_matrix)")
A11<- CreateSeuratObject(counts = ABC.data, project = "HS44c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\46\230404HSc\outs\filtered_feature_bc_matrix)")
A12<- CreateSeuratObject(counts = ABC.data, project = "HS46c", min.cells = 3, min.features = 200)


#HSn+c*8######################################
ABC.data <- Read10X(data.dir = r"(G:\filterd data\30\230203HSc)")
A1 <- CreateSeuratObject(counts = ABC.data, project = "230203c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\30\230203HSn)")
A2 <- CreateSeuratObject(counts = ABC.data, project = "230203n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\32\230214HSn)")
A3 <- CreateSeuratObject(counts = ABC.data, project = "230214n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\33\230217HSc)")
A4 <- CreateSeuratObject(counts = ABC.data, project = "230217c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\33\230217HSn)")
A5 <- CreateSeuratObject(counts = ABC.data, project = "230217n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\36\230223HSc)")
A6 <- CreateSeuratObject(counts = ABC.data, project = "230223c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\36\230223HSn)")
A7 <- CreateSeuratObject(counts = ABC.data, project = "230223n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"()")
A8 <- CreateSeuratObject(counts = ABC.data, project = "", min.cells = 3, min.features = 200)

#HS*ALL##################
ABC.data <- Read10X(data.dir = r"(G:\filterd data\58\KA1\filtered_feature_bc_matrix)")
A1<- CreateSeuratObject(counts = ABC.data, project = "KA1", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\58\KA2\filtered_feature_bc_matrix)")
A2<- CreateSeuratObject(counts = ABC.data, project = "KA2", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\58\KA3\filtered_feature_bc_matrix)")
A3<- CreateSeuratObject(counts = ABC.data, project = "KA3", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\58\N_1\filtered_feature_bc_matrix)")
A4<- CreateSeuratObject(counts = ABC.data, project = "N_1", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\58\N_2\filtered_feature_bc_matrix)")
A5<- CreateSeuratObject(counts = ABC.data, project = "N_2", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\58\N_4\filtered_feature_bc_matrix)")
A6<- CreateSeuratObject(counts = ABC.data, project = "N_4", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\61\KYSY-11967_221218_result\result\TLE_1\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A7<- CreateSeuratObject(counts = ABC.data, project = "TLE_1", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\61\KYSY-11967_221218_result\result\TLE_2\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A8<- CreateSeuratObject(counts = ABC.data, project = "TLE_2", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\61\KYSY-11967_221218_result\result\TLE_3\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A9<- CreateSeuratObject(counts = ABC.data, project = "TLE_3", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\61\KYSY-11967_221218_result\result\Ctrl_1\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A10<- CreateSeuratObject(counts = ABC.data, project = "Ctrl_1", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\61\KYSY-11967_221218_result\result\Ctrl_2\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A11<- CreateSeuratObject(counts = ABC.data, project = "Ctrl_2", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\61\KYSY-11967_221218_result\result\Ctrl_3\2.Basic_analysis\2.2.filtered_feature_bc_matrix)")
A12<- CreateSeuratObject(counts = ABC.data, project = "Ctrl_3", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\15\HS220506_1\2.2.filtered_feature_bc_matrix)")
A13<- CreateSeuratObject(counts = ABC.data, project = "HS15.1", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\15\HS220506_2\2.2.filtered_feature_bc_matrix)")
A14<- CreateSeuratObject(counts = ABC.data, project = "HS15.2", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\16\HS_220510_1\2.2.filtered_feature_bc_matrix)")
A15<- CreateSeuratObject(counts = ABC.data, project = "HS16.1", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\16\HS_220510_2\2.2.filtered_feature_bc_matrix)")
A16<- CreateSeuratObject(counts = ABC.data, project = "HS16.2", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\17\HS1_220804\2.2.filtered_feature_bc_matrix)")
A17<- CreateSeuratObject(counts = ABC.data, project = "HS17.1", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\20\220822_HS1\2.2.filtered_feature_bc_matrix)")
A18<- CreateSeuratObject(counts = ABC.data, project = "HS20.1", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\20\220822_HS2\2.2.filtered_feature_bc_matrix)")
A19<- CreateSeuratObject(counts = ABC.data, project = "HS20.2", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\27\HS1\filtered_feature_bc_matrix)")
A20<- CreateSeuratObject(counts = ABC.data, project = "HS27.1", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\27\HS2\filtered_feature_bc_matrix)")
A21<- CreateSeuratObject(counts = ABC.data, project = "HS27.2", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\28\230109HSc\filtered_feature_bc_matrix)")
A22<- CreateSeuratObject(counts = ABC.data, project = "HS28c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\28\230109HSn\filtered_feature_bc_matrix)")
A23<- CreateSeuratObject(counts = ABC.data, project = "HS28n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\29\HSC_cellranger\filtered_feature_bc_matrix)")
A24<- CreateSeuratObject(counts = ABC.data, project = "HS29c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\29\HSN_cellranger\filtered_feature_bc_matrix)")
A25<- CreateSeuratObject(counts = ABC.data, project = "HS29n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\30\230203HSc)")
A26<- CreateSeuratObject(counts = ABC.data, project = "HS30c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\30\230203HSn)")
A27<- CreateSeuratObject(counts = ABC.data, project = "HS30n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\32\230214HSn)")
A28<- CreateSeuratObject(counts = ABC.data, project = "HS32.1", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\33\230217HSc)")
A29<- CreateSeuratObject(counts = ABC.data, project = "HS33c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\33\230217HSn)")
A30<- CreateSeuratObject(counts = ABC.data, project = "HS33n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\36\230223HSc)")
A31<- CreateSeuratObject(counts = ABC.data, project = "HS36c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\36\230223HSn)")
A32<- CreateSeuratObject(counts = ABC.data, project = "HS36n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\37\230224HSc\filtered_feature_bc_matrix)")
A33<- CreateSeuratObject(counts = ABC.data, project = "HS37c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\37\230224HSn\filtered_feature_bc_matrix)")
A34<- CreateSeuratObject(counts = ABC.data, project = "HS37n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\38\230303HSc\filtered_feature_bc_matrix)")
A35<- CreateSeuratObject(counts = ABC.data, project = "HS38c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\39\230307HSc\filtered_feature_bc_matrix)")
A36<- CreateSeuratObject(counts = ABC.data, project = "HS39c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\39\230307HSn\filtered_feature_bc_matrix)")
A37<- CreateSeuratObject(counts = ABC.data, project = "HS39n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\40\2303091HSc\filtered_feature_bc_matrix)")
A38<- CreateSeuratObject(counts = ABC.data, project = "HS40c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\40\2303091HSn\filtered_feature_bc_matrix)")
A39<- CreateSeuratObject(counts = ABC.data, project = "HS40n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\41\230314HSc\filtered_feature_bc_matrix)")
A40<- CreateSeuratObject(counts = ABC.data, project = "HS41c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\41\230314HSn\filtered_feature_bc_matrix)")
A41<- CreateSeuratObject(counts = ABC.data, project = "HS41n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\44\230324HSc\outs\filtered_feature_bc_matrix)")
A42<- CreateSeuratObject(counts = ABC.data, project = "HS44c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\44\230324HSn\outs\filtered_feature_bc_matrix)")
A43<- CreateSeuratObject(counts = ABC.data, project = "HS44n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\45\230328HSn\outs\filtered_feature_bc_matrix)")
A44<- CreateSeuratObject(counts = ABC.data, project = "HS45n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\46\230404HSc\outs\filtered_feature_bc_matrix)")
A45<- CreateSeuratObject(counts = ABC.data, project = "HS46c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\46\230404HSn\outs\filtered_feature_bc_matrix)")
A46<- CreateSeuratObject(counts = ABC.data, project = "HS46n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\47\230413HSc\filtered_feature_bc_matrix)")
A47<- CreateSeuratObject(counts = ABC.data, project = "HS47c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\47\230413HSn\filtered_feature_bc_matrix)")
A48<- CreateSeuratObject(counts = ABC.data, project = "HS47n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\48\230420HSn\filtered_feature_bc_matrix)")
A49<- CreateSeuratObject(counts = ABC.data, project = "HS48n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\49\230423HSc\filtered_feature_bc_matrix)")
A50<- CreateSeuratObject(counts = ABC.data, project = "HS49c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\49\230423HSn\filtered_feature_bc_matrix)")
A51<- CreateSeuratObject(counts = ABC.data, project = "HS49n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\50\230506HSC\filtered_feature_bc_matrix)")
A52<- CreateSeuratObject(counts = ABC.data, project = "HS50c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\50\230506HSN\filtered_feature_bc_matrix)")
A53<- CreateSeuratObject(counts = ABC.data, project = "HS50n", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\53\230526HSc\filtered_feature_bc_matrix)")
A54<- CreateSeuratObject(counts = ABC.data, project = "HS53c", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = r"(G:\filterd data\53\230526HSn\filtered_feature_bc_matrix)")
A55<- CreateSeuratObject(counts = ABC.data, project = "HS53n", min.cells = 3, min.features = 200)



#A* to brain list()####################################
rm(ABC.data)

gc()

brain = list()
brain[[1]]=A1
brain[[2]]=A2
brain[[3]]=A3
brain[[4]]=A4
brain[[5]]=A5
brain[[6]]=A6
brain[[7]]=A7
brain[[8]]=A8
brain[[9]]=A9
brain[[10]]=A10
brain[[11]]=A11
brain[[12]]=A12
brain[[13]]=A13
brain[[14]]=A14
brain[[15]]=A15
brain[[16]]=A16
brain[[17]]=A17
brain[[18]]=A18
brain[[19]]=A19
brain[[20]]=A20
brain[[21]]=A21
brain[[22]]=A22
brain[[23]]=A23
brain[[24]]=A24
brain[[25]]=A25
brain[[26]]=A26
brain[[27]]=A27
brain[[28]]=A28
brain[[29]]=A29
brain[[30]]=A30
brain[[31]]=A31
brain[[32]]=A32
brain[[33]]=A33
brain[[34]]=A34
brain[[35]]=A35
brain[[36]]=A36
brain[[37]]=A37
brain[[38]]=A38
brain[[39]]=A39
brain[[40]]=A40
brain[[41]]=A41
brain[[42]]=A42
brain[[43]]=A43
brain[[44]]=A44
brain[[45]]=A45
brain[[46]]=A46
brain[[47]]=A47
brain[[48]]=A48
brain[[49]]=A49
brain[[50]]=A50
brain[[51]]=A51
brain[[52]]=A52
brain[[53]]=A53
brain[[54]]=A54
brain[[55]]=A55
brain[[56]]=A56
brain[[57]]=A57
brain[[58]]=A58
brain[[59]]=A59


rm(list = ls()[1:length(brain)])

#merge,&harmony########################
Find_doublet <- function(data){
  sweep.res.list <- paramSweep_v3(data, PCs = 1:50, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  nExp_poi <- round(0.05*ncol(data))
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  data <- doubletFinder_v3(data, PCs = 1:50, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  return(data)
}
for (i in 1:length(brain)) {
  brain[[i]][["percent.mt"]] <- PercentageFeatureSet(brain[[i]], pattern = "^MT-")
  brain[[i]] <- subset(brain[[i]], subset = nFeature_RNA > 200 & percent.mt < 10)
  brain[[i]] <- NormalizeData(brain[[i]])
  brain[[i]] <- FindVariableFeatures(brain[[i]], selection.method = "vst", nfeatures = 3000)
  brain[[i]] <- ScaleData(brain[[i]])
  brain[[i]] <- RunPCA(brain[[i]])
  brain[[i]] <- RunUMAP(brain[[i]], dims = 1:50)
  brain[[i]]<-Find_doublet(brain[[i]])
  brain[[i]]<-subset(brain[[i]],subset=doublet_info=="Singlet")
  c <- grep("pANN_",colnames(brain[[i]]@meta.data))
  brain[[i]]@meta.data <- brain[[i]]@meta.data[,-c]
}
list = brain
list[[1]]=NULL
brain.merge = merge(brain[[1]],y = list)

brain.merge <- NormalizeData(brain.merge)
brain.merge <- FindVariableFeatures(brain.merge, selection.method = "vst", nfeatures = 5000)
brain.merge <- ScaleData(brain.merge, verbose = FALSE)
brain.merge <- RunPCA(brain.merge, npcs = 50, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, reduction = "pca", dims = 1:50)
brain.merge <- FindClusters(brain.merge, resolution = 0.2,verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, reduction = "pca", dims = 1:50)
brain.merge <- RunTSNE(brain.merge, reduction = "pca", dims = 1:50)

library(harmony)
brain.merge = RunHarmony(brain.merge , group.by.vars = 'orig.ident',plot_convergence = T)

brain.merge <- FindNeighbors(brain.merge, reduction = "harmony", dims = 1:50)
brain.merge <- FindClusters(brain.merge, resolution = 0.2,verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, reduction = "harmony", dims = 1:50)
brain.merge <- RunTSNE(brain.merge, reduction = "harmony", dims = 1:50)
DimPlot(brain.merge, reduction = "umap", group.by = c("seurat_clusters", "orig.ident"),label = T,repel = T)/DimPlot(brain.merge, reduction = "tsne", group.by = c("seurat_clusters", "orig.ident"),label = T,repel = T)
saveRDS(brain.merge,"merge_harmony.rds") 
#sample merge per assay by patient#############

gc()
tempt = list()
for (i in 1:(length(brain)/2)) {
  tempt[[i]] = merge(brain[[2*i-1]],brain[[2*i]], add.cell.ids = c(brain[[2*i-1]]@meta.data$orig.ident[1],brain[[2*i]]@meta.data$orig.ident[1]), project = paste('HS',as.character(i),sep = '_'), merge.data = TRUE)
}
#old################
for (i in 1:length(brain)) {
  brain[[i]][["percent.mt"]] <- PercentageFeatureSet(brain[[i]], pattern = "^MT-")
  brain[[i]] <- subset(brain[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
  DefaultAssay(brain[[i]])='RNA'
  brain[[i]] = NormalizeData(brain[[i]])
  brain[[i]] = FindVariableFeatures(brain[[i]] , selection.method = "vst" , nfeatures = 2000)
}
features <- SelectIntegrationFeatures(object.list = brain,nfeatures = 10000,fvf.nfeatures = 2000)
brain <- lapply(X = brain, FUN = function(x) {
  DefaultAssay(x) = 'RNA'
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
brain.anchors <- FindIntegrationAnchors(object.list = brain, anchor.features = features,reduction = "rpca")
brain.merge <- IntegrateData(anchorset = brain.anchors)


#or
for (i in 1:length(brain)) {
  brain[[i]][["percent.mt"]] <- PercentageFeatureSet(brain[[i]], pattern = "^MT-")
  brain[[i]] <- subset(brain[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
  brain[[i]] <- SCTransform(brain[[i]], assay = "RNA",variable.features.n = 2000, return.only.var.genes = FALSE, verbose = FALSE)
}
memory.limit(100000000000)
gc()
features <- SelectIntegrationFeatures(object.list = brain,nfeatures = 8000)
gc()
brain = PrepSCTIntegration(object.list = brain, anchor.features = features)
gc()
brain.anchors <- FindIntegrationAnchors(object.list = brain, normalization.method = 'SCT' ,anchor.features = features)
gc()
rm(brain)
brain.merge <- IntegrateData(anchorset = brain.anchors,normalization.method = 'SCT')
gc()

#Seurat V5################
library(Seurat)
library(BPCells)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
options(future.globals.maxSize = 64 * 1000 * 1024^2)
options(Seurat.object.assay.version = "v5")

Find_doublet <- function(data){
  sweep.res.list <- paramSweep_v3(data, PCs = 1:50, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  nExp_poi <- round(0.05*ncol(data))
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  data <- doubletFinder_v3(data, PCs = 1:50, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  return(data)
}
for (i in 1:length(brain)) {
  brain[[i]][["percent.mt"]] <- PercentageFeatureSet(brain[[i]], pattern = "^MT-")
  brain[[i]] <- subset(brain[[i]], subset = nFeature_RNA > 200 & percent.mt < 10)
  brain[[i]] <- NormalizeData(brain[[i]])
  brain[[i]] <- FindVariableFeatures(brain[[i]], selection.method = "vst", nfeatures = 3000)
  brain[[i]] <- ScaleData(brain[[i]])
  brain[[i]] <- RunPCA(brain[[i]])
  brain[[i]] <- RunUMAP(brain[[i]], dims = 1:50)
  brain[[i]]<-Find_doublet(brain[[i]])
  brain[[i]]<-subset(brain[[i]],subset=doublet_info=="Singlet")
  c <- grep("pANN_",colnames(brain[[i]]@meta.data))
  brain[[i]]@meta.data <- brain[[i]]@meta.data[,-c]
  print('##########################')
  print(i)
}
list = brain
list[[1]]=NULL
brain.merge = merge(brain[[1]],y = list)
rm(list)
rm(brain)
brain.merge <- NormalizeData(brain.merge)
# split assay into nums(orig.ident) layers
brain.merge[["RNA"]] <- JoinLayers(brain.merge[["RNA"]])

brain.merge[["RNA"]] <- split(brain.merge[["RNA"]], f = brain.merge$orig.ident)
brain.merge <- FindVariableFeatures(brain.merge, verbose = FALSE)
brain.merge <- SketchData(object = brain.merge, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")
brain.merge

DefaultAssay(brain.merge) <- "sketch"
brain.merge <- FindVariableFeatures(brain.merge, verbose = F)
brain.merge <- ScaleData(brain.merge, verbose = F)
brain.merge <- RunPCA(brain.merge, verbose = F)
# integrate the datasets
brain.merge <- IntegrateLayers(brain.merge, 
                               method = RPCAIntegration, 
                               orig = "pca", 
                               new.reduction = "integrated.rpca",
                               dims = 1:30, 
                               k.anchor = 20, #reference = which(Layers(brain.merge, search = "data") %in% c("data.H_3060")),
                               verbose = F)
# cluster the integrated data
brain.merge <- FindNeighbors(brain.merge, reduction = "integrated.rpca", dims = 1:30)
brain.merge <- FindClusters(brain.merge, resolution = 2)

brain.merge <- RunUMAP(brain.merge, reduction = "integrated.rpca", dims = 1:30, return.model = T, verbose = F)
# you can now rejoin the layers in the sketched assay this is required to perform differential
# expression
brain.merge[["sketch"]] <- JoinLayers(brain.merge[["sketch"]])


#c10_markers <- FindMarkers(object = brain.merge, ident.1 = 10, max.cells.per.ident = 500, only.pos = TRUE)
##head(c10_markers)
# You can now annotate clusters using marker genes.  We performed this step, and include the
# results in the 'sketch.celltype' metadata column

plot.s1 <- DimPlot(brain.merge, group.by = "orig.ident", reduction = "umap")
plot.s2 <- DimPlot(brain.merge, group.by = "orig.ident", reduction = "umap")
plot.s1 + plot.s2 + plot_layout(ncol = 1)

# resplit the sketched cell assay into layers this is required to project the integration onto
# all cells
brain.merge[["sketch"]] <- split(brain.merge[["sketch"]], f = brain.merge$orig.ident)
brain.merge <- ProjectIntegration(object = brain.merge, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")
brain.merge <- ProjectData(object = brain.merge, sketched.assay = "sketch", assay = "RNA", 
                           sketched.reduction = "integrated.rpca.full",full.reduction = "integrated.rpca.full", #refdata = list(celltype.full = "celltype.manual", )
                           dims = 1:30)
object <- ProjectData(object = object, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
    full.reduction = "integrated.rpca.full", dims = 1:30, refdata = list(celltype.full = "celltype.manual"))
brain.merge <- RunUMAP(brain.merge, reduction = "integrated.rpca.full", dims = 1:30, reduction.name = "umap.full",reduction.key = "UMAP_full_")
p1 <- DimPlot(brain.merge, reduction = "umap.full", group.by = "orig.ident", alpha = 0.1)
p2 <- DimPlot(brain.merge, reduction = "umap.full", group.by = "orig.ident", alpha = 0.1)
p1 + p2 + plot_layout(ncol = 1)

brain.merge[["RNA"]] <- as(brain.merge[["RNA"]], "Assay")
DefaultAssay(brain.merge) = 'RNA'
saveRDS(brain.merge,'merge_v5.rds')
#large assay test##########

tampt = brain[[1]]

for (i in 2:length(brain)){
  tamptlist = c(tampt,brain[[i]])
  features <- SelectIntegrationFeatures(object.list = tamptlist,nfeatures = 8000,fvf.nfeatures = 2000)
  tamptlist <- lapply(X = tamptlist, FUN = function(x) {
    DefaultAssay(x) = 'RNA'
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  brain.anchors <- FindIntegrationAnchors(object.list = tamptlist, anchor.features = features,reduction = "rpca")
  brain.merge <- IntegrateData(anchorset = brain.anchors)
  tampt = brain.merge
}

#
brain.merge = brain[[1]]

for (i in 2:length(brain)){
  DefaultAssay(brain.merge) = 'RNA'
  brain.merge <- NormalizeData(brain.merge)
  brain.merge <- FindVariableFeatures(brain.merge, selection.method = "vst", nfeatures = 3000)
  brain.merge <- ScaleData(brain.merge)
  brain.merge <- RunPCA(brain.merge)
  brain.merge <- SCTransform(brain.merge, assay = "RNA",variable.features.n = 3000, return.only.var.genes = FALSE, verbose = FALSE)
  print(i)
  print('first')
  levels(as.factor(brain.merge@meta.data$orig.ident))
  brain[[i]] <- NormalizeData(brain[[i]])
  brain[[i]] <- FindVariableFeatures(brain[[i]], selection.method = "vst", nfeatures = 3000)
  brain[[i]] <- ScaleData(brain[[i]])
  brain[[i]] <- RunPCA(brain[[i]])
  brain[[i]] <- SCTransform(brain[[i]], assay = "RNA",variable.features.n = 3000, return.only.var.genes = FALSE, verbose = FALSE)
  tamptlist = c(brain.merge,brain[[i]])
  tamptlist
  features <- SelectIntegrationFeatures(object.list = tamptlist,nfeatures = 8000)

  tamptlist = PrepSCTIntegration(object.list = tamptlist, anchor.features = features)

  brain.anchors <- FindIntegrationAnchors(object.list = tamptlist, normalization.method = 'SCT' ,anchor.features = features)

  brain.merge <- IntegrateData(anchorset = brain.anchors,normalization.method = 'SCT')
  

}



#doubletFinder#################
Find_doublet <- function(data){
  sweep.res.list <- paramSweep_v3(data, PCs = 1:50, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  nExp_poi <- round(0.05*ncol(data))
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  data <- doubletFinder_v3(data, PCs = 1:50, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  return(data)
}
for (i in 1:length(brain)) {
  brain[[i]][["percent.mt"]] <- PercentageFeatureSet(brain[[i]], pattern = "^MT-")
  brain[[i]] <- subset(brain[[i]], subset = nFeature_RNA > 200 & percent.mt < 10)
  brain[[i]] <- NormalizeData(brain[[i]])
  brain[[i]] <- FindVariableFeatures(brain[[i]], selection.method = "vst", nfeatures = 3000)
  brain[[i]] <- ScaleData(brain[[i]])
  brain[[i]] <- RunPCA(brain[[i]])
  brain[[i]] <- RunUMAP(brain[[i]], dims = 1:50)
  brain[[i]]<-Find_doublet(brain[[i]])
  #brain[[i]]<-subset(brain[[i]],subset=doublet_info=="Singlet")
  c <- grep("pANN_",colnames(brain[[i]]@meta.data))
  brain[[i]]@meta.data <- brain[[i]]@meta.data[,-c]
  brain[[i]] <- SCTransform(brain[[i]], assay = "RNA",variable.features.n = 2000, return.only.var.genes = FALSE, verbose = FALSE)
  cat(i);
  #DefaultAssay(brain[[i]]) = 'RNA'
}
memory.limit(100000000000)
gc()
features <- SelectIntegrationFeatures(object.list = brain,nfeatures = 8000)
gc()
brain = PrepSCTIntegration(object.list = brain, anchor.features = features)
gc()
brain.anchors <- FindIntegrationAnchors(object.list = brain, normalization.method = 'SCT' ,anchor.features = features)
gc()
#rm(brain)
brain.merge <- IntegrateData(anchorset = brain.anchors,normalization.method = 'SCT')
gc()




#tampt ancher save and load########
brain.anchors = readRDS('anchors.rds')
saveRDS(brain.anchors,'anchors.rds')

brain.anchors <- FindIntegrationAnchors(object.list = brain,anchor.features = features,dims=1:50)
brain.merge <- IntegrateData(anchorset = brain.anchors,dims=1:50)
###########

dim(brain.merge)
DefaultAssay(brain.merge) <- "integrated"
brain.merge <- ScaleData(brain.merge, verbose = FALSE)
brain.merge <- RunPCA(brain.merge, npcs = 50, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, reduction = "pca", dims = 1:50)
brain.merge <- FindClusters(brain.merge, resolution = 1.0,verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, reduction = "pca", dims = 1:50)
brain.merge <- RunTSNE(brain.merge, reduction = "pca", dims = 1:50)
pix = DimPlot(brain.merge, raster=FALSE,reduction = "umap", group.by = c("seurat_clusters", "orig.ident"),label = T,repel = T)/DimPlot(brain.merge,raster=FALSE, reduction = "tsne", group.by = c("seurat_clusters", "orig.ident"),label = T,repel = T)
ggsave('merge_merge1.png',pix,width=15,height=15)
saveRDS(brain.merge,"merge_merge.rds") 


group_by = c('seurat_clusters','orig.ident','Group','No.','manucelltype','manusubcelltype')
pix = DimPlot(brain.merge, reduction = "umap",group.by = group_by, label = TRUE ,repel = TRUE,ncol = length(group_by))/DimPlot(brain.merge, reduction = "tsne",group.by = group_by, label = TRUE ,repel = TRUE,ncol = length(group_by))
ggsave('merge_merge2.png',pix,width=40,height=10)


brain.markers <- FindAllMarkers(brain.merge, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
brain.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
write.csv(top10,"merge_diff.csv")

list = SplitObject(list[[2]], split.by = 'orig.ident')

DoHeatmap(brain.merge, features = top10$gene)




