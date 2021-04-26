#Import mapq values
inputMapQUV1 <- read.table(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L2_Pool_1_PA_UV/mapq.txt")
inputMapQUV2 <- read.table(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L3_Pool_2_PA_UV/mapq.txt")
inputMapQUV3 <- read.table(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L4_Pool_3_PA_UV/mapq.txt")
inputMapQVIS1 <- read.table(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L2_Pool_1_PA_VIS/mapq.txt")
inputMapQVIS2 <- read.table(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L3_Pool_2_PA_VIS/mapq.txt")
inputMapQVIS3 <- read.table(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L4_Pool_3_PA_VIS/mapq.txt")

#Count the number of reads with a MapQ value of 60
length(which(inputMapQUV1$V1 == "60"))
length(which(inputMapQUV2$V1 == "60"))
length(which(inputMapQUV3$V1 == "60"))
length(which(inputMapQVIS1$V1 == "60"))
length(which(inputMapQVIS2$V1 == "60"))
length(which(inputMapQVIS3$V1 == "60"))

#Count the number of reads with a MapQ value of 1
length(which(inputMapQUV1$V1 == "1"))
length(which(inputMapQUV2$V1 == "1"))
length(which(inputMapQUV3$V1 == "1"))
length(which(inputMapQVIS1$V1 == "1"))
length(which(inputMapQVIS2$V1 == "1"))
length(which(inputMapQVIS3$V1 == "1"))

#Count the number of reads with a MapQ value of 0
length(which(inputMapQUV1$V1 == "0"))
length(which(inputMapQUV2$V1 == "0"))
length(which(inputMapQUV3$V1 == "0"))
length(which(inputMapQVIS1$V1 == "0"))
length(which(inputMapQVIS2$V1 == "0"))
length(which(inputMapQVIS3$V1 == "0"))
