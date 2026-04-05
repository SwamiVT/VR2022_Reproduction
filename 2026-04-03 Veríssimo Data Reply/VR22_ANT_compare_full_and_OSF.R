# This script compares the full dataset sent to Holcombe & Tuladhar (on 3 April 2026) ...
# ... to the ANT-data.csv publicly available at the OSF at https://osf.io/59er2/files/osfstorage
# João Veríssimo
# 3 April 2026

# Read data
head(ant <- read.table("VR22_ANT_full_dataset_for_Holcombe_Tuladhar.csv", sep=",", header=T, fileEncoding="utf-8"))

# How many participants and trials?
length(unique(ant$Participant_ID))
nrow(ant)

# Exclude participants with disorders
ant[ant$Disorders == 1, "Participant_ID"] |> unique() |> length()
ant <- ant[ant$Disorders == 0,] |> droplevels()
length(unique(ant$Participant_ID))

# Exclude participants with accuracy lower than 50%
ant[ant$AccuracyAtLeast50 == 0, "Participant_ID"] |> unique() |> length()
ant <- ant[ant$AccuracyAtLeast50 == 1,] |> droplevels()
length(unique(ant$Participant_ID))
nrow(ant)

# Exclude timeouts, wrong button presses, and short outliers (< 100ms)
ant[ant$Timeout == 1 | ant$Wrong_Button_Press == 1 | ant$Short_Outlier == 1,] |> nrow()
ant <- ant[ant$Timeout == 0 & ant$Wrong_Button_Press == 0 & ant$Short_Outlier == 0,]
nrow(ant)

# Read previously available OSF data
head(ant_osf <- read.table("VR22_ANT_OSF_data.csv", sep=",", header=T, fileEncoding="utf-8"))

# Participant numbers anonymised to match identifiers in OSF file
ant$Participant_ID <- gsub("New", "", ant$Participant_ID) |> as.numeric()
ant$Participant_ID <- as.factor(ant$Participant_ID) |> as.numeric()

# Order as in the OSF file (by participant x trial combination)
order_osf   <- interaction(ant_osf$Participant_ID, ant_osf$Trial, drop = TRUE)
order_new   <- interaction(ant$Participant_ID, ant$Trial, drop = TRUE)
ant <- ant[match(order_osf, order_new), ]
rownames(ant) <- NULL

# Keep only those columns in OSF file
ant <- ant[, colnames(ant_osf)]

# All equal?
all.equal(ant, ant_osf)
