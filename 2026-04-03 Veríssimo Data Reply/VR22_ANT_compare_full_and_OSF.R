# This script compares the full dataset sent to Holcombe & Tuladhar (on 3 April 2026) ...
# ... to the ANT-data.csv publicly available at the OSF at https://osf.io/59er2/files/osfstorage
# João Veríssimo
# 3 April 2026

# Read data
ant <- read.table("2026-04-03 Veríssimo Data Reply/VR22_ANT_full_dataset_for_Holcombe_Tuladhar.csv", sep = ",", header = TRUE, fileEncoding = "utf-8")
head(ant)

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

# Age breakdown of participants with lower than 50% Accuracy
ant_under50 <- ant[ant$AccuracyAtLeast50 == 0,] |> droplevels()
length(unique(ant_under50$Participant_ID))
ant_under50$Participant_ID <- gsub("New", "", ant_under50$Participant_ID) 
ant_under50$Participant_ID <- as.factor(ant_under50$Participant_ID) |> as.numeric()
age_break <- ant_under50[, c("Participant_ID", "Age")] |>
  distinct(Participant_ID, .keep_all = TRUE)
age_break |>
  mutate(
    age_group = cut(
      Age,
      breaks = c(0, 60, 65, 70, 75, 80, 85, Inf),
      labels = c("<60", "60-65", "65-70", "70-75", "75-80", "80-85", "85<"),
      right = FALSE
    )
  ) |>
  count(age_group, sort = TRUE) |>
  mutate(pct = round(n / sum(n) * 100, 1))


# Exclude timeouts, wrong button presses, and short outliers (< 100ms)
ant[ant$Timeout == 1 | ant$Wrong_Button_Press == 1 | ant$Short_Outlier == 1,] |> nrow()
ant <- ant[ant$Timeout == 0 & ant$Wrong_Button_Press == 0 & ant$Short_Outlier == 0,]
nrow(ant)

# Read previously available OSF data
head(ant_osf <- read.table("2026-04-03 Veríssimo Data Reply/VR22_ANT_OSF_data.csv", sep=",", header=T, fileEncoding="utf-8"))
length(unique(ant_osf$Participant_ID))

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

osf_age <- ant_osf[ant_osf$AccuracyAtLeast75 == 0,] |> droplevels()
osf_age <- osf_age[, c("Participant_ID", "Age")] |> distinct(Participant_ID, .keep_all = TRUE)
osf_age |>
  mutate(
    age_group = cut(
      Age,
      breaks = c(0, 60, 65, 70, 75, 80, 85, Inf),
      labels = c("<60", "60-65", "65-70", "70-75", "75-80", "80-85", "85<"),
      right = FALSE
    )
  ) |>
  count(age_group, sort = TRUE) |>
  mutate(pct = round(n / sum(n) * 100, 1))

