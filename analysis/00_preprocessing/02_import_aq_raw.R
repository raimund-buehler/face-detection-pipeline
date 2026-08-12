## Manuscript section: preprocessing
## Analysis family: AQ questionnaire import and scoring
## Original source path: legacy/data_archive/original_layout/data_questionnaire/AQ.r
## Primary input dataset(s): legacy/data_archive/original_layout/data_questionnaire/AQ.csv, legacy/data_archive/original_layout/data_questionnaire/scoring_key.csv
## Primary output(s): data/derived/questionnaire_scoring/AQ_final_data.csv
## Known TODOs: script remains largely provenance-oriented; import logic is preserved and not yet fully standardized
## Scientific logic unchanged from source except for canonical file-path cleanup

# This script reads a CSV file in GNU R.
# While reading this file, comments will be created for all variables.
# The comments for values will be stored as attributes (attr) as well.

df_file <- "legacy/data_archive/original_layout/data_questionnaire/AQ.csv"

options(encoding = "UTF-8")
df <- read.delim(
    file = df_file, encoding = "UTF-8", fileEncoding = "UTF-8",
    header = FALSE, sep = "\t", quote = "\"",
    dec = ".", row.names = NULL,
    col.names = c(
        "CASE", "SERIAL", "REF", "QUESTNNR", "MODE", "STARTED", "AQ01", "AQ02", "AQ03", "AQ04",
        "AQ05", "AQ06", "AQ07", "AQ08", "AQ09", "AQ10", "AQ11", "AQ12", "AQ13", "AQ14", "AQ15",
        "AQ16", "AQ17", "AQ18", "AQ19", "AQ20", "AQ21", "AQ22", "AQ23", "AQ24", "AQ25", "AQ26",
        "AQ27", "AQ28", "AQ29", "AQ30", "AQ31", "AQ32", "AQ33", "SF29_01", "SF73", "SF74", "SF75",
        "SF76", "TIME001", "TIME002", "TIME003", "TIME_SUM", "MAILSENT", "LASTDATA",
        "FINISHED", "Q_VIEWER", "LASTPAGE", "MAXPAGE", "MISSING", "MISSREL", "TIME_RSI"
    ),
    as.is = TRUE,
    colClasses = c(
        CASE = "numeric", SERIAL = "character", REF = "character", QUESTNNR = "character",
        MODE = "factor", STARTED = "POSIXct", AQ01 = "numeric", AQ02 = "numeric",
        AQ03 = "numeric", AQ04 = "numeric", AQ05 = "numeric", AQ06 = "numeric",
        AQ07 = "numeric", AQ08 = "numeric", AQ09 = "numeric", AQ10 = "numeric",
        AQ11 = "numeric", AQ12 = "numeric", AQ13 = "numeric", AQ14 = "numeric",
        AQ15 = "numeric", AQ16 = "numeric", AQ17 = "numeric", AQ18 = "numeric",
        AQ19 = "numeric", AQ20 = "numeric", AQ21 = "numeric", AQ22 = "numeric",
        AQ23 = "numeric", AQ24 = "numeric", AQ25 = "numeric", AQ26 = "numeric",
        AQ27 = "numeric", AQ28 = "numeric", AQ29 = "numeric", AQ30 = "numeric",
        AQ31 = "numeric", AQ32 = "numeric", AQ33 = "numeric", SF29_01 = "character",
        SF73 = "numeric", SF74 = "numeric", SF75 = "numeric", SF76 = "numeric",
        TIME001 = "integer", TIME002 = "integer", TIME003 = "integer", TIME_SUM = "integer",
        MAILSENT = "POSIXct", LASTDATA = "POSIXct", FINISHED = "logical",
        Q_VIEWER = "logical", LASTPAGE = "numeric", MAXPAGE = "numeric",
        MISSING = "numeric", MISSREL = "numeric", TIME_RSI = "numeric"
    ),
    skip = 1,
    check.names = TRUE, fill = TRUE,
    strip.white = FALSE, blank.lines.skip = TRUE,
    comment.char = "",
    na.strings = ""
)

row.names(df) <- df$CASE

rm(df_file)

attr(df, "project") <- "oxy_sept"
attr(df, "description") <- "Oxy_Sept"
attr(df, "date") <- "2024-09-06 13:32:22"
attr(df, "server") <- "https://www.soscisurvey.de"

# Variable und Value Labels
df$AQ01 <- factor(df$AQ01,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ02 <- factor(df$AQ02,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ03 <- factor(df$AQ03,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ04 <- factor(df$AQ04,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ05 <- factor(df$AQ05,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ06 <- factor(df$AQ06,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ07 <- factor(df$AQ07,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ08 <- factor(df$AQ08,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ09 <- factor(df$AQ09,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ10 <- factor(df$AQ10,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ11 <- factor(df$AQ11,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ12 <- factor(df$AQ12,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ13 <- factor(df$AQ13,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ14 <- factor(df$AQ14,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ15 <- factor(df$AQ15,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ16 <- factor(df$AQ16,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ17 <- factor(df$AQ17,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ18 <- factor(df$AQ18,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ19 <- factor(df$AQ19,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ20 <- factor(df$AQ20,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ21 <- factor(df$AQ21,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ22 <- factor(df$AQ22,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ23 <- factor(df$AQ23,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ24 <- factor(df$AQ24,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ25 <- factor(df$AQ25,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ26 <- factor(df$AQ26,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ27 <- factor(df$AQ27,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ28 <- factor(df$AQ28,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ29 <- factor(df$AQ29,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ30 <- factor(df$AQ30,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ31 <- factor(df$AQ31,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ32 <- factor(df$AQ32,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$AQ33 <- factor(df$AQ33,
    levels = c("1", "2", "3", "4", "-9"),
    labels = c("ich stimme eindeutig zu", "ich stimme ein wenig zu", "ich stimme eher nicht zu", "ich stimme überhaupt nicht zu", "[NA] Not answered"),
    ordered = FALSE
)
df$SF73 <- factor(df$SF73,
    levels = c("1", "2", "3", "4", "5", "6", "7", "-9"),
    labels = c("überhaupt nicht [1]", "[2]", "[3]", "[4]", "[5]", "[6]", "sehr [7]", "[NA] Not answered"),
    ordered = FALSE
)
df$SF74 <- factor(df$SF74,
    levels = c("1", "2", "3", "4", "5", "6", "7", "-9"),
    labels = c("überhaupt nicht [1]", "[2]", "[3]", "[4]", "[5]", "[6]", "sehr [7]", "[NA] Not answered"),
    ordered = FALSE
)
df$SF75 <- factor(df$SF75,
    levels = c("1", "2", "3", "4", "5", "6", "7", "-9"),
    labels = c("überhaupt nicht [1]", "[2]", "[3]", "[4]", "[5]", "[6]", "sehr [7]", "[NA] Not answered"),
    ordered = FALSE
)
df$SF76 <- factor(df$SF76,
    levels = c("1", "2", "3", "4", "5", "6", "7", "-9"),
    labels = c("überhaupt nicht [1]", "[2]", "[3]", "[4]", "[5]", "[6]", "sehr [7]", "[NA] Not answered"),
    ordered = FALSE
)
attr(df$FINISHED, "F") <- "Canceled"
attr(df$FINISHED, "T") <- "Finished"
attr(df$Q_VIEWER, "F") <- "Respondent"
attr(df$Q_VIEWER, "T") <- "Spectator"
comment(df$SERIAL) <- "Serial number (if provided)"
comment(df$REF) <- "Reference (if provided in link)"
comment(df$QUESTNNR) <- "Questionnaire that has been used in the interview"
comment(df$MODE) <- "Interview mode"
comment(df$STARTED) <- "Time the interview has started (Europe/Berlin)"
comment(df$AQ01) <- "1 (AQ-K)"
comment(df$AQ02) <- "2 (AQ-K)"
comment(df$AQ03) <- "3 (AQ-K)"
comment(df$AQ04) <- "4 (AQ-K)"
comment(df$AQ05) <- "5 (AQ-K)"
comment(df$AQ06) <- "6 (AQ-K)"
comment(df$AQ07) <- "7 (AQ-K)"
comment(df$AQ08) <- "8 (AQ-K)"
comment(df$AQ09) <- "9 (AQ-K)"
comment(df$AQ10) <- "10 (AQ-K)"
comment(df$AQ11) <- "11 (AQ-K)"
comment(df$AQ12) <- "12 (AQ-K)"
comment(df$AQ13) <- "13 (AQ-K)"
comment(df$AQ14) <- "14 (AQ-K)"
comment(df$AQ15) <- "15 (AQ-K)"
comment(df$AQ16) <- "16 (AQ-K)"
comment(df$AQ17) <- "17 (AQ-K)"
comment(df$AQ18) <- "18 (AQ-K)"
comment(df$AQ19) <- "19 (AQ-K)"
comment(df$AQ20) <- "20 (AQ-K)"
comment(df$AQ21) <- "21 (AQ-K)"
comment(df$AQ22) <- "22 (AQ-K)"
comment(df$AQ23) <- "23 (AQ-K)"
comment(df$AQ24) <- "24 (AQ-K)"
comment(df$AQ25) <- "25 (AQ-K)"
comment(df$AQ26) <- "26 (AQ-K)"
comment(df$AQ27) <- "27 (AQ-K)"
comment(df$AQ28) <- "28 (AQ-K)"
comment(df$AQ29) <- "29 (AQ-K)"
comment(df$AQ30) <- "30 (AQ-K)"
comment(df$AQ31) <- "31 (AQ-K)"
comment(df$AQ32) <- "32 (AQ-K)"
comment(df$AQ33) <- "33 (AQ-K)"
comment(df$SF29_01) <- "Codeeingabe: Code"
comment(df$SF73) <- "SYM1"
comment(df$SF74) <- "SYM2"
comment(df$SF75) <- "SYM3"
comment(df$SF76) <- "SYM4"
comment(df$TIME001) <- "Time spent on page 1"
comment(df$TIME002) <- "Time spent on page 2"
comment(df$TIME003) <- "Time spent on page 3"
comment(df$TIME_SUM) <- "Time spent overall (except outliers)"
comment(df$MAILSENT) <- "Time when the invitation mailing was sent (personally identifiable recipients, only)"
comment(df$LASTDATA) <- "Time when the data was most recently updated"
comment(df$FINISHED) <- "Has the interview been finished (reached last page)?"
comment(df$Q_VIEWER) <- "Did the respondent only view the questionnaire, omitting mandatory questions?"
comment(df$LASTPAGE) <- "Last page that the participant has handled in the questionnaire"
comment(df$MAXPAGE) <- "Hindmost page handled by the participant"
comment(df$MISSING) <- "Missing answers in percent"
comment(df$MISSREL) <- "Missing answers (weighted by relevance)"
comment(df$TIME_RSI) <- "Completion Speed (relative)"



# Assure that the comments are retained in subsets
as.data.frame.avector <- as.data.frame.vector
`[.avector` <- function(x, i, ...) {
    r <- NextMethod("[")
    mostattributes(r) <- attributes(x)
    r
}
df_tmp <- data.frame(
    lapply(df, function(x) {
        structure(x, class = c("avector", class(x)))
    })
)
mostattributes(df_tmp) <- attributes(df)
df <- df_tmp
rm(df_tmp)

# SozSpont = 1,7,8,10,11,13,14,20,24,28,31
# Imag = 3,5,6,9,16,17,18,22,23,26,32,33
# CommRec = 2,4,12,15,19,21,25,27,29,30
library(tidyverse)


scoring_key <- read.csv("legacy/data_archive/original_layout/data_questionnaire/scoring_key.csv")

# Subset only the AQ response columns and the ID column (SF29_01)
aq_data <- df %>%
    select(SF29_01, starts_with("AQ")) # Include the ID column

# Convert the responses into numeric categories (4 = "ich stimme eindeutig zu", etc.)
# Ensure that the responses are recoded into 1 to 4 based on their meaning
recoded_aq_data <- aq_data %>%
    mutate_at(vars(starts_with("AQ")), ~ recode(.,
        "ich stimme eindeutig zu" = 1,
        "ich stimme ein wenig zu" = 2,
        "ich stimme eher nicht zu" = 3,
        "ich stimme überhaupt nicht zu" = 4
    ))

# Reshape the scoring key for easier mapping
# The scoring_key has columns that correspond to the recoded response numbers, so we will map these
scoring_key_long <- scoring_key %>%
    gather(key = "Response", value = "Score", -X)

# Apply the scoring key
for (i in 1:nrow(scoring_key)) {
    item <- scoring_key[i, 1]
    scores <- as.numeric(scoring_key[i, -1]) # Convert the scores to a numeric vector

    recoded_aq_data[[item]] <- recode(recoded_aq_data[[item]],
        "4" = scores[1],
        "3" = scores[2],
        "2" = scores[3],
        "1" = scores[4]
    )
}

# Create subscales based on the provided item mapping
subscales <- recoded_aq_data %>%
    mutate(
        SozSpont = AQ01 + AQ07 + AQ08 + AQ10 + AQ11 + AQ13 + AQ14 + AQ20 + AQ24 + AQ28 + AQ31,
        Imag = AQ03 + AQ05 + AQ06 + AQ09 + AQ16 + AQ17 + AQ18 + AQ22 + AQ23 + AQ26 + AQ32 + AQ33,
        CommRec = AQ02 + AQ04 + AQ12 + AQ15 + AQ19 + AQ21 + AQ25 + AQ27 + AQ29 + AQ30
    )

subscales <- subscales %>%
    mutate(OverallAQ = rowSums(.[2:34], na.rm = TRUE))

# Filter out rows with NA values in any of the subscales (SozSpont, Imag, CommRec)
filtered_data <- subscales %>%
    filter(!is.na(SozSpont) & !is.na(Imag) & !is.na(CommRec))

# Select the ID and the calculated subscales and scores
final_data <- filtered_data %>%
    select(SF29_01, SozSpont, Imag, CommRec, OverallAQ) %>%
    rename(ID = SF29_01) %>%
    arrange(ID) %>%
    mutate(Soz_Scales = SozSpont + CommRec)

# View the result with overall score and subscales
print(subscales)

# Save the resulting data to a new CSV file
write.csv(final_data, "data/derived/questionnaire_scoring/AQ_final_data.csv", row.names = FALSE)
