# Code for documenting the datasets


#' Sample Data
#'
#' A dataset containing subsample from a single tenplex
#'
#' @format A single dataframe containing the
#' same 14 columns.  The first three rows of each dataframe constitute a data
#' header that will be used to specify the correct model to use in the data
#' analysis.  The 14 columns must all be present and must have the names seen in
#' this example.
#' \describe{
#'  \item{Protein}{A string that identifies the protein or protein group for
#'  which relative estimates will be computed.  The first three rows of this
#'  column must have 0 or 1 entries with a 1 denoting the use of the relevant
#'  header row.  The first row is used to specify conditions for the signal to
#'  noise columns.  Since conditions must always be present, it does not matter
#'  if a 0 or 1 is used here.  The indicator in the second row lets us know if
#'  biological replicates are to be defined in the second row (they might
#'  alternatively be defined in column 3).  The third and final header row is
#'  used to identify post-translational modifications.  If these are used a 1 must
#'  be entered into row 3 column 1. }
#'  \item{Peptide}{Unique peptide identifier.  Only used in PTM analyses but the
#'  column must be present.  Header entries here are irrelevant.}
#'  \item{bioID}{This column may be used to denote biological replicates.  Either
#'  strings or numbers can be used, to key point being that any entries with the
#'  same value will be treated as biological replicates.  Depending on the data
#'  source, specifying this in a column header may not be practical, which is why
#'  the option to specify biological replicates within each row has been made
#'  available.  Entries to the header rows in this column are irrelevant with the
#'  exception of Row 1 which indicates whether or not the column is being used.}
#'  \item{Covariate}{This column provides a continuous covariate to be used in a
#'  compositional non-linear regression.  Typically these entries will either be
#'  summed signal-to-noise values or isolation specificities.  If a covariate is
#'  to be used, a one must be entered in the first row of this column.  The
#'  second and third rows do nothing.}
#'  \item{varCat}{This column provides a way to allow for additional variance
#'  components.  Separate variance components are always generated for each
#'  plex/tag combination.  However, it might be of interest to further
#'  separate variance accoriding to some other factor that varies within a
#'  tag, e.g. species.  This is accomplished by adding a categorical variable
#'  for variance groups.  Observations beloning to the same category should
#'  have the same category entry (either an integer or string).  These
#'  categories are only intended for protein analysis and have no bearing on
#'  PTM's (beyond altering the protein results). A 1 must be entered into the
#'  first row of theheader when using this column.}
#'  \item{tag1-tagN}{These columns contain signal to noise intensity ratios.
#'  The column names are case sensitive and must include the string 'tag'
#'  followed by a number.  The entries in the first three rows of these columns
#'  are used to tell the program which columns represent various conditions,
#'  biological replicates, and PTMs (not currently in use).  Columns with the
#'  same entry in row 1 are considered to be from the same condition (this
#'  determines what comparisons will be made).  Similarly equivalent values in
#'  row two will determine which columns are biological replicates. The same
#'  value in the third column denotes the same type of PTM, e.g. 1 =
#'  phosphorylation and 2 = ubiqutination. Ten columns are not necessary as the program accepts
#'  an arbitrary number of tag columns.  Values are assumed to represent signal
#'  to noise measurements.  Accordingly any values that are missing or less than
#'  one will be replaced with 1's during data processing.}
#' }
#'
#'
"sampleDat"
