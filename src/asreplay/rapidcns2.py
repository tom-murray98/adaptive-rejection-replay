"""Apply the experiment-local modification-code fix to SCARLET RapidCNS2."""
FILTER = '''# Keep methylation observations before dropping the modification-code column.
# ROBIN also filters bedMethyl to mod_code == "m". The experiment BAM has
# already converted h to m, so its residual h pileup rows are not observations
# of additional methylation and must not compete with m at the same CpG.
if (ncol(meth) < 11L) stop("Expected bedMethyl with at least 11 columns")
message("Input modification codes: ", paste(names(table(meth[[4]])),
        as.integer(table(meth[[4]])), collapse = "; "))
meth <- meth[as.character(meth[[4]]) == "m", , drop = FALSE]
if (nrow(meth) == 0L) stop("No methylation (m) rows in bedMethyl input")
message("Methylation rows retained: ", nrow(meth))
'''


def repair(text):
    if '# Keep methylation observations before dropping' in text:
        return text
    anchor = 'meth <- read.delim(opt$in_file,header=FALSE)'
    if text.count(anchor) != 1:
        raise ValueError('Expected one known SCARLET bedMethyl input statement')
    return text.replace(anchor, anchor + '\n\n' + FILTER)


if __name__ == '__main__':
    import sys
    from pathlib import Path
    source, dest = map(Path, sys.argv[1:])
    dest.write_text(repair(source.read_text()))
