def samplots(samparse_tsv):

    print("Plotting...")

    samparse_df = pd.read_csv(samparse_tsv, sep="\t")

    # remove unaligned query names
    samparse_df = samparse_df[samparse_df["reference_name"] != "*"]

    # log log length plot
    samparse_df["log_aln_length"] = np.log10(samparse_df["aln_length"])
    pyplot.figure(figsize=(12, 12))
    pyplot.hist(samparse_df["log_aln_length"], log=True, histtype="step")
    pyplot.xlabel("log aln length")
    pyplot.ylabel("log count")
    pyplot.savefig("log_aln_length.png")

    # NM plot
    pyplot.figure(figsize=(12, 12))
    pyplot.hist(
        samparse_df[["reference_name", "NM"]].groupby("reference_name").mean()["NM"],
        log=True,
        histtype="step",
    )
    pyplot.xlabel("Mean mismatches per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("NM_per_contig_mean.png")

    # M plot
    pyplot.figure(figsize=(12, 12))
    pyplot.hist(
        samparse_df[["reference_name", "M"]].groupby("reference_name").mean()["M"],
        log=True,
        histtype="step",
    )
    pyplot.xlabel("Mean matches per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("M_per_contig_mean.png")

    # i plot
    pyplot.figure(figsize=(12, 12))
    pyplot.hist(
        samparse_df[["reference_name", "I"]].groupby("reference_name").mean()["I"],
        log=True,
        histtype="step",
    )
    pyplot.xlabel("Mean insertions per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("I_per_contig_mean.png")

    # d plot
    pyplot.figure(figsize=(12, 12))
    pyplot.hist(
        samparse_df[["reference_name", "D"]].groupby("reference_name").mean()["D"],
        log=True,
        histtype="step",
    )
    pyplot.xlabel("Mean deletions per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("D_per_contig_mean.png")

    # N mean plot
    pyplot.figure(figsize=(12, 12))
    d_max_contig = pyplot.hist(
        samparse_df[["reference_name", "N"]].groupby("reference_name").mean()["N"],
        log=True,
        histtype="step",
    )
    pyplot.xlabel("Mean intron length per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("N_per_contig_mean.png")

    # N max plot
    pyplot.figure(figsize=(12, 12))
    d_max_contig = pyplot.hist(
        samparse_df[["reference_name", "N"]].groupby("reference_name").max()["N"],
        log=True,
        histtype="step",
    )
    pyplot.xlabel("Max intron length per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("N_per_contig_max.png")

    # mean identity per contig plot
    pyplot.figure(figsize=(12, 12))
    identity_contig = pyplot.hist(
        samparse_df[["reference_name", "identity"]].groupby("reference_name").mean()["identity"],
        log=True,
        histtype="step",
    )
    pyplot.xlabel("Mean identity per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("identity_per_contig_mean.png")

    # mean query coverage per contig plot
    pyplot.figure(figsize=(12, 12))
    query_coverage_contig = pyplot.hist(
        samparse_df[["reference_name", "query_coverage"]].groupby("reference_name").mean()["query_coverage"],
        log=True,
        histtype="step",
    )
    pyplot.xlabel("Mean query coverage per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("qcov_per_contig_mean.png")

    print("Done")

