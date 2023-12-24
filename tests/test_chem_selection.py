from deemian.chem.selection import mol_dataframe_selection


def test_mol_dataframe_selection(n1_df, spike_df, vps4_df):
    amino_acids = [
        "ALA",
        "CYS",
        "ASP",
        "GLU",
        "PHE",
        "GLY",
        "HIS",
        "ILE",
        "LYS",
        "LEU",
        "MET",
        "ASN",
        "PRO",
        "GLN",
        "ARG",
        "SER",
        "THR",
        "VAL",
        "TRP",
        "TYR",
        "PYL",
        "SEC",
        "NVA",
    ]

    protein_A = [("chain", "A"), ("and", "protein")]
    oseltamivir = [("chain", "A"), ("and", "resname", "G39")]

    ace2_hotspot31 = [("chain", "A"), ("and", "protein"), ("and", "resid_range", "1", "55")]
    ace2_hotspot353 = [("chain", "A"), ("and", "protein"), ("and", "resid_range", "341", "364")]
    spike_rbm = [("chain", "E"), ("and", "protein"), ("and", "resid_range", "470", "510")]

    vps4 = [("protein",), ("and", "chain", "A")]
    vps4_20_to_23 = [("protein",), ("and", "chain", "A"), ("and", "resid", "20", "21", "22", "23")]
    chmp6 = [("protein",), ("and", "chain", "B")]
    chmp6_too = [("protein",), ("and", "not", "chain", "A")]

    protein_A_result = mol_dataframe_selection(protein_A, n1_df)
    oseltamivir_result = mol_dataframe_selection(oseltamivir, n1_df)

    ace2_hotspot31_result = mol_dataframe_selection(ace2_hotspot31, spike_df)
    ace2_hostpot353_result = mol_dataframe_selection(ace2_hotspot353, spike_df)
    spike_rbm_result = mol_dataframe_selection(spike_rbm, spike_df)

    vps4_result = mol_dataframe_selection(vps4, vps4_df)
    vps4_20_to_23_result = mol_dataframe_selection(vps4_20_to_23, vps4_df)
    chmp6_result = mol_dataframe_selection(chmp6, vps4_df)
    chmp6_too_result = mol_dataframe_selection(chmp6_too, vps4_df)

    protein_A_expected = n1_df[(n1_df["chain_id"] == "A") & (n1_df["residue_name"].isin(amino_acids))]
    oseltamivir_expected = n1_df[(n1_df["chain_id"] == "A") & (n1_df["residue_name"] == "G39")]

    h31_id_range = range(1, 56)
    h353_id_range = range(341, 365)
    rbm_id_range = range(470, 511)
    ace2_hotspot31_expected = spike_df[
        (spike_df["chain_id"] == "A")
        & (spike_df["residue_name"].isin(amino_acids))
        & (spike_df["residue_number"].isin(h31_id_range))
    ]
    ace2_hotspot353_expected = spike_df[
        (spike_df["chain_id"] == "A")
        & (spike_df["residue_name"].isin(amino_acids))
        & (spike_df["residue_number"].isin(h353_id_range))
    ]
    spike_rbm_expected = spike_df[
        (spike_df["chain_id"] == "E")
        & (spike_df["residue_name"].isin(amino_acids))
        & (spike_df["residue_number"].isin(rbm_id_range))
    ]

    id_20_to23 = [20, 21, 22, 23]
    vps4_expected = vps4_df[(vps4_df["chain_id"] == "A") & (vps4_df["residue_name"].isin(amino_acids))]
    vps4_20_to_23_expected = vps4_df[
        (vps4_df["chain_id"] == "A")
        & (vps4_df["residue_name"].isin(amino_acids))
        & (vps4_df["residue_number"].isin(id_20_to23))
    ]
    chmp6_expected = vps4_df[(vps4_df["chain_id"] == "B") & (vps4_df["residue_name"].isin(amino_acids))]
    chmp6_too_expected = vps4_df[(vps4_df["chain_id"] != "A") & (vps4_df["residue_name"].isin(amino_acids))]

    assert protein_A_result.equals(protein_A_expected)
    assert oseltamivir_result.equals(oseltamivir_expected)

    assert ace2_hotspot31_result.equals(ace2_hotspot31_expected)
    assert ace2_hostpot353_result.equals(ace2_hotspot353_expected)
    assert spike_rbm_result.equals(spike_rbm_expected)

    assert vps4_result.equals(vps4_expected)
    assert vps4_20_to_23_result.equals(vps4_20_to_23_expected)
    assert chmp6_result.equals(chmp6_expected)
    assert chmp6_too_result.equals(chmp6_too_expected)
