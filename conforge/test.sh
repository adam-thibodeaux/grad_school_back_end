/Applications/ccp4-9/bin/phaser << eof
    MODE MR_AUTO
    HKLIN /Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha.mtz
    ENSEMBLE par PDBFILE /Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha_frag_4_correct_conformation.pdb IDENTITY 0.44
    ENSEMBLE par HETATM ON
    FORMFACTORS ELECTRON
    ELLG TARGET 100
    PACK CUTOFF 100
    PACK KEEP HIGH TFZ ON
    PURGE ROT ENABLE OFF
    PURGE TRA ENABLE OFF
    PURGE RNP ENABLE OFF
    SGALTERNATIVE SELECT ALL
    XYZOUT ON ENSEMBLE ON
    XYZOUT ON PACKING ON 
    TOPFILES 5
    ZSCORE USE OFF
    ZSCORE SOLVED 2
    ROOT /Users/adam/Downloads/outputs_from_molec_replac/PAR_FRAG/PAR_FRAG_4/pariteprevir_frag_4
    SEARCH ENSEMBLE par
    eof