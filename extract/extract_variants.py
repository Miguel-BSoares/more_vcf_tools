def parse_coords(region):
    # Updated region calculation
    start = max(1, pos + int(left_flank) - 1)  # start = pos + left_flank - 1
    end = pos + int(right_flank)

    return start, end

# ... Other code ... 

for variant in variants:
    # Updated AF handling logic
    afs = variant.INFO.get('AF')
    if afs is not None:
        if isinstance(afs, (list, tuple)):
            af_val = afs[0]
        else:
            af_val = afs
        try:
            allele_freq_alt = float(af_val)
            allele_freq_ref = 1 - allele_freq_alt
        except Exception:
            allele_freq_alt = 'NA'
            allele_freq_ref = 'NA'
    else:
        # ... Fallback logic ...