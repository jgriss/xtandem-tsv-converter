package uk.ac.ebi.pride.tools.xtandemtsvconverter.util;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Created by jg on 16.07.15.
 */
public class PSM implements Comparable<PSM>{
    private final String sequence;
    private final String proteinAccession;
    private final boolean isDecoy;
    private final double expect;
    private final int spectrumIndex;
    private final List<PTM> ptms;

    public PSM(int spectrumIndex, String sequence, String proteinAccession, boolean isDecoy, double expect, List<PTM> ptms) {
        this.spectrumIndex = spectrumIndex;
        this.sequence = sequence;
        this.proteinAccession = proteinAccession;
        this.isDecoy = isDecoy;
        this.expect = expect;
        this.ptms = ptms;
    }

    public int getSpectrumIndex() {
        return spectrumIndex;
    }

    public String getSequence() {
        return sequence;
    }

    public String getProteinAccession() {
        return proteinAccession;
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public double getExpect() {
        return expect;
    }

    public List<PTM> getPtms() {
        return Collections.unmodifiableList(ptms);
    }

    @Override
    public int compareTo(PSM o) {
        int result = Double.compare(expect, o.getExpect());

        if (result != 0)
            return result;

        return String.CASE_INSENSITIVE_ORDER.compare(sequence, o.getSequence());
    }

    @Override
    public String toString() {
        return String.format("[Spectrum = %d, expect = %f, isDecoy = %s, sequence = %s, protein = %s]", spectrumIndex, expect, String.valueOf(isDecoy), sequence, proteinAccession);
    }
}
