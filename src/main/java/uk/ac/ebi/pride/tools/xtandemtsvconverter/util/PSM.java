package uk.ac.ebi.pride.tools.xtandemtsvconverter.util;

/**
 * Created by jg on 16.07.15.
 */
public class PSM implements Comparable<PSM>{
    private final String sequence;
    private final String proteinAccession;
    private final boolean isDecoy;
    private final double expect;
    private final int spectrumIndex;

    public PSM(int spectrumIndex, String sequence, String proteinAccession, boolean isDecoy, double expect) {
        this.spectrumIndex = spectrumIndex;
        this.sequence = sequence;
        this.proteinAccession = proteinAccession;
        this.isDecoy = isDecoy;
        this.expect = expect;
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

    @Override
    public int compareTo(PSM o) {
        return Double.compare(expect, o.getExpect());
    }

    @Override
    public String toString() {
        return String.format("[Spectrum = %d, expect = %f, isDecoy = %s, sequence = %s, protein = %s]", spectrumIndex, expect, String.valueOf(isDecoy), sequence, proteinAccession);
    }
}
