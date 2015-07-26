package uk.ac.ebi.pride.tools.xtandemtsvconverter.util;

/**
 * Created by jg on 24.07.15.
 */
public class PTM implements Comparable<PTM> {
    private final String aa;
    private final int position;
    private final double delta;

    public PTM(String aa, int position, double delta) {
        this.aa = aa;
        this.position = position;
        this.delta = delta;
    }

    public String getAa() {
        return aa;
    }

    public int getPosition() {
        return position;
    }

    public double getDelta() {
        return delta;
    }

    @Override
    public String toString() {
        return "[" + String.valueOf(delta) + "@" + String.valueOf(position) + "]";
    }

    @Override
    public int compareTo(PTM o) {
        int compare = Integer.compare(this.getPosition(), o.getPosition());

        if (compare != 0)
            return compare;

        return Double.compare(this.getDelta(), o.getDelta());
    }
}
