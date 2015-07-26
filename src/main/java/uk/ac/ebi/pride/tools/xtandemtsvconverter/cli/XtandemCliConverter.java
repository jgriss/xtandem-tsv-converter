package uk.ac.ebi.pride.tools.xtandemtsvconverter.cli;

import uk.ac.ebi.pride.tools.xtandemtsvconverter.util.PSM;
import uk.ac.ebi.pride.tools.xtandemtsvconverter.util.PTM;
import uk.ac.ebi.pride.tools.xtandemtsvconverter.util.SimpleXTandemParser;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by jg on 16.07.15.
 */
public class XtandemCliConverter {
    public static void main(String[] args) {
        try {
            if (args.length != 2 && args.length != 3) {
                printUsage();
                return;
            }

            File inputFile = new File(args[0]);
            File outputFile = new File(args[0] + ".tsv");
            Double fdr = new Double(args[1]);

            File mgfFile = null;
            if (args.length >= 3)
                mgfFile = new File(args[2]);

            // make sure the input file and output file does not exist
            if (!inputFile.exists())
                throw new Exception("Cannot find input file: " + args[1]);

            if (outputFile.exists())
                throw new Exception("Outputfile '" + outputFile.toString() + "' already exists");

            // load the mgf titles is set
            List<String> titles = null;
            if (mgfFile != null) {
                titles = loadTitlesFromMgf(mgfFile);
            }

            // load the sequences
            List<String> sequences = null;
            if (mgfFile != null) {
                sequences = loadSequencesFromFile(mgfFile);
            }

            // convert the file
            convertFile(inputFile, outputFile, fdr, titles, sequences);

            System.out.println("Results written to " + outputFile);
        }
        catch (Exception e) {
            e.printStackTrace();
            System.out.println("\nError: " + e.getMessage());
            System.exit(1);
        }
    }

    private static List<String> loadSequencesFromFile(File mgfFile) throws Exception {
        List<String> sequences = loadPropertyFromMgf(mgfFile, "SEQ=");

        if (sequences.size() < 1)
            sequences = loadPropertyFromMgf(mgfFile, "SEQUENCE=");

        return sequences;
    }

    private static List<String> loadTitlesFromMgf(File mgfFile) throws Exception {
        return loadPropertyFromMgf(mgfFile, "TITLE=");
    }

    private static List<String> loadPropertyFromMgf(File mgfFile, String property) throws Exception {
        if (!mgfFile.exists())
            throw new Exception("Cannot find MGF file: " + mgfFile);

        BufferedReader reader = new BufferedReader(new FileReader(mgfFile));
        List<String> titles = new ArrayList<String>();
        String line;

        while ((line = reader.readLine()) != null) {
            if (!line.startsWith(property))
                continue;

            titles.add(line.substring(property.length()).trim());
        }

        reader.close();

        return titles;
    }

    protected static void convertFile(File inputFile, File outputFile, Double fdr, List<String> titles, List<String> sequences) throws Exception {
        List<PSM> psms = extractPsms(inputFile);
        Collections.sort(psms);

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));

        // write the header
        writer.write("spectrum_id\tprotein\tsequence\tptms\texpect\tis_decoy\tspectrum_sequence\tspectrum_title\n");

        int nTarget = 0, nDecoy = 0;

        for (PSM psm : psms) {
            if (psm.isDecoy())
                nDecoy++;
            else
                nTarget++;

            double currentFdr = (double) nDecoy / nTarget;

            if (currentFdr > fdr)
                break;

            // write the psm
            writer.write(String.valueOf(psm.getSpectrumIndex()) + "\t");
            writer.write(psm.getProteinAccession() + "\t");
            writer.write(psm.getSequence() + "\t");
            writer.write(buidPtmString(psm.getPtms()) + "\t");
            writer.write(String.valueOf(psm.getExpect()) + "\t");
            writer.write(String.valueOf(psm.isDecoy()) + "\t");

            String sequence = "";
            if (sequences != null) {
                if (psm.getSpectrumIndex() <= sequences.size()) {
                    sequence = sequences.get(psm.getSpectrumIndex() - 1);
                }
            }

            String title = "";
            if (titles != null) {
                if (psm.getSpectrumIndex() <= titles.size()) {
                    title = titles.get(psm.getSpectrumIndex() - 1);
                }
            }

            writer.write(sequence + "\t");
            writer.write(title + "\n");
        }

        writer.close();
    }

    private static String buidPtmString(List<PTM> ptms) {
        StringBuilder ptmString = new StringBuilder();

        for (int i = 0; i < ptms.size(); i++) {
            if (i > 0)
                ptmString.append(",");

            PTM ptm = ptms.get(i);
            ptmString.append(ptm.getDelta() + "@" + ptm.getPosition());
        }

        return ptmString.toString();
    }

    protected static List<PSM> extractPsms(File inputFile) throws Exception {
        SimpleXTandemParser simpleXTandemParser = new SimpleXTandemParser(inputFile);
        return new ArrayList<PSM>(simpleXTandemParser.getPsms());
    }

    private static void printUsage() {
        System.out.println("X!Tandem to TSV converter version 1.0");
        System.out.println("Usage: java -jar xtandem-tsv-converter.jar [input file] [FDR] [(optional) MGF file]");
    }
}
