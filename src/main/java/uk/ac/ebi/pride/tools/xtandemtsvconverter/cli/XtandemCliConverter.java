package uk.ac.ebi.pride.tools.xtandemtsvconverter.cli;

import de.proteinms.xtandemparser.xtandem.*;
import uk.ac.ebi.pride.tools.xtandemtsvconverter.util.PSM;

import java.io.*;
import java.util.*;

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

            // convert the file
            convertFile(inputFile, outputFile, fdr, titles);

            System.out.println("Results written to " + outputFile);
        }
        catch (Exception e) {
            e.printStackTrace();
            System.out.println("\nError: " + e.getMessage());
            System.exit(1);
        }
    }

    private static List<String> loadTitlesFromMgf(File mgfFile) throws Exception {
        if (!mgfFile.exists())
            throw new Exception("Cannot find MGF file: " + mgfFile);

        BufferedReader reader = new BufferedReader(new FileReader(mgfFile));
        List<String> titles = new ArrayList<String>();
        String line;

        while ((line = reader.readLine()) != null) {
            if (!line.startsWith("TITLE="))
                continue;

            titles.add(line.substring("TITLE=".length()).trim());
        }

        reader.close();

        return titles;
    }

    protected static void convertFile(File inputFile, File outputFile, Double fdr, List<String> titles) throws Exception {
        List<PSM> psms = extractPsms(inputFile);
        Collections.sort(psms);

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));

        // write the header
        writer.write("spectrum_id\tprotein\tsequence\texpect\tis_decoy\tspectrum_title\n");

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
            writer.write(String.valueOf(psm.getExpect()) + "\t");
            writer.write(String.valueOf(psm.isDecoy()) + "\t");

            String title = "";
            if (titles != null) {
                if (psm.getSpectrumIndex() <= titles.size()) {
                    title = titles.get(psm.getSpectrumIndex() - 1);
                }
            }

            writer.write(title + "\n");
        }

        writer.close();
    }

    protected static List<PSM> extractPsms(File inputFile) throws Exception {
        XTandemFile xtandemFile = new XTandemFile(inputFile.getAbsolutePath());
        PeptideMap peptideMap = xtandemFile.getPeptideMap();
        ProteinMap proteinMap = xtandemFile.getProteinMap();
        List<PSM> psms = new ArrayList<PSM>();

        if (xtandemFile.getSpectraNumber() < 1)
            throw new Exception("X!Tandem result file does not contain any results");

        // loop through all spectra
        int nReliableTarget = 0, nReliableDecoy = 0;
        Map<String, Integer> sequenceCounts = new HashMap<String, Integer>();

        for (int spectrumIndex = 1; spectrumIndex <= xtandemFile.getSpectraNumber(); spectrumIndex++) {
            // all peptides for one spectrum
            List<Peptide> peptides = peptideMap.getAllPeptides(spectrumIndex);

            // rank-1 identification
            Domain bestDomain = getBestDomainForPeptides(peptides);

            // test whether the identification is a decoy ID
            Protein protein = proteinMap.getProtein(bestDomain.getProteinKey());

            boolean isDecoy = protein.getDescription().contains("###REV###") || protein.getDescription().contains("###RND###");

            // extract the spectrum index
            Integer fileSpectrumIndex = new Integer(bestDomain.getDomainID().substring(0, bestDomain.getDomainID().indexOf('.')));

            PSM psm = new PSM(fileSpectrumIndex, bestDomain.getDomainSequence(), protein.getDescription(), isDecoy, bestDomain.getDomainExpect());
            psms.add(psm);
        }

        return psms;
    }

    private static void printUsage() {
        System.out.println("X!Tandem to TSV converter version 1.0");
        System.out.println("Usage: java -jar xtandem-tsv-converter.jar [input file] [FDR] [(optional) MGF file]");
    }

    private static Domain getBestDomainForPeptides(List<Peptide> peptides) {
        double bestExpectScore = Double.MAX_VALUE;
        Domain bestDomain = null;

        for (Peptide peptide : peptides) {
            for (Domain domain : peptide.getDomains()) {
                if (domain.getDomainExpect() < bestExpectScore) {
                    bestExpectScore = domain.getDomainExpect();
                    bestDomain = domain;
                }
            }


        }

        return bestDomain;
    }
}
