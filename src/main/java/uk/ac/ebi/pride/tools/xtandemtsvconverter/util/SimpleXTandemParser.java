package uk.ac.ebi.pride.tools.xtandemtsvconverter.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.*;

/**
 * Created by jg on 18.07.15.
 */
public class SimpleXTandemParser {
    private String decoyTag = "###REV###";
    private List<PSM> psms = new ArrayList<PSM>();

    public SimpleXTandemParser(File resultFile) throws Exception {
        BufferedReader reader = new BufferedReader(new FileReader(resultFile));
        String line;

        // save the status of the current protein
        Map<Integer, List<PSM>> psmsPerSpectrum = new HashMap<Integer, List<PSM>>();
        int currentGroupId = 0;
        String currentAccession = null;
        boolean isDecoy = false;

        while ((line = reader.readLine()) != null) {
            line = line.trim();

            if (line.startsWith("<protein")) {
                String accession = getFieldValue("label", line);

                isDecoy = accession.contains(decoyTag);
                currentAccession = accession;
            }
            else if (line.startsWith("<group")) {
                if (!line.contains("type=\"model\""))
                    continue;
                currentGroupId = new Integer(getFieldValue("id", line));
            }
            else if (line.startsWith("<domain")) {
                PSM psm = new PSM(
                        currentGroupId,
                        getFieldValue("seq", line),
                        currentAccession,
                        isDecoy,
                        new Double(getFieldValue("expect", line)));

                if (!psmsPerSpectrum.containsKey(currentGroupId))
                    psmsPerSpectrum.put(currentGroupId, new ArrayList<PSM>());
                psmsPerSpectrum.get(currentGroupId).add(psm);
            }
        }

        reader.close();

        // only keep one PSM per spectrum
        for (int specId : psmsPerSpectrum.keySet()) {
            List<PSM> specPsms = psmsPerSpectrum.get(specId);
            Collections.sort(specPsms);
            psms.add(specPsms.get(0)); // lowest expect is best hit
        }
    }

    private String getFieldValue(String field, String xmlLine) {
        int pos = xmlLine.indexOf(field + "=");

        if (pos < 0)
            return null;

        int endIndex = xmlLine.indexOf("\"", pos + field.length() + 2);

        return xmlLine.substring(pos + field.length() + 2, endIndex);
    }

    public List<PSM> getPsms() {
        return Collections.unmodifiableList(psms);
    }
}
